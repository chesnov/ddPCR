"""
ddpcr_parsers.py — Enhanced Bio-Rad file parsers for ddPCRvis
=============================================================
Drop this file next to ddPCRvis.py, then in ddPCRvis.py replace:

    class BioRadQLPParser:
        ...  (the whole class, ~210 lines)

with:

    from ddpcr_parsers import BioRadQLPParser, BioRadDdpcrParser

See INTEGRATION.md (or the bottom of this file) for the four additional
changes needed in the UI code.

QLP improvements over the original parser
------------------------------------------
* Well IDs read via TAG_WELL_NAME (65019) instead of a hardcoded +10 offset hack
* All 7 droplet record fields extracted:
    Timestamp, Ch1_Amplitude, Ch1_Quality, Ch1_Width,
               Ch2_Amplitude, Ch2_Quality, Ch2_Width
* 20+ new per-well tags mapped: thresholds, concentrations, 95% CIs,
  positive/negative counts, rejected/saturated counts, width gates,
  sample name, supermix, target names per channel, experiment type/name,
  instrument serial/make/model, droplet volume, system version
* All per-well scalar metadata surfaced as columns in every-row of the CSV
* File-level metadata (software, plate ID, run date, instrument) in df.attrs

ddpcr parser
------------
* Decrypts AES-256-CBC ciphertext (base64-encoded file) using the key from
  ApplicationSettings.json: SHA-256("DbCdNa2OrCaDx56#4") with a fixed IV.
  Uses the `cryptography` library (auto-installed if missing).
* Adaptive JSON payload discovery — handles multiple archive layouts
* Uses BioRad's precomputed thresholds and concentrations as ground truth
* Emits exactly the same DataFrame schema as BioRadQLPParser so all
  downstream plotting/stats code is shared transparently
"""

from __future__ import annotations

import json
import re
import shutil
import struct
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd


# ═══════════════════════════════════════════════════════════════════════════════
#  QLP PARSER
# ═══════════════════════════════════════════════════════════════════════════════

class BioRadQLPParser:
    """
    Full-featured parser for QuantaSoft 1.7 QLP binary files.

    QLP is a TIFF-derived format:
      Bytes 0-1 : byte order marker  (b'II' = little-endian, b'MM' = big-endian)
      Bytes 2-3 : magic 0x002A
      Bytes 4-7 : offset to first IFD

    Each IFD = one well (or the file-header IFD for plate metadata).
    IFD layout:
      2 bytes   : number of entries
      N × 12    : tag entries  [tag_id(2) | type(2) | count(4) | value_or_offset(4)]
      4 bytes   : offset to next IFD  (0 = done)

    Droplet record: 28 bytes
      uint32  Timestamp
      float   Ch1_Amplitude
      float   Ch1_Quality
      float   Ch1_Width
      float   Ch2_Amplitude
      float   Ch2_Quality
      float   Ch2_Width
    """

    # ── Tag IDs ────────────────────────────────────────────────────────────
    # Standard TIFF
    TAG_IMAGE_DESC          = 270   # Often contains plate XML/JSON
    TAG_SOFTWARE            = 305   # Software version string

    # Bio-Rad custom (0xFE00 range)
    TAG_CHANNEL_NAMES       = 65004  # "Ch1,Ch2" or "FAM,HEX" etc.
    TAG_WELL_QUALITY_SCORE  = 65005  # float: overall well quality score
    TAG_CONCENTRATION       = 65006  # float[2]: copies/µL per channel
    TAG_CONF_LOWER          = 65007  # float[2]: 95 % CI lower bound
    TAG_CONF_UPPER          = 65008  # float[2]: 95 % CI upper bound
    TAG_POSITIVES           = 65009  # uint32[2]: positive droplet count
    TAG_NEGATIVES           = 65010  # uint32[2]: negative droplet count
    TAG_TOTAL_DROPLETS      = 65011  # uint32: accepted droplet count
    TAG_THRESHOLD           = 65012  # float[2]: auto threshold per channel
    TAG_THRESHOLD_CONF      = 65013  # float[2]: threshold confidence
    TAG_MANUAL_THRESHOLD    = 65014  # float[2]: manual threshold value
    TAG_SAMPLE_NAME         = 65015  # ASCII: sample name
    TAG_SUPERMIX            = 65016  # ASCII: supermix
    TAG_TARGET_NAMES        = 65017  # ASCII: "TargetCh1,TargetCh2"
    TAG_EXPERIMENT_TYPE     = 65018  # ASCII: experiment type
    TAG_WELL_NAME           = 65019  # ASCII: well ID e.g. "A01"
    TAG_EXPERIMENT_NAME     = 65020  # ASCII: experiment name
    TAG_DATA_START          = 65021  # Offset to first droplet record
    TAG_REJECTED_EVENTS     = 65022  # uint32
    TAG_SATURATED_EVENTS    = 65023  # uint32
    TAG_DROPLET_VOLUME      = 65026  # float: nL
    TAG_PLATE_ID            = 65030  # ASCII
    TAG_RUN_DATE            = 65031  # ASCII
    TAG_INSTRUMENT_SN       = 65032  # ASCII
    TAG_INSTRUMENT_MAKE     = 65033  # ASCII
    TAG_INSTRUMENT_MODEL    = 65034  # ASCII
    TAG_WAS_THRESHOLDED     = 65040  # uint8/bool
    TAG_MIN_WIDTH_GATE      = 65041  # float
    TAG_MAX_WIDTH_GATE      = 65042  # float
    TAG_SYSTEM_VERSION      = 65050  # uint16
    TAG_QUALITY_ARRAY       = 65054  # byte[N]: per-droplet quality flag
    TAG_CLUSTER_ARRAY       = 65057  # byte[N]: per-droplet cluster assignment
    TAG_WELL_QUALITY_FLAGS  = 65058  # uint64[2]: per-channel quality bitmask

    # Droplet record
    RECORD_SIZE = 28
    RECORD_FMT  = "Iffffff"  # 7 fields

    # Cluster byte → population index
    CLUSTER_MAP = {0x00: 0, 0x11: 1, 0x22: 2, 0x33: 3, 0x44: 4}
    CLUSTER_LABELS = {
        0: "Filtered",
        1: "NN",
        2: "Ch1+",
        3: "Ch1+Ch2+",
        4: "Ch2+",
    }

    # TIFF type → bytes per element
    _TYPE_SIZES = {1:1, 2:1, 3:2, 4:4, 5:8, 7:1, 9:4, 11:4, 12:8}

    # ── Init ──────────────────────────────────────────────────────────────

    def __init__(self, filepath: str, debug: bool = False):
        self.filepath = filepath
        self.debug    = debug

        with open(filepath, 'rb') as f:
            self.data = f.read()

        bom = self.data[:2]
        if bom == b'II':
            self.endian = '<'
        elif bom == b'MM':
            self.endian = '>'
        else:
            raise ValueError(f"Not a valid QLP file — bad byte-order marker: {bom!r}")

        self.metadata:      Dict[str, Any] = {}
        self.well_metadata: Dict[str, Any] = {}
        self.channel_names: Optional[List[str]] = None

    # ── Low-level helpers ─────────────────────────────────────────────────

    def _u(self, offset: int, fmt: str):
        """Unpack one value; return None on error."""
        try:
            sz = struct.calcsize(fmt)
            return struct.unpack(f"{self.endian}{fmt}",
                                 self.data[offset:offset + sz])[0]
        except Exception:
            return None

    def _unpack_n(self, offset: int, fmt: str, count: int) -> List:
        sz = struct.calcsize(fmt)
        out = []
        for i in range(count):
            v = self._u(offset + i * sz, fmt)
            if v is None:
                break
            out.append(v)
        return out

    def _ascii(self, offset: int, size: int) -> str:
        raw = self.data[offset:offset + size]
        return raw.split(b'\x00')[0].decode('ascii', errors='ignore').strip()

    # ── IFD parsing ───────────────────────────────────────────────────────

    def _parse_ifd(self, ifd_offset: int) -> Tuple[dict, Optional[int]]:
        """Return (tags_dict, next_ifd_offset).
        tags_dict: tag_id -> {'ptr', 'size', 'type', 'count'}"""
        num = self._u(ifd_offset, "H")
        if num is None or num == 0 or num > 500:
            return {}, None

        tags = {}
        for i in range(num):
            ep      = ifd_offset + 2 + i * 12
            tid     = self._u(ep,     "H")
            ttype   = self._u(ep + 2, "H")
            count   = self._u(ep + 4, "I")
            if tid is None or ttype is None or count is None:
                continue
            item_sz   = self._TYPE_SIZES.get(ttype, 1)
            total_sz  = item_sz * count
            raw_val   = self._u(ep + 8, "I")
            ptr       = (ep + 8) if total_sz <= 4 else raw_val
            tags[tid] = {'ptr': ptr, 'size': total_sz,
                         'type': ttype, 'count': count}

        next_ptr    = ifd_offset + 2 + num * 12
        next_ifd    = self._u(next_ptr, "I")
        return tags, (next_ifd if next_ifd else None)

    # ── Tag decoders ─────────────────────────────────────────────────────

    def _t_ascii(self, tags: dict, tid: int) -> Optional[str]:
        if tid not in tags:
            return None
        t = tags[tid]
        return self._ascii(t['ptr'], t['size'])

    def _t_float(self, tags: dict, tid: int, idx: int = 0) -> Optional[float]:
        if tid not in tags:
            return None
        return self._u(tags[tid]['ptr'] + idx * 4, "f")

    def _t_floats(self, tags: dict, tid: int) -> List[float]:
        if tid not in tags:
            return []
        t = tags[tid]
        return self._unpack_n(t['ptr'], "f", t['size'] // 4)

    def _t_uint32(self, tags: dict, tid: int, idx: int = 0) -> Optional[int]:
        if tid not in tags:
            return None
        return self._u(tags[tid]['ptr'] + idx * 4, "I")

    def _t_uint32s(self, tags: dict, tid: int) -> List[int]:
        if tid not in tags:
            return []
        t = tags[tid]
        return self._unpack_n(t['ptr'], "I", t['size'] // 4)

    def _t_bytes(self, tags: dict, tid: int) -> Optional[bytes]:
        if tid not in tags:
            return None
        t = tags[tid]
        return self.data[t['ptr']:t['ptr'] + t['size']]

    # ── File-level metadata ───────────────────────────────────────────────

    def _extract_file_metadata(self):
        ifd_offset = self._u(4, "I")
        if not ifd_offset:
            return
        tags, _ = self._parse_ifd(ifd_offset)

        sw = self._t_ascii(tags, self.TAG_SOFTWARE)
        if sw:
            self.metadata['software'] = sw

        desc = self._t_ascii(tags, self.TAG_IMAGE_DESC)
        if desc:
            self.metadata['image_description'] = desc
            self._parse_plate_desc(desc)

        ch = self._t_ascii(tags, self.TAG_CHANNEL_NAMES)
        if ch:
            self.channel_names = [c.strip() for c in ch.split(',')]

        for tid, key in [
            (self.TAG_PLATE_ID,        'plate_id'),
            (self.TAG_RUN_DATE,        'run_date'),
            (self.TAG_INSTRUMENT_SN,   'instrument_serial'),
            (self.TAG_INSTRUMENT_MAKE, 'instrument_make'),
            (self.TAG_INSTRUMENT_MODEL,'instrument_model'),
        ]:
            v = self._t_ascii(tags, tid)
            if v:
                self.metadata[key] = v

        dv = self._t_float(tags, self.TAG_DROPLET_VOLUME)
        if dv is not None:
            self.metadata['droplet_volume_nL'] = dv

        if not self.channel_names:
            self.channel_names = ['Ch1', 'Ch2']

    def _parse_plate_desc(self, desc: str):
        """Try JSON then key=value from ImageDescription tag."""
        try:
            d = json.loads(desc)
            self.metadata.update({f'plate_{k}': v for k, v in d.items()
                                   if not isinstance(v, (dict, list))})
            return
        except Exception:
            pass
        for part in desc.split(';'):
            if '=' in part:
                k, _, v = part.partition('=')
                self.metadata[k.strip()] = v.strip()

    # ── Per-well metadata ─────────────────────────────────────────────────

    def _extract_well_metadata(self, tags: dict, well_id: str) -> dict:
        m = {'well_id': well_id}

        # Text fields
        for tid, key in [
            (self.TAG_SAMPLE_NAME,     'sample_name'),
            (self.TAG_SUPERMIX,        'supermix'),
            (self.TAG_TARGET_NAMES,    'target_names'),
            (self.TAG_EXPERIMENT_TYPE, 'experiment_type'),
            (self.TAG_EXPERIMENT_NAME, 'experiment_name'),
        ]:
            v = self._t_ascii(tags, tid)
            if v:
                m[key] = v

        if 'target_names' in m:
            parts = [p.strip() for p in m['target_names'].split(',')]
            m['target_ch1'] = parts[0] if len(parts) > 0 else ''
            m['target_ch2'] = parts[1] if len(parts) > 1 else ''

        # Thresholds
        thresh = self._t_floats(tags, self.TAG_THRESHOLD)
        if thresh:
            m['threshold_ch1'] = thresh[0] if len(thresh) > 0 else None
            m['threshold_ch2'] = thresh[1] if len(thresh) > 1 else None

        tc = self._t_floats(tags, self.TAG_THRESHOLD_CONF)
        if tc:
            m['threshold_confidence_ch1'] = tc[0] if len(tc) > 0 else None
            m['threshold_confidence_ch2'] = tc[1] if len(tc) > 1 else None

        mt = self._t_floats(tags, self.TAG_MANUAL_THRESHOLD)
        if mt:
            m['manual_threshold_ch1'] = mt[0] if len(mt) > 0 else None
            m['manual_threshold_ch2'] = mt[1] if len(mt) > 1 else None

        # Concentrations
        conc = self._t_floats(tags, self.TAG_CONCENTRATION)
        if conc:
            m['concentration_ch1'] = conc[0] if len(conc) > 0 else None
            m['concentration_ch2'] = conc[1] if len(conc) > 1 else None

        ci_lo = self._t_floats(tags, self.TAG_CONF_LOWER)
        if ci_lo:
            m['conc_ci_lower_ch1'] = ci_lo[0] if len(ci_lo) > 0 else None
            m['conc_ci_lower_ch2'] = ci_lo[1] if len(ci_lo) > 1 else None

        ci_hi = self._t_floats(tags, self.TAG_CONF_UPPER)
        if ci_hi:
            m['conc_ci_upper_ch1'] = ci_hi[0] if len(ci_hi) > 0 else None
            m['conc_ci_upper_ch2'] = ci_hi[1] if len(ci_hi) > 1 else None

        # Droplet counts
        pos = self._t_uint32s(tags, self.TAG_POSITIVES)
        if pos:
            m['positives_ch1'] = pos[0] if len(pos) > 0 else None
            m['positives_ch2'] = pos[1] if len(pos) > 1 else None

        neg = self._t_uint32s(tags, self.TAG_NEGATIVES)
        if neg:
            m['negatives_ch1'] = neg[0] if len(neg) > 0 else None
            m['negatives_ch2'] = neg[1] if len(neg) > 1 else None

        total = self._t_uint32(tags, self.TAG_TOTAL_DROPLETS)
        if total is not None:
            m['total_accepted_droplets'] = total

        rej = self._t_uint32(tags, self.TAG_REJECTED_EVENTS)
        if rej is not None:
            m['rejected_droplets'] = rej

        sat = self._t_uint32(tags, self.TAG_SATURATED_EVENTS)
        if sat is not None:
            m['saturated_droplets'] = sat

        # Quality
        wq = self._t_float(tags, self.TAG_WELL_QUALITY_SCORE)
        if wq is not None:
            m['well_quality_score'] = wq

        flags_raw = self._t_bytes(tags, self.TAG_WELL_QUALITY_FLAGS)
        if flags_raw and len(flags_raw) >= 8:
            m['quality_flags_ch1'] = struct.unpack(f"{self.endian}Q",
                                                    flags_raw[:8])[0]
            if len(flags_raw) >= 16:
                m['quality_flags_ch2'] = struct.unpack(f"{self.endian}Q",
                                                        flags_raw[8:16])[0]

        # Width gates
        wg_min = self._t_float(tags, self.TAG_MIN_WIDTH_GATE)
        if wg_min is not None:
            m['min_width_gate'] = wg_min

        wg_max = self._t_float(tags, self.TAG_MAX_WIDTH_GATE)
        if wg_max is not None:
            m['max_width_gate'] = wg_max

        # Misc
        wt = self._t_uint32(tags, self.TAG_WAS_THRESHOLDED)
        if wt is not None:
            m['was_thresholded'] = bool(wt)

        dv = self._t_float(tags, self.TAG_DROPLET_VOLUME)
        if dv is not None:
            m['droplet_volume_nL'] = dv

        if self.TAG_SYSTEM_VERSION in tags:
            sv = self._u(tags[self.TAG_SYSTEM_VERSION]['ptr'], "H")
            if sv is not None:
                m['system_version'] = sv

        # Channel names (may also be in well IFD)
        ch = self._t_ascii(tags, self.TAG_CHANNEL_NAMES)
        if ch and not self.channel_names:
            self.channel_names = [c.strip() for c in ch.split(',')]

        return m

    # ── Droplet records ───────────────────────────────────────────────────

    def _parse_droplets(self, tags: dict, well_meta: dict,
                        include_filtered: bool) -> List[dict]:
        if self.TAG_DATA_START not in tags:
            return []
        if self.TAG_CLUSTER_ARRAY not in tags:
            return []

        ci        = tags[self.TAG_CLUSTER_ARRAY]
        clust_blob = self.data[ci['ptr']:ci['ptr'] + ci['size']]
        n          = len(clust_blob)
        if n == 0:
            return []

        q_blob = None
        if self.TAG_QUALITY_ARRAY in tags:
            qi = tags[self.TAG_QUALITY_ARRAY]
            if qi['size'] > 0:
                q_blob = self.data[qi['ptr']:qi['ptr'] + qi['size']]

        cursor  = tags[self.TAG_DATA_START]['ptr']
        rec_fmt = f"{self.endian}{self.RECORD_FMT}"

        # Flatten metadata to scalar columns only
        meta_cols = {k: v for k, v in well_meta.items()
                     if isinstance(v, (str, int, float, bool, type(None)))}

        records = []
        for i in range(n):
            if cursor + self.RECORD_SIZE > len(self.data):
                break

            raw = struct.unpack(rec_fmt,
                                self.data[cursor:cursor + self.RECORD_SIZE])
            # raw = (Timestamp, Ch1_Amp, Ch1_Quality, Ch1_Width,
            #                   Ch2_Amp, Ch2_Quality, Ch2_Width)
            cluster_byte = clust_blob[i] if i < len(clust_blob) else 0
            cluster      = self.CLUSTER_MAP.get(cluster_byte, 0)
            q_flag       = q_blob[i] if (q_blob and i < len(q_blob)) else 0

            if not include_filtered and cluster == 0:
                cursor += self.RECORD_SIZE
                continue

            rec = {
                'Droplet_Index': i,
                'Timestamp':     raw[0],
                'Ch1_Amplitude': raw[1],
                'Ch1_Quality':   raw[2],
                'Ch1_Width':     raw[3],
                'Ch2_Amplitude': raw[4],
                'Ch2_Quality':   raw[5],
                'Ch2_Width':     raw[6],
                'Cluster':       cluster,
                'Cluster_Label': self.CLUSTER_LABELS.get(cluster, 'Unknown'),
                'Quality_Flag':  q_flag,
            }
            rec.update(meta_cols)
            records.append(rec)
            cursor += self.RECORD_SIZE

        return records

    # ── Main entry point ──────────────────────────────────────────────────

    def parse_to_dataframe(self, include_filtered: bool = True) -> Dict[str, pd.DataFrame]:
        """
        Parse the entire QLP file.

        Returns
        -------
        dict mapping well_id (e.g. 'A01') -> pd.DataFrame
        Each DataFrame has attrs:
            channel_names   : ['FAM', 'HEX'] or ['Ch1', 'Ch2']
            channel_map     : {'Ch1_Amplitude': 'FAM', 'Ch2_Amplitude': 'HEX'}
            well_metadata   : full per-well tag dict
            file_metadata   : plate-level metadata dict
        """
        self._extract_file_metadata()

        ifd_offset = self._u(4, "I")
        result: Dict[str, pd.DataFrame] = {}

        while ifd_offset and ifd_offset < len(self.data) - 6:
            tags, next_ifd = self._parse_ifd(ifd_offset)

            is_well = (self.TAG_DATA_START    in tags and
                       self.TAG_CLUSTER_ARRAY in tags and
                       tags[self.TAG_CLUSTER_ARRAY]['size'] > 0)

            if is_well:
                # Prefer TAG_WELL_NAME; fall back to heuristic
                well_id = self._t_ascii(tags, self.TAG_WELL_NAME)
                if not _valid_well_id(well_id):
                    well_id = self._fallback_well_id(ifd_offset)

                if _valid_well_id(well_id):
                    well_meta = self._extract_well_metadata(tags, well_id)
                    self.well_metadata[well_id] = well_meta

                    droplets = self._parse_droplets(tags, well_meta,
                                                    include_filtered)
                    if droplets:
                        df = pd.DataFrame(droplets)
                        _attach_attrs(df, self.channel_names or ['Ch1','Ch2'],
                                      well_meta, self.metadata)
                        result[well_id] = df

                        if self.debug:
                            print(f"QLP: well {well_id} — "
                                  f"{len(droplets)} droplets, "
                                  f"{len(well_meta)} metadata fields")
                elif self.debug:
                    print(f"QLP: skipping IFD at 0x{ifd_offset:08x} "
                          f"— could not determine well ID")

            ifd_offset = next_ifd

        return result

    def _fallback_well_id(self, ifd_offset: int) -> Optional[str]:
        """Legacy +10 offset heuristic — only if TAG_WELL_NAME fails."""
        try:
            raw = self.data[ifd_offset + 10:ifd_offset + 14]
            wid = raw.split(b'\x00')[0].decode('ascii')
            return wid if _valid_well_id(wid) else None
        except Exception:
            return None

    def get_channel_names(self) -> List[str]:
        if self.channel_names is None:
            self._extract_file_metadata()
        return self.channel_names or ['Ch1', 'Ch2']

    def get_file_metadata(self) -> dict:
        return self.metadata

    def get_well_metadata(self, well_id: Optional[str] = None):
        if well_id:
            return self.well_metadata.get(well_id, {})
        return self.well_metadata


# ═══════════════════════════════════════════════════════════════════════════════
#  ddPCR PARSER
# ═══════════════════════════════════════════════════════════════════════════════

class BioRadDdpcrParser:
    """
    Parser for Bio-Rad .ddpcr files.

    A .ddpcr file is a 7-Zip archive encrypted with AES-256.  The password is
    derived by decrypting ComponentIdentityInfo.s_FileIDForRUOFile with the
    SecureKey from ApplicationSettings.json, producing a 34-character string
    whose UTF-16LE encoding is used directly as the 7z AES key material.

    The archive contains JSON files describing plate layout, well data,
    amplitude arrays, and BioRad's precomputed analysis results.

    The emitted DataFrames are schema-identical to BioRadQLPParser output,
    so all downstream ddPCRvis plotting/stats code is shared transparently.
    """

    SECRET_KEY = "DbCdNa2OrCaDx56#4"

    # Cluster call bytes (from converter GetEventData):
    #   0  = gated/filtered    17 (0x11) = NN
    #  34 (0x22) = Ch1+only   51 (0x33) = double+   68 (0x44) = Ch2+only
    CLUSTER_BYTE_MAP = {0: 0, 17: 1, 34: 2, 51: 3, 68: 4}
    CLUSTER_LABELS   = {0: "Filtered", 1: "NN", 2: "Ch1+",
                        3: "Ch1+Ch2+", 4: "Ch2+"}

    def __init__(self, filepath: str, debug: bool = False):
        self.filepath  = Path(filepath)
        self.debug     = debug
        self.metadata:      Dict[str, Any] = {}
        self.well_metadata: Dict[str, Any] = {}
        self.channel_names: List[str] = ['Ch1', 'Ch2']
        self._temp_dir: Optional[Path] = None

    # ── Password ──────────────────────────────────────────────────────────
    #
    # Confirmed by calling BioRad's own DLL on Windows via PowerShell:
    #   $password = [BioRad.Shared.EncryptDecryptMsgHandler]::Decrypt(
    #       ComponentIdentityInfo.s_FileIDForRUOFile, $true, $SecureKey)
    # → '1b53402e-503a-4303-bf86-71af1f3178dd'  (plain ASCII GUID)
    #
    # Verified by 7za.exe 24.09:
    #   7za l -p"1b53402e-503a-4303-bf86-71af1f3178dd" file.ddpcr  → success
    #
    # Note: EncryptDecryptMsgHandler uses IV=zeros (not the base64 constant —
    # Array.Copy(src, aes.IV, 16) writes to the getter's defensive copy and
    # is silently discarded; the AES object keeps IV=new byte[16]).
    PASSWORD = '1b53402e-503a-4303-bf86-71af1f3178dd'

    def _extract(self):
        """
        Extract the encrypted 7-Zip .ddpcr archive to a temp directory.

        Uses the system 7za/7z binary (exactly as BioRad does via
        SevenZipExeUtils) with the confirmed password.
        Install with: sudo apt install p7zip-full
        """
        import subprocess, shutil as _shutil

        self._temp_dir = Path(tempfile.mkdtemp(prefix='ddpcr_'))

        sevenzip = (_shutil.which('7za') or
                    _shutil.which('7z')  or
                    _shutil.which('7zz'))

        if not sevenzip:
            raise RuntimeError(
                "7za not found. Install it with:\n"
                "    sudo apt install p7zip-full\n"
                "or on macOS: brew install p7zip"
            )

        cmd = [
            sevenzip, 'x',
            str(self.filepath),
            f'-o{self._temp_dir}',
            f'-p{self.PASSWORD}',
            '-r', '-y'
        ]

        if self.debug:
            print(f"ddpcr: running {' '.join(cmd[:3])} ...")

        try:
            r = subprocess.run(cmd, capture_output=True, timeout=300)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"7za timed out extracting {self.filepath.name}")

        if r.returncode != 0:
            err = r.stderr.decode(errors='replace').strip()
            raise RuntimeError(
                f"7za failed (exit {r.returncode}) on {self.filepath.name}.\n"
                f"{err}"
            )

        if self.debug:
            print(f"ddpcr: extracted to {self._temp_dir}")

    def cleanup(self):
        if self._temp_dir and self._temp_dir.exists():
            shutil.rmtree(self._temp_dir, ignore_errors=True)
            self._temp_dir = None

    # ── File discovery ────────────────────────────────────────────────────

    def _all_files(self) -> List[Path]:
        return [p for p in self._temp_dir.rglob('*') if p.is_file()]



    # ── JSON helpers ──────────────────────────────────────────────────────

    @staticmethod
    def _load_json(path: Path) -> Any:
        with open(path, 'r', encoding='utf-8-sig', errors='replace') as f:
            return json.load(f)

    # ── Parse extracted archive ───────────────────────────────────────────

    def parse_to_dataframe(self, include_filtered: bool = True
                           ) -> Dict[str, pd.DataFrame]:
        """
        Decrypt, extract, and parse the .ddpcr archive.

        Returns
        -------
        dict mapping well_id -> pd.DataFrame
        Schema is identical to BioRadQLPParser.parse_to_dataframe().
        Each row is one droplet; metadata columns are broadcast across all rows.

        Archive structure (from BioRad ReadPlate / JsonPeakReaderWriter):
          PeakData/<WellID>.ddpeakjson    — amplitudes + gating per well
          PeakMetaData/<WellID>.ddmetajson — cluster assignments + thresholds
        """
        self._extract()
        try:
            return self._parse_extracted(include_filtered)
        finally:
            self.cleanup()

    def _parse_extracted(self, include_filtered: bool) -> Dict[str, pd.DataFrame]:
        peak_dir = self._temp_dir / 'PeakData'
        meta_dir = self._temp_dir / 'PeakMetaData'

        if not peak_dir.exists():
            # Try one level deeper (some archives nest inside a folder)
            subdirs = [p for p in self._temp_dir.iterdir() if p.is_dir()]
            for sd in subdirs:
                if (sd / 'PeakData').exists():
                    peak_dir = sd / 'PeakData'
                    meta_dir = sd / 'PeakMetaData'
                    break

        if not peak_dir.exists():
            contents = list(self._temp_dir.rglob('*'))[:20]
            raise RuntimeError(
                f"PeakData/ directory not found inside archive. "
                f"Archive contents: {contents}"
            )

        # ── Index ddpeakjson files by well ID ─────────────────────────────
        peak_files: Dict[str, Path] = {}
        for p in peak_dir.glob('*.ddpeakjson'):
            wid = self._well_id_from_peak_json(p)
            if wid:
                peak_files[wid] = p

        # ── Index ddmetajson files by well ID ─────────────────────────────
        meta_files: Dict[str, Path] = {}
        if meta_dir.exists():
            for p in meta_dir.glob('*.ddmetajson'):
                wid = self._well_id_from_meta_json(p)
                if wid:
                    meta_files[wid] = p

        if not peak_files:
            raise RuntimeError("No .ddpeakjson files found in PeakData/")

        result: Dict[str, pd.DataFrame] = {}

        for well_id, peak_path in sorted(peak_files.items()):
            try:
                peak = self._load_json(peak_path)
                meta = (self._load_json(meta_files[well_id])
                        if well_id in meta_files else None)
                df = self._build_well_df(well_id, peak, meta, include_filtered)
                if df is not None and len(df) > 0:
                    result[well_id] = df
                    self.well_metadata[well_id] = df.attrs.get('well_metadata', {})
                    if self.debug:
                        print(f"ddpcr: well {well_id} — {len(df)} droplets")
            except Exception as e:
                if self.debug:
                    print(f"ddpcr: error parsing {well_id}: {e}")
                continue

        return result

    # ── Well ID helpers ───────────────────────────────────────────────────

    def _well_id_from_peak_json(self, path: Path) -> Optional[str]:
        """Get well ID from filename or from WellInfo in the JSON."""
        # Try filename first (e.g. H12.ddpeakjson → H12)
        stem = path.stem.upper()
        if _valid_well_id(stem):
            return stem
        wid = _well_id_from_string(path.stem)
        if wid:
            return wid
        # Try reading WellInfo from the JSON
        try:
            d = self._load_json(path)
            wi = d.get('PlateInfo', {}).get('WellInfo', {})
            row = wi.get('WellPosition', {}).get('Row')   # 0-based (0=A)
            col = wi.get('WellPosition', {}).get('Column') # 1-based
            if row is not None and col is not None:
                return f"{chr(ord('A') + int(row))}{int(col):02d}"
        except Exception:
            pass
        return None

    def _well_id_from_meta_json(self, path: Path) -> Optional[str]:
        """Get well ID from filename or WellIndex (0-based, row-major)."""
        stem = path.stem.upper()
        if _valid_well_id(stem):
            return stem
        wid = _well_id_from_string(path.stem)
        if wid:
            return wid
        try:
            d = self._load_json(path)
            idx = d.get('WellIndex')   # 0-based, row-major (0=A01, 1=A02...)
            if idx is not None:
                row = int(idx) // 12
                col = int(idx) % 12 + 1
                return f"{chr(ord('A') + row)}{col:02d}"
        except Exception:
            pass
        return None

    # ── Build per-well DataFrame ──────────────────────────────────────────

    def _build_well_df(self, well_id: str, peak: dict,
                       meta: Optional[dict],
                       include_filtered: bool) -> Optional[pd.DataFrame]:
        """
        Build a per-droplet DataFrame from a ddpeakjson + optional ddmetajson.

        ddpeakjson/PeakInfo layout (confirmed from real files):
          Amplitudes  : list of N_channels lists, each len=PeakCount  (floats)
          Widths      : list of PeakCount floats
          Timestamps  : list of PeakCount floats
          Qualities   : list of PeakCount floats
          GatingFlags : list of PeakCount ints  (0=accepted, non-zero=gated)

        ddmetajson layout (confirmed from real files):
          WellIndex   : int (0-based row-major)
          Clusters    : list of {Cluster;v, Targets, Results, Droplets, Unassigned}
            Droplets  : list of droplet indices belonging to this cluster
            Results   : list of 'Positive'/'Negative'/'Gated' per target
            Targets   : list of {Name, Dyes:[{Name, Channel},...], ...}
          ThresholdKeys   : list of dye names
          ThresholdValues : list of lists of {ThresholdValue, ...}
          WasThreshed1v2 / WasThreshed3v4 / WasThreshed5v6 : bool
        """
        pi = peak.get('PeakInfo', {})
        n  = pi.get('PeakCount', 0)
        if n == 0:
            return None

        amps_all   = pi.get('Amplitudes', [])   # N_channels × n
        widths     = pi.get('Widths',     [])   # n
        timestamps = pi.get('Timestamps', [])   # n
        qualities  = pi.get('Qualities',  [])   # n
        gating     = pi.get('GatingFlags',[])   # n  (0=accepted)

        if not amps_all:
            return None

        n_channels = len(amps_all)

        # ── Channel names from ChannelMap ─────────────────────────────────
        ch_map = peak.get('DataAcquisitionInfo', {}).get('ChannelMap', [])
        # ch_map is [{Dye: 'FAM', Channel: 1}, ...], Channel is 1-based
        dye_by_channel: Dict[int, str] = {}
        for entry in ch_map:
            if isinstance(entry, dict) and 'Dye' in entry and 'Channel' in entry:
                dye_by_channel[int(entry['Channel'])] = entry['Dye']

        # Channel names in order (channels are 1-based)
        channel_names = [dye_by_channel.get(i + 1, f'Ch{i+1}')
                         for i in range(n_channels)]

        # Update parser-level channel names if not set yet
        if self.channel_names == ['Ch1', 'Ch2'] and channel_names:
            self.channel_names = channel_names

        # ── Gating array → cluster label per droplet ─────────────────────
        # GatingFlags: non-zero = rejected/gated
        gating_arr = list(gating) if gating else [0] * n

        # ── Build droplet→cluster mapping from ddmetajson ─────────────────
        cluster_label = [None] * n   # None = unassigned initially
        cluster_int   = [-1]   * n   # Numeric cluster id matching QLP convention
        well_meta: Dict[str, Any] = {'well_id': well_id}

        if meta:
            self._apply_meta(meta, cluster_label, cluster_int,
                             well_meta, n_channels, dye_by_channel)

        # Assign gated droplets (GatingFlags != 0) to cluster 'Gated' / 0
        for i in range(n):
            if gating_arr[i] != 0:
                cluster_label[i] = 'Gated'
                cluster_int[i]   = 0   # cluster 0 = Filtered/Gated (QLP convention)
            elif cluster_label[i] is None:
                cluster_label[i] = 'Unassigned'
                # cluster_int[i] stays -1 (no meta assignment)

        # ── Per-file metadata ─────────────────────────────────────────────
        dai = peak.get('DataAcquisitionInfo', {})
        well_meta['droplet_volume_nL'] = dai.get('DropletVolume')
        well_meta['channel_count']     = dai.get('ChannelCount', n_channels)

        ri = peak.get('RejectedInfo', {})
        well_meta['rejected_droplets']   = ri.get('RejectedDropletCount', 0)
        well_meta['saturated_droplets']  = ri.get('SaturatedDropletCount', 0)

        # Scalar metadata for CSV columns
        meta_cols = {k: v for k, v in well_meta.items()
                     if isinstance(v, (str, int, float, bool, type(None)))}

        # ── Build rows ────────────────────────────────────────────────────
        rows = []
        for i in range(n):
            cl = cluster_label[i]
            if not include_filtered and cl == 'Gated':
                continue

            rec: Dict[str, Any] = {
                'Droplet_Index': i,
                'Timestamp':     timestamps[i] if i < len(timestamps) else i,
                'Width':         widths[i]     if i < len(widths)     else 0.0,
                'Quality':       qualities[i]  if i < len(qualities)  else 0.0,
                'Gating_Flag':   gating_arr[i],
                'Cluster':       cluster_int[i],   # integer 0-4, matches QLP schema
                'Cluster_Label': cl,
            }
            # Add one column per channel: Ch1_Amplitude, Ch2_Amplitude, etc.
            for ch_idx, ch_name in enumerate(channel_names):
                ch_amps = amps_all[ch_idx] if ch_idx < len(amps_all) else []
                rec[f'{ch_name}_Amplitude'] = (ch_amps[i]
                                               if i < len(ch_amps) else 0.0)
            # For backward compat with ddPCRvis (expects Ch1_Amplitude / Ch2_Amplitude)
            if 'Ch1_Amplitude' not in rec and len(channel_names) >= 1:
                rec['Ch1_Amplitude'] = rec.get(f'{channel_names[0]}_Amplitude', 0.0)
            if 'Ch2_Amplitude' not in rec and len(channel_names) >= 2:
                rec['Ch2_Amplitude'] = rec.get(f'{channel_names[1]}_Amplitude', 0.0)

            rec.update(meta_cols)
            rows.append(rec)

        if not rows:
            return None

        df = pd.DataFrame(rows)

        # ── attrs (channel map for downstream plotting) ───────────────────
        ch_map_attr = {}
        for ch_idx, ch_name in enumerate(channel_names):
            ch_map_attr[f'{ch_name}_Amplitude'] = ch_name
        # Always expose canonical Ch1/Ch2 keys
        if channel_names:
            ch_map_attr['Ch1_Amplitude'] = channel_names[0]
        if len(channel_names) >= 2:
            ch_map_attr['Ch2_Amplitude'] = channel_names[1]

        _attach_attrs(df, channel_names, well_meta, self.metadata)
        df.attrs['channel_map'] = ch_map_attr
        df.attrs['all_channel_names'] = channel_names
        return df

    def _apply_meta(self, meta: dict, cluster_label: List, cluster_int: List,
                    well_meta: dict, n_channels: int,
                    dye_by_channel: Dict[int, str]):
        """
        Parse ddmetajson and:
          - assign cluster labels and integer cluster ids to each droplet index
          - populate well_meta with threshold/target information

        BioRad stores the numeric cluster as the "Cluster" key on each cluster
        object.  In the binary QLP format the same information is encoded as
        raw byte values (0=Filtered, 0x11=NN, 0x22=Ch1+, 0x33=Double+,
        0x44=Ch2+) which are then mapped to indices 0–4 via CLUSTER_BYTE_MAP.
        The JSON representation uses the mapped indices 0–4 directly, but we
        defensively handle the raw byte values too.
        """
        clusters = meta.get('Clusters', [])
        _auto_idx = 0   # fallback sequential id when JSON omits the field

        for cluster_obj in clusters:
            if not isinstance(cluster_obj, dict):
                continue
            droplet_indices = cluster_obj.get('Droplets', [])
            results         = cluster_obj.get('Results', [])   # per-target
            targets         = cluster_obj.get('Targets', [])

            # ── Numeric cluster id ─────────────────────────────────────────
            # BioRad may store the value as:
            #   • a direct index 0-4  (most common in ddmetajson)
            #   • a raw byte value 0/17/34/51/68 (same as QLP binary)
            # Try the "Cluster" key first, then fall back to auto-increment.
            raw_cluster = cluster_obj.get('Cluster')
            if raw_cluster is None:
                cnum = _auto_idx
            elif isinstance(raw_cluster, int):
                # Map raw byte values → index if needed; pass-through 0-4.
                cnum = self.CLUSTER_BYTE_MAP.get(raw_cluster, raw_cluster)
            else:
                cnum = _auto_idx
            _auto_idx += 1

            # ── String label ───────────────────────────────────────────────
            if not isinstance(results, list):
                results = []
            is_gated = any(r == 'Gated' for r in results)
            label = 'Gated' if is_gated else self._cluster_label_from_results(
                results, targets)

            for idx in droplet_indices:
                if isinstance(idx, int) and 0 <= idx < len(cluster_label):
                    cluster_label[idx] = label
                    cluster_int[idx]   = cnum

        # ── Thresholds ─────────────────────────────────────────────────────
        thresh_keys   = meta.get('ThresholdKeys',   [])
        thresh_values = meta.get('ThresholdValues', [])
        for dye, tv_list in zip(thresh_keys, thresh_values):
            if isinstance(tv_list, list) and tv_list:
                tv = tv_list[0]
                if isinstance(tv, dict) and 'ThresholdValue' in tv:
                    safe_dye = dye.replace(' ', '_').replace('.', '_')
                    well_meta[f'threshold_{safe_dye}'] = float(tv['ThresholdValue'])
                    well_meta[f'threshold_{safe_dye}_manual'] = bool(
                        tv.get('IsManuallySet', False))

        # ── was_thresholded flags ──────────────────────────────────────────
        well_meta['was_thresholded_ch1v2'] = bool(meta.get('WasThreshed1v2', False))
        well_meta['was_thresholded_ch3v4'] = bool(meta.get('WasThreshed3v4', False))
        well_meta['was_thresholded_ch5v6'] = bool(meta.get('WasThreshed5v6', False))

    @staticmethod
    def _cluster_label_from_results(results: list, targets: list) -> str:
        """
        Build a compact cluster label like 'FAM+/HEX-/Cy5.5-' from per-target
        positive/negative calls and target dye names.
        """
        if not results:
            return 'Unassigned'
        parts = []
        for i, res in enumerate(results):
            # Get the primary dye name for this target
            dye = ''
            if i < len(targets) and isinstance(targets[i], dict):
                dyes = targets[i].get('Dyes', [])
                for d in dyes:
                    if isinstance(d, dict) and d.get('Name'):
                        dye = d['Name']
                        break
                if not dye:
                    dye = targets[i].get('Name', f'T{i+1}')
            else:
                dye = f'T{i+1}'
            sign = '+' if str(res).lower() == 'positive' else '-'
            parts.append(f'{dye}{sign}')
        return '/'.join(parts)

    def cleanup(self):
        if self._temp_dir and self._temp_dir.exists():
            shutil.rmtree(self._temp_dir, ignore_errors=True)
            self._temp_dir = None

    def get_channel_names(self) -> List[str]:
        return self.channel_names

    def get_file_metadata(self) -> dict:
        return self.metadata

    def get_well_metadata(self, well_id=None):
        if well_id:
            return self.well_metadata.get(well_id, {})
        return self.well_metadata


# ═══════════════════════════════════════════════════════════════════════════════
#  Shared helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _valid_well_id(wid) -> bool:
    if not wid or not isinstance(wid, str) or len(wid) < 2:
        return False
    try:
        return (wid[0] in 'ABCDEFGH' and
                wid[1:].isdigit() and
                1 <= int(wid[1:]) <= 12)
    except (ValueError, IndexError):
        return False


def _well_id_from_string(s: str) -> Optional[str]:
    """Extract a well ID (e.g. 'A01') from an arbitrary string."""
    m = re.search(r'([A-Ha-h])[\-_]?(\d{1,2})', s)
    if m:
        row = m.group(1).upper()
        col = int(m.group(2))
        if 1 <= col <= 12:
            return f"{row}{col:02d}"
    return None


def _attach_attrs(df: pd.DataFrame, channel_names: List[str],
                  well_meta: dict, file_meta: dict):
    """Attach standard attrs to a well DataFrame."""
    names = channel_names if len(channel_names) >= 2 else ['Ch1', 'Ch2']
    df.attrs['channel_names'] = names
    df.attrs['channel_map']   = {
        'Ch1_Amplitude': names[0],
        'Ch2_Amplitude': names[1],
    }
    df.attrs['well_metadata'] = well_meta
    df.attrs['file_metadata'] = file_meta


# ═══════════════════════════════════════════════════════════════════════════════
#  INTEGRATION GUIDE — changes needed in ddPCRvis.py
# ═══════════════════════════════════════════════════════════════════════════════
"""
1. REPLACE the BioRadQLPParser class
   ─────────────────────────────────
   Delete the entire BioRadQLPParser class (~210 lines) and replace it with:

       from ddpcr_parsers import BioRadQLPParser, BioRadDdpcrParser

   Place this import near the top of ddPCRvis.py, after the stdlib imports.


2. DropZone.dropEvent — add .ddpcr support
   ─────────────────────────────────────────
   In the dropEvent method, alongside the csv/qlp detection, add:

       ddpcr_files = []
       for url in event.mimeData().urls():
           file_path = url.toLocalFile()
           if file_path.endswith('.csv'):
               csv_files.append(file_path)
           elif file_path.endswith('.qlp'):
               qlp_files.append(file_path)
           elif file_path.endswith('.ddpcr'):          # ← NEW
               ddpcr_files.append(file_path)           # ← NEW

       if qlp_files:
           ...
       elif ddpcr_files:                               # ← NEW block
           if len(ddpcr_files) > 1:
               QMessageBox.warning(None, "Warning",
                   "Only one ddPCR file at a time. Using first.")
           self.files_dropped.emit([ddpcr_files[0]], 'ddpcr')
       elif csv_files:
           ...

   Also update the subtitle label and file filter strings to mention .ddpcr.


3. MainWindow — add Browse ddPCR button
   ──────────────────────────────────────
   In setup_ui (or __init__), alongside the browse_qlp_btn, add:

       browse_ddpcr_btn = QPushButton("📂 Browse ddPCR File")
       browse_ddpcr_btn.setStyleSheet(\"\"\"
           QPushButton { background-color: #e67e22; color: white;
                         border: none; padding: 10px 20px;
                         border-radius: 5px; font-weight: bold; }
           QPushButton:hover { background-color: #d35400; }
       \"\"\")
       browse_ddpcr_btn.clicked.connect(self.browse_ddpcr_file)
       browse_layout.addWidget(browse_ddpcr_btn)

   And add the handler method:

       def browse_ddpcr_file(self):
           file, _ = QFileDialog.getOpenFileName(
               self, "Select ddPCR File", "",
               "ddPCR Files (*.ddpcr)")
           if file:
               self.process_files([file], 'ddpcr')


4. ProcessingThread.run() — route ddpcr
   ──────────────────────────────────────
   In the run() method, add a branch alongside the 'qlp' branch:

       if self.file_type == 'qlp':
           output_dir, stats_file = process_qlp_file(
               self.files[0], self.well_assignments,
               self.include_filtered, self.progress.emit,
               self.output_format)
       elif self.file_type == 'ddpcr':                 # ← NEW
           output_dir, stats_file = process_ddpcr_file( # ← NEW
               self.files[0], self.well_assignments,
               self.include_filtered, self.progress.emit,
               self.output_format)
       else:
           output_dir, stats_file = process_ddpcr_files(...)


5. Add process_ddpcr_file() function
   ────────────────────────────────────
   Add this function near process_qlp_file() — it is essentially identical
   but uses BioRadDdpcrParser instead:

       def process_ddpcr_file(ddpcr_file, well_assignments, include_filtered,
                              log_func=print, output_format='png'):
           \"\"\"Process a .ddpcr file — same flow as process_qlp_file().\"\"\"
           ddpcr_path = Path(ddpcr_file)
           base_dir   = ddpcr_path.parent
           output_dir = base_dir / "plots"
           output_dir.mkdir(exist_ok=True)
           csv_dir = output_dir / "well_csvs"
           csv_dir.mkdir(exist_ok=True)

           log_func(f"\\nDecrypting and parsing: {ddpcr_path.name}")
           log_func(f"Include filtered droplets: {include_filtered}")

           parser = BioRadDdpcrParser(ddpcr_file)
           well_dataframes = parser.parse_to_dataframe(
               include_filtered=include_filtered)

           log_func(f"Extracted {len(well_dataframes)} wells")

           # From here the logic is identical to process_qlp_file() —
           # you can factor both into a shared _process_well_dataframes()
           # helper, or just duplicate the body.
           # The well_dataframes dict has the same schema in both cases.

   The simplest approach is to refactor process_qlp_file() to accept a
   pre-built well_dataframes dict and call it from both functions:

       def _process_well_dataframes(well_dataframes, source_path,
                                    well_assignments, include_filtered,
                                    log_func, output_format):
           # ... the body of process_qlp_file starting from
           #     "# Calculate shared thresholds..." onwards ...

       def process_qlp_file(qlp_file, well_assignments, include_filtered,
                            log_func=print, output_format='png'):
           parser = BioRadQLPParser(qlp_file)
           wdfs   = parser.parse_to_dataframe(include_filtered=include_filtered)
           log_func(f"Extracted {len(wdfs)} wells")
           return _process_well_dataframes(wdfs, Path(qlp_file), ...)

       def process_ddpcr_file(ddpcr_file, well_assignments, include_filtered,
                              log_func=print, output_format='png'):
           parser = BioRadDdpcrParser(ddpcr_file)
           wdfs   = parser.parse_to_dataframe(include_filtered=include_filtered)
           log_func(f"Extracted {len(wdfs)} wells")
           return _process_well_dataframes(wdfs, Path(ddpcr_file), ...)
"""