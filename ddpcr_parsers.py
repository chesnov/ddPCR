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

    Droplet record: 28 bytes  (confirmed from ManagedQLEvent struct layout)
      uint32  Timestamp
      float   Amplitude0  → Ch1_Amplitude
      float   Width0      → Ch1_Width
      float   Quality0    → Ch1_Quality
      float   Amplitude1  → Ch2_Amplitude
      float   Width1      → Ch2_Width
      float   Quality1    → Ch2_Quality

    Note: Width0==Width1 and Quality0==Quality1 — they are single per-droplet
    scalars that BioRad writes into both channel slots of the record.
    """

    # ── Tag IDs ────────────────────────────────────────────────────────────
    # Standard TIFF
    TAG_IMAGE_DESC          = 270   # Often contains plate XML/JSON
    TAG_EQUIPMENT_MAKE      = 271   # ASCII: instrument make (WriteEquipmentData)
    TAG_EQUIPMENT_MODEL     = 272   # ASCII: instrument model (WriteEquipmentData)
    TAG_SOFTWARE            = 305   # ASCII: software version string

    # Bio-Rad custom — confirmed from ProcessedWriter source
    TAG_CHANNEL_NAMES       = 65004  # ASCII: "Ch1,Ch2" or "FAM,HEX" etc.
    TAG_FLOW_RATE           = 65005  # FLOAT: droplet flow rate (µL/hr)
    # 65006 = ColorCompensationMatrix (FLOAT[4]) — not extracted
    # 65007–65017 = not standalone tags; data lives in blob structs
    TAG_SETUP_BLOB          = 65018  # UNDEFINED blob: QLBFile2ChannelSetupData × N_wells
                                     # Contains: Well, SampleName, ExperimentType,
                                     # ExperimentName, ExperimentComment, Type0/1, Target0/1
    TAG_WELL_NAME           = 65019  # ASCII: well ID e.g. "A01" (inline, ≤4 bytes)
    TAG_DATA_START          = 65021  # UNDEFINED: offset to first droplet record
    # Tag 65031 — QuantitationDataOffset
    # UNDEFINED[NbChannels]: ValueOffset = absolute file ptr to QLBFileQuantitationData[NbChannels]
    # Struct (Pack=1, 12 bytes): float Concentration, ConfidenceLowerBound, ConfidenceUpperBound
    # NumberOfRejectedPeaks (uint32) and NumberOfSaturatedPeaks (uint32) are written
    # immediately before this array (i.e. at offset - 8).
    TAG_QUANT_DATA          = 65031
    # Tag 65032 — QuantitationProcessingDetailsOffset
    # UNDEFINED[NbChannels]: ValueOffset = absolute file ptr to QLBFileQuantitationProcessingDetail[NbChannels]
    # Struct (Pack=1, 54 bytes): 13 floats then 2 bools — see _parse_quantitation_blobs
    TAG_QUANT_PROC_DETAIL   = 65032
    TAG_EQUIPMENT_SERIAL    = 65033  # ASCII: instrument serial (max 15 chars)
    TAG_EQUIPMENT_DROPLET_VOL = 65049 # FLOAT (inline): equipment-level droplet volume (nL)
    TAG_EVENT_GATING_FLAGS  = 65054  # UNDEFINED[NbChannels×N]: per-channel per-droplet gating flags (uint32)
    TAG_CLUSTER_ARRAY       = 65057  # BYTE[N]: per-droplet cluster assignment
    TAG_CLUSTER_MODES       = 65058  # UNDEFINED[NbChannels]: cluster algorithm mode flags
    TAG_REACTION_VOLUME     = 65063  # FLOAT (inline): reaction volume
    TAG_DILUTION_FACTOR     = 65064  # FLOAT (inline): dilution factor
    TAG_EVENT_CLUSTER_CONF  = 65065  # FLOAT (inline): single-well cluster events confidence
    TAG_MULTIWELL_CLUST_CONF= 65066  # FLOAT (inline): multi-well cluster events confidence
    TAG_WAS_THRESHOLDED     = 65067  # BYTE (inline): was-thresholded flag (bool)
    TAG_SYSTEM_VERSION      = 65074  # SHORT (inline): system version uint16
    TAG_COLOR_COMP_GRP_IDX  = 65075  # LONG (inline): color compensation matrix group index
    TAG_DROPLET_CARTRIDGE   = 65076  # ASCII: droplet generator cartridge name
    TAG_SUPERMIX            = 65077  # ASCII: supermix name
    TAG_DROPLET_VOLUME      = 65078  # FLOAT (inline): per-well droplet volume (nL)
    TAG_WELL_COMP_MATRIX    = 65079  # SHORT (inline): well compensation matrix selected mask
    # 65080 = SecondColorCompensationMatrices — not extracted
    TAG_PRODUCTION_GATING   = 65081  # LONG (inline): production gating mask

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

    def _parse_quantitation_blobs(self, tags: dict, m: dict):
        """
        Read QLBFileQuantitationData and QLBFileQuantitationProcessingDetail
        blobs for this well.

        Both tags (65031, 65032) are UNDEFINED type with ValueCount = NbChannels
        (1 or 2).  Since total_sz = NbChannels bytes ≤ 4 the value is stored
        inline in the IFD entry: the uint32 at tags[TAG]['ptr'] IS the absolute
        file offset to the struct array.

        Struct layouts (LayoutKind.Sequential, Pack=1) confirmed from source:

        QLBFileQuantitationData (12 bytes per channel):
          float Concentration
          float ConfidenceLowerBound
          float ConfidenceUpperBound

        QLBFileQuantitationProcessingDetail (54 bytes per channel):
          float WidthGateSigmaMultiplier
          float MinWidthGate
          float MinWidthGateConfidence
          float MaxWidthGate
          float MaxWidthGateConfidence
          float MinQualityGate
          float MinQualityGateConfidence
          float Threshold
          float ThresholdConfidence
          float ManualThreshold
          float MultiWellThreshold
          float MultiWellThresholdConfidence
          float MultiWellManualThreshold
          bool  UseAutoThreshold           (1 byte)
          bool  UseSingleWell              (1 byte)

        NumberOfRejectedPeaks (uint32) and NumberOfSaturatedPeaks (uint32)
        are written immediately before the QuantitationData array — confirmed
        from WriteWellData: UserInterfaceOffset = QuantitationDataOffset +
        NbChannels * sizeof(QLBFileQuantitationData) + 8.
        """
        ch_labels = ['ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6']

        # ── QuantitationData (tag 65031) ──────────────────────────────────
        if self.TAG_QUANT_DATA in tags:
            t      = tags[self.TAG_QUANT_DATA]
            n_ch   = t['count']                  # ValueCount = NbChannels
            f_off  = self._u(t['ptr'], 'I')      # inline uint32 = absolute file offset
            QD_FMT = f'{self.endian}3f'
            QD_SZ  = 12
            if f_off and f_off + n_ch * QD_SZ <= len(self.data):
                for i in range(n_ch):
                    try:
                        conc, ci_lo, ci_hi = struct.unpack_from(
                            QD_FMT, self.data, f_off + i * QD_SZ)
                        ch = ch_labels[i] if i < len(ch_labels) else f'ch{i+1}'
                        m[f'concentration_{ch}'] = conc
                        m[f'conc_ci_lower_{ch}'] = ci_lo
                        m[f'conc_ci_upper_{ch}'] = ci_hi
                    except struct.error:
                        break

                # Rejected/saturated counts sit 8 bytes before QuantitationData
                rej_off = f_off - 8
                if rej_off >= 0 and rej_off + 8 <= len(self.data):
                    try:
                        rej, sat = struct.unpack_from(
                            f'{self.endian}II', self.data, rej_off)
                        m['rejected_droplets']  = rej
                        m['saturated_droplets'] = sat
                    except struct.error:
                        pass

        # ── QuantitationProcessingDetail (tag 65032) ──────────────────────
        if self.TAG_QUANT_PROC_DETAIL in tags:
            t      = tags[self.TAG_QUANT_PROC_DETAIL]
            n_ch   = t['count']
            f_off  = self._u(t['ptr'], 'I')
            # Pack=1: 13 floats (52 bytes) + 2 bools (2 bytes) = 54 bytes
            QP_FMT = f'{self.endian}13f2B'
            QP_SZ  = 54
            if f_off and f_off + n_ch * QP_SZ <= len(self.data):
                for i in range(n_ch):
                    try:
                        fields = struct.unpack_from(
                            QP_FMT, self.data, f_off + i * QP_SZ)
                        ch = ch_labels[i] if i < len(ch_labels) else f'ch{i+1}'
                        (wg_sigma, min_wg, min_wg_conf, max_wg, max_wg_conf,
                         min_qual, min_qual_conf, thresh, thresh_conf,
                         manual_thresh, mw_thresh, mw_thresh_conf,
                         mw_manual_thresh, use_auto, use_single) = fields
                        m[f'threshold_{ch}']                   = thresh
                        m[f'threshold_confidence_{ch}']        = thresh_conf
                        m[f'manual_threshold_{ch}']            = manual_thresh
                        m[f'multi_well_threshold_{ch}']        = mw_thresh
                        m[f'multi_well_manual_threshold_{ch}'] = mw_manual_thresh
                        m[f'min_width_gate_{ch}']              = min_wg
                        m[f'max_width_gate_{ch}']              = max_wg
                        m[f'min_quality_gate_{ch}']            = min_qual
                        m[f'use_auto_threshold_{ch}']          = bool(use_auto)
                        m[f'use_single_well_{ch}']             = bool(use_single)
                    except struct.error:
                        break

    def _extract_well_metadata(self, tags: dict, well_id: str) -> dict:
        """
        Extract metadata from a well data IFD (the IFD that has TAG_DATA_START).

        Only reads tags confirmed to exist in well data IFDs from
        ProcessedWriter.WriteWellData source. Setup metadata (sample name,
        targets, experiment type, supermix, etc.) is NOT in well data IFDs —
        it comes from separate setup IFDs and is merged in parse_to_dataframe.
        """
        m = {'well_id': well_id}

        # Confirmed from WriteWellData: FlowRateEntry → tag 65005
        fr = self._t_float(tags, self.TAG_FLOW_RATE)
        if fr is not None:
            m['flow_rate'] = fr

        # Confirmed from WriteWellData: WasThresholdedEntry → tag 65067
        wt = self._t_uint32(tags, self.TAG_WAS_THRESHOLDED)
        if wt is not None:
            m['was_thresholded'] = bool(wt)

        # Confirmed from WriteWellData: DropletVolumeEntry → tag 65078
        dv = self._t_float(tags, self.TAG_DROPLET_VOLUME)
        if dv is not None:
            m['droplet_volume_nL'] = dv

        # Confirmed from WriteWellData: SystemVersionEntry → tag 65074
        if self.TAG_SYSTEM_VERSION in tags:
            sv = self._u(tags[self.TAG_SYSTEM_VERSION]['ptr'], "H")
            if sv is not None:
                m['system_version'] = sv

        # Concentrations, CI bounds, thresholds, width/quality gates,
        # and rejected/saturated counts from quantitation blob structs
        self._parse_quantitation_blobs(tags, m)

        # Channel names (also present in well IFDs in some versions)
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

        gating_blob = None
        if self.TAG_EVENT_GATING_FLAGS in tags:
            gi = tags[self.TAG_EVENT_GATING_FLAGS]
            if gi['size'] > 0:
                gating_blob = self.data[gi['ptr']:gi['ptr'] + gi['size']]

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
            # raw = (Timestamp, Amplitude0, Width0, Quality0,
            #                   Amplitude1, Width1, Quality1)
            # Width0==Width1 and Quality0==Quality1 (single scalar, duplicated)
            cluster_byte = clust_blob[i] if i < len(clust_blob) else 0
            cluster      = self.CLUSTER_MAP.get(cluster_byte, 0)

            # EventGatingFlags: uint32[NbChannels][N]; read channel-0 flag for droplet i
            gating_flag = 0
            if gating_blob:
                gf_offset = i * 4
                if gf_offset + 4 <= len(gating_blob):
                    gating_flag = struct.unpack_from(f'{self.endian}I',
                                                     gating_blob, gf_offset)[0]

            if not include_filtered and cluster == 0:
                cursor += self.RECORD_SIZE
                continue

            rec = {
                'Droplet_Index': i,
                'Timestamp':     raw[0],
                'Ch1_Amplitude': raw[1],
                'Ch1_Width':     raw[2],
                'Ch1_Quality':   raw[3],
                'Ch2_Amplitude': raw[4],
                'Ch2_Width':     raw[5],
                'Ch2_Quality':   raw[6],
                'Cluster':       cluster,
                'Cluster_Label': self.CLUSTER_LABELS.get(cluster, 'Unknown'),
                'Gating_Flag':   gating_flag,
            }
            rec.update(meta_cols)
            records.append(rec)
            cursor += self.RECORD_SIZE

        return records

    # ── Setup IFD parsers ─────────────────────────────────────────────────

    def _extract_equipment_ifd(self, tags: dict):
        """
        Extract instrument metadata from the equipment IFD written by
        ProcessedWriter.WriteEquipmentData.
        Tags confirmed: 271 (Make), 272 (Model), 65033 (Serial), 65049 (DropletVol).
        """
        make = self._t_ascii(tags, self.TAG_EQUIPMENT_MAKE)
        if make:
            self.metadata['equipment_make'] = make

        model = self._t_ascii(tags, self.TAG_EQUIPMENT_MODEL)
        if model:
            self.metadata['equipment_model'] = model

        serial = self._t_ascii(tags, self.TAG_EQUIPMENT_SERIAL)
        if serial:
            self.metadata['equipment_serial'] = serial

        dv = self._t_float(tags, self.TAG_EQUIPMENT_DROPLET_VOL)
        if dv is not None:
            self.metadata['equipment_droplet_volume_nL'] = dv

    def _parse_setup_blob(self, tags: dict, setup_by_well: dict):
        """
        Parse the QLBFile2ChannelSetupData struct blob at tag 65018.

        Struct layout (LayoutKind.Sequential, Pack=1) confirmed from source:
          uint32    NextSetupDataOffset       4
          byte      Format                    1
          bool      SaveRawData               1
          char[4]   Well                      4   (ByValTStr, SizeConst=4)
          char[256] SampleName              256
          char[32]  ExperimentType           32
          char[256] ExperimentName          256
          char[256] ExperimentComment       256
          char[32]  Type0                    32
          char[256] Target0                 256
          char[32]  Type1                    32   (2-channel only)
          char[256] Target1                 256   (2-channel only)
        Total 2ch: 1386 bytes  |  1ch: 1098 bytes (no Type1/Target1)
        """
        if self.TAG_SETUP_BLOB not in tags:
            return
        t        = tags[self.TAG_SETUP_BLOB]
        n_wells  = t['count']
        offset   = t['ptr']

        FMT_2CH  = f'{self.endian}IBB4s256s32s256s256s32s256s32s256s'
        FMT_1CH  = f'{self.endian}IBB4s256s32s256s256s32s256s'
        SZ_2CH   = struct.calcsize(FMT_2CH)   # 1386
        SZ_1CH   = struct.calcsize(FMT_1CH)   # 1098

        if n_wells <= 0:
            return

        # Infer channel count from blob size; fall back to 2-channel
        blob_size = t['size']
        per_well  = blob_size // n_wells if n_wells else SZ_2CH
        if per_well == SZ_1CH:
            fmt, sz, is_2ch = FMT_1CH, SZ_1CH, False
        else:
            fmt, sz, is_2ch = FMT_2CH, SZ_2CH, True

        def decode(b: bytes) -> str:
            return b.split(b'\x00')[0].decode('ascii', errors='ignore').strip()

        cursor = offset
        for _ in range(n_wells):
            if cursor + sz > len(self.data):
                break
            try:
                fields = struct.unpack_from(fmt, self.data, cursor)
            except struct.error:
                break

            # Unpack positional fields
            # (NextOffset, Format, SaveRaw, Well, SampleName, ExpType, ExpName,
            #  ExpComment, Type0, Target0[, Type1, Target1])
            well_raw    = decode(fields[3])
            sample_name = decode(fields[4])
            exp_type    = decode(fields[5])
            exp_name    = decode(fields[6])
            exp_comment = decode(fields[7])
            type0       = decode(fields[8])
            target0     = decode(fields[9])
            type1       = decode(fields[10]) if is_2ch else ''
            target1     = decode(fields[11]) if is_2ch else ''

            if _valid_well_id(well_raw):
                setup_by_well[well_raw] = {
                    'sample_name':        sample_name,
                    'experiment_type':    exp_type,
                    'experiment_name':    exp_name,
                    'experiment_comment': exp_comment,
                    'target_ch1':         target0,
                    'target_ch1_type':    type0,
                    'target_ch2':         target1,
                    'target_ch2_type':    type1,
                }
            cursor += sz

    def _parse_additional_setup_ifd(self, tags: dict, setup_by_well: dict):
        """
        Merge additional per-well fields from an IFD written by
        ProcessedWriter.WriteAdditionalSetupInfo.
        Tags confirmed: 65019 (well name), 65077 (supermix), 65076 (cartridge),
                        65063 (reaction volume), 65064 (dilution factor).
        """
        well_id = self._t_ascii(tags, self.TAG_WELL_NAME)
        if not _valid_well_id(well_id):
            return
        if well_id not in setup_by_well:
            setup_by_well[well_id] = {}

        supermix = self._t_ascii(tags, self.TAG_SUPERMIX)
        if supermix:
            setup_by_well[well_id]['supermix'] = supermix

        cartridge = self._t_ascii(tags, self.TAG_DROPLET_CARTRIDGE)
        if cartridge:
            setup_by_well[well_id]['droplet_generator_cartridge'] = cartridge

        rv = self._t_float(tags, self.TAG_REACTION_VOLUME)
        if rv is not None:
            setup_by_well[well_id]['reaction_volume'] = rv

        df_val = self._t_float(tags, self.TAG_DILUTION_FACTOR)
        if df_val is not None:
            setup_by_well[well_id]['dilution_factor'] = df_val

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
            well_metadata   : full per-well metadata dict
            file_metadata   : plate-level metadata dict

        The IFD chain written by ProcessedWriter contains four IFD types:
          1. File header IFD (constructor) — software, channel names
          2. Setup IFD (WriteSetupData) — tag 65018 blob with sample/target/experiment data
          3. Additional setup IFDs × N (WriteAdditionalSetupInfo) — supermix, cartridge, etc.
          4. Equipment IFD (WriteEquipmentData) — make, model, serial, droplet volume
          5. Well data IFDs × N (WriteWellData) — droplet records, cluster calls, thresholds
        """
        self._extract_file_metadata()

        # ── Pass 1: collect setup and equipment data ───────────────────────
        setup_by_well: Dict[str, dict] = {}
        ifd_offset = self._u(4, "I")
        while ifd_offset and ifd_offset < len(self.data) - 6:
            tags, next_ifd = self._parse_ifd(ifd_offset)

            # Equipment IFD: has standard TIFF Make (271) or Equipment Serial (65033),
            # and does NOT have a droplet data blob (TAG_DATA_START / TAG_SETUP_BLOB)
            if ((self.TAG_EQUIPMENT_MAKE in tags or self.TAG_EQUIPMENT_SERIAL in tags) and
                    self.TAG_DATA_START not in tags and self.TAG_SETUP_BLOB not in tags):
                self._extract_equipment_ifd(tags)

            # Setup IFD: has the per-well struct blob (tag 65018)
            elif self.TAG_SETUP_BLOB in tags:
                self._parse_setup_blob(tags, setup_by_well)

            # Additional setup IFD: has supermix (65077) or cartridge (65076),
            # no data blob and no setup blob
            elif ((self.TAG_SUPERMIX in tags or self.TAG_DROPLET_CARTRIDGE in tags) and
                    self.TAG_DATA_START not in tags and self.TAG_SETUP_BLOB not in tags):
                self._parse_additional_setup_ifd(tags, setup_by_well)

            ifd_offset = next_ifd

        # ── Pass 2: build per-well DataFrames ─────────────────────────────
        result: Dict[str, pd.DataFrame] = {}
        ifd_offset = self._u(4, "I")
        while ifd_offset and ifd_offset < len(self.data) - 6:
            tags, next_ifd = self._parse_ifd(ifd_offset)

            is_well = (self.TAG_DATA_START    in tags and
                       self.TAG_CLUSTER_ARRAY in tags and
                       tags[self.TAG_CLUSTER_ARRAY]['size'] > 0)

            if is_well:
                well_id = self._t_ascii(tags, self.TAG_WELL_NAME)
                if not _valid_well_id(well_id):
                    well_id = self._fallback_well_id(ifd_offset)

                if _valid_well_id(well_id):
                    well_meta = self._extract_well_metadata(tags, well_id)
                    # Merge setup data collected in pass 1
                    if well_id in setup_by_well:
                        well_meta.update(setup_by_well[well_id])
                    self.well_metadata[well_id] = well_meta

                    droplets = self._parse_droplets(tags, well_meta,
                                                    include_filtered)
                    if droplets:
                        df = pd.DataFrame(droplets)
                        _attach_attrs(df, self.channel_names or ['Ch1', 'Ch2'],
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
        self._well_samples: Dict[int, dict] = {}   # well_index → parsed WellSample

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

        Archive structure (from JsonPlateFileReaderWriter.ReadPlate):
          <stem>.ddplt                   — plate setup (samples, targets, experiment)
          PersistableHeader.json         — file header (operator name, timestamps)
          PlateInfo.json                 — instrument / run acquisition info
          RunInfo.json                   — run-level metadata
          QX600.SignalProcessing.log     — instrument signal processing log
          QXMgrStandard.Run.log          — instrument run log
          <stem>.csv                     — BioRad's own analysis export
          PeakData/<WellID>.ddpeakjson   — per-droplet amplitudes + gating
          PeakMetaData/<WellID>.ddmetajson — cluster assignments + thresholds
          (plus chart option JSONs, audit log, calibration files — not extracted)

        Plate-level files are stored as raw JSON in df.attrs['file_metadata']
        under the keys 'plate_setup', 'header', 'plate_info', and 'run_info'.
        """
        self._extract()
        try:
            return self._parse_extracted(include_filtered)
        finally:
            self.cleanup()

    def _load_plate_files(self, root: Path):
        """
        Scan the archive root and load plate-level files, mirroring the logic
        in JsonPlateFileReaderWriter.ReadPlate() which iterates
        Directory.GetFiles(pathName) and dispatches by extension/name.

        Loaded files are stored as raw parsed JSON in self.metadata under:
          'plate_setup'  — from *.ddplt / *.ddplate  (PlateSetup.ToJObject)
          'header'       — from PersistableHeader.json
          'plate_info'   — from PlateInfo.json        (DNA2PlateInfoFile)
          'run_info'     — from RunInfo.json           (RunData.Load / RunData.Save)

        No field-level interpretation is done here because the schemas of these
        files are defined by BioRad classes (PlateSetup, PersistableHeader,
        DNA2PlateInfoFile, RunInfo) whose source is not yet available.
        """
        for p in root.iterdir():
            if not p.is_file():
                continue
            name_lower = p.name.lower()

            # PlateSetup — *.ddplt or *.ddplate
            if name_lower.endswith('.ddplt') or name_lower.endswith('.ddplate'):
                try:
                    self.metadata['plate_setup'] = self._load_json(p)
                    if self.debug:
                        print(f"ddpcr: loaded plate_setup from {p.name}")
                except Exception as e:
                    if self.debug:
                        print(f"ddpcr: failed to load plate_setup: {e}")

            # PersistableHeader — contains at minimum ModifiedByUserName
            elif name_lower == 'persistableheader.json':
                try:
                    self.metadata['header'] = self._load_json(p)
                    if self.debug:
                        print(f"ddpcr: loaded header from {p.name}")
                except Exception as e:
                    if self.debug:
                        print(f"ddpcr: failed to load header: {e}")

            # DNA2PlateInfoFile — instrument / acquisition run info
            elif name_lower == 'plateinfo.json':
                try:
                    self.metadata['plate_info'] = self._load_json(p)
                    if self.debug:
                        print(f"ddpcr: loaded plate_info from {p.name}")
                except Exception as e:
                    if self.debug:
                        print(f"ddpcr: failed to load plate_info: {e}")

            # RunInfo — run-level metadata written by RunData.Save / read by RunData.Load
            elif name_lower == 'runinfo.json':
                try:
                    self.metadata['run_info'] = self._load_json(p)
                    if self.debug:
                        print(f"ddpcr: loaded run_info from {p.name}")
                except Exception as e:
                    if self.debug:
                        print(f"ddpcr: failed to load run_info: {e}")

        # Index per-well sample data from plate setup
        if 'plate_setup' in self.metadata:
            self._index_well_samples(self.metadata['plate_setup'])

        # Extract scalar plate-level metadata from header and run_info
        self._extract_scalar_metadata()

    def _index_well_samples(self, plate_setup: dict):
        """
        Build self._well_samples: {well_index → dict} from the WellSamples
        array in PlateSetup.ToJObject() / GetJObjectFromWellSample().

        Field names confirmed from source:
          PlateSetup.GetJObjectFromWellSample  — top-level WellSample fields
          PlateSetup.GetJObjectFromTarget      — TargetName, TargetType,
                                                 IsUIReferenceTarget, Dye/Dye2…
          PlateSetup.GetJObjectFromDye         — DyeName, Channel (0-based)
          PlateSetup.WriteDyeObject            — writes Dye1–6 + DyeAmount1–6
        """
        for ws in plate_setup.get('WellSamples', []):
            if not isinstance(ws, dict):
                continue
            idx = ws.get('WellIndex')
            if idx is None:
                continue

            # Sample name: SampleIds joined with "-"
            # Mirrors GetCombinedSampleID: joins wellSample.SampleIds with "-"
            sample_ids = ws.get('SampleIds', [])
            sample_name = '-'.join(s for s in sample_ids if s) if sample_ids else ''

            # Targets from Panel (GetJObjectFromAssay / GetJObjectFromTarget)
            panel   = ws.get('Panel', {}) if isinstance(ws.get('Panel'), dict) else {}
            # targets_by_channel: {0-based channel index → target dict}
            targets_by_channel: Dict[int, dict] = {}
            for t in panel.get('Targets', []):
                if not isinstance(t, dict):
                    continue
                # Primary dye determines the channel for this target.
                # GetJObjectFromDye writes: {"DyeName": ..., "Channel": ChannelIndex}
                # Channel is 0-based (dye.ChannelIndex).
                dye_obj = t.get('Dye')
                if isinstance(dye_obj, dict):
                    ch = dye_obj.get('Channel')  # 0-based
                    dye_name = dye_obj.get('DyeName', '')
                else:
                    ch = None
                    dye_name = ''

                target = {
                    'name':              t.get('TargetName', ''),
                    'type':              t.get('TargetType', ''),
                    'is_reference':      t.get('TargetType') == 'Reference',
                    'is_ui_reference':   bool(t.get('IsUIReferenceTarget', False)),
                    'dye_name':          dye_name,
                    'channel':           ch,   # 0-based; None if no dye set
                }

                if ch is not None:
                    targets_by_channel[int(ch)] = target
                elif targets_by_channel is not None:
                    # No dye assigned — store by insertion order as fallback
                    fallback_ch = len(targets_by_channel)
                    targets_by_channel[fallback_ch] = target

            self._well_samples[int(idx)] = {
                'sample_name':       sample_name,
                'experiment_name':   ws.get('ExperimentName', ''),
                'experiment_type':   ws.get('OrcaExperimentType', ''),
                'plex_mode':         ws.get('PlexMode', ''),
                'assay_name':        panel.get('Name', ''),
                'targets_by_channel': targets_by_channel,
            }

    def _extract_scalar_metadata(self):
        """
        Promote a flat set of scalar fields from PersistableHeader.json and
        RunInfo.json into self.metadata for inclusion in df.attrs['file_metadata'].

        Field names confirmed from PersistableHeader.ToJObject() and
        RunInfo.FromJObject() source.
        """
        # PersistableHeader fields
        header = self.metadata.get('header') or {}
        for dest, src in [
            ('created_by_user',     'CreatedByUser'),
            ('created_date',        'CreatedDate'),
            ('created_app_version', 'CreatedByAppVersion'),
            ('modified_by_user',    'ModifiedByUserName'),
            ('modified_date',       'ModifiedDate'),
        ]:
            v = header.get(src)
            if v:
                self.metadata[dest] = v

        # RunInfo fields
        run_info = self.metadata.get('run_info') or {}
        for dest, src in [
            ('run_start_date',   'RunStartDate'),
            ('run_end_date',     'RunEndDate'),
            ('plate_file_name',  'PlateFileName'),
            ('data_file_name',   'DataFileName'),
            ('run_type',         'RunType'),
            ('software_edition', 'SoftwareEdition'),
        ]:
            v = run_info.get(src)
            if v:
                self.metadata[dest] = v

    def _parse_extracted(self, include_filtered: bool) -> Dict[str, pd.DataFrame]:
        # ── Locate PeakData directory (handle optional nesting) ───────────
        root     = self._temp_dir
        peak_dir = root / 'PeakData'
        meta_dir = root / 'PeakMetaData'

        if not peak_dir.exists():
            subdirs = [p for p in root.iterdir() if p.is_dir()]
            for sd in subdirs:
                if (sd / 'PeakData').exists():
                    root     = sd
                    peak_dir = sd / 'PeakData'
                    meta_dir = sd / 'PeakMetaData'
                    break

        if not peak_dir.exists():
            contents = list(self._temp_dir.rglob('*'))[:20]
            raise RuntimeError(
                f"PeakData/ directory not found inside archive. "
                f"Archive contents: {contents}"
            )

        # ── Load plate-level files from archive root ──────────────────────
        self._load_plate_files(root)

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
                # Resolve 0-based row-major well index from well_id (e.g. "B03" → 14)
                well_index  = (ord(well_id[0]) - ord('A')) * 12 + (int(well_id[1:]) - 1)
                well_sample = self._well_samples.get(well_index)

                peak = self._load_json(peak_path)
                meta = (self._load_json(meta_files[well_id])
                        if well_id in meta_files else None)
                df = self._build_well_df(well_id, peak, meta, include_filtered,
                                         well_sample)
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
                       include_filtered: bool,
                       well_sample: Optional[dict] = None) -> Optional[pd.DataFrame]:
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

        # ── Plate-setup metadata (from .ddplt WellSamples) ────────────────
        # Field names confirmed from PlateSetup.GetJObjectFromWellSample,
        # GetJObjectFromTarget, and GetJObjectFromDye source.
        if well_sample:
            well_meta['sample_name']     = well_sample.get('sample_name', '')
            well_meta['experiment_name'] = well_sample.get('experiment_name', '')
            well_meta['experiment_type'] = well_sample.get('experiment_type', '')
            well_meta['plex_mode']       = well_sample.get('plex_mode', '')
            well_meta['assay_name']      = well_sample.get('assay_name', '')
            # targets_by_channel: {0-based channel index → target dict}
            # Dye.Channel is 0-based (confirmed from GetJObjectFromDye).
            tbc = well_sample.get('targets_by_channel', {})
            for ch_idx, t in tbc.items():
                ch_label = f'ch{ch_idx + 1}'   # convert to 1-based label
                well_meta[f'target_{ch_label}']              = t.get('name', '')
                well_meta[f'target_{ch_label}_type']         = t.get('type', '')
                well_meta[f'target_{ch_label}_is_reference'] = t.get('is_reference', False)
                well_meta[f'target_{ch_label}_dye']          = t.get('dye_name', '')

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

        # ── Thresholds (user-facing) ───────────────────────────────────────
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

        # ── Auto-thresholds (BioRad's computed values, separate from user thresholds)
        # Confirmed from PersistedWellAnalysisResults.ToJObject: AutoThresholdKeys/AutoThresholdValues
        auto_keys   = meta.get('AutoThresholdKeys',   [])
        auto_values = meta.get('AutoThresholdValues', [])
        for dye, tv_list in zip(auto_keys, auto_values):
            if isinstance(tv_list, list) and tv_list:
                tv = tv_list[0]
                if isinstance(tv, dict) and 'ThresholdValue' in tv:
                    safe_dye = dye.replace(' ', '_').replace('.', '_')
                    well_meta[f'auto_threshold_{safe_dye}'] = float(tv['ThresholdValue'])

        # ── was_thresholded flags ──────────────────────────────────────────
        well_meta['was_thresholded_ch1v2'] = bool(meta.get('WasThreshed1v2', False))
        well_meta['was_thresholded_ch3v4'] = bool(meta.get('WasThreshed3v4', False))
        well_meta['was_thresholded_ch5v6'] = bool(meta.get('WasThreshed5v6', False))

        # ── Cluster quality scores ─────────────────────────────────────────
        # Confirmed from PersistedWellAnalysisResults.ToJObject: ClusterQual1v2, ClusterQual3v4
        cq12 = meta.get('ClusterQual1v2')
        if cq12 is not None:
            well_meta['cluster_quality_ch1v2'] = cq12
        cq34 = meta.get('ClusterQual3v4')
        if cq34 is not None:
            well_meta['cluster_quality_ch3v4'] = cq34

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
