# ddPCRvis

Parse and visualize droplet amplitude data from Bio-Rad ddPCR instruments.
Supports `.qlp` binary files (QuantaSoft 1.7), `.ddpcr` encrypted archives
(Droplet Analysis Software), and Bio-Rad CSV exports.

---

## Installation

```bash
pip install pandas numpy matplotlib PyQt5 scikit-learn scipy
```

To open `.ddpcr` files you also need 7-Zip:

```bash
# Ubuntu / Debian
sudo apt install p7zip-full

# macOS
brew install p7zip
```

---

## Quick Start

```bash
python ddPCRvis.py
```

Drop files onto the window, or use the Browse buttons. The tool accepts:

| File type | How it is handled |
|-----------|-------------------|
| `.csv` | Plotted individually and combined |
| `.qlp` | Prompts for well-to-condition assignments |
| `.ddpcr` | Decrypted automatically, then same workflow as QLP |

---

## Programmatic Usage

```python
from ddpcr_parsers import BioRadQLPParser, BioRadDdpcrParser

# QLP
parser = BioRadQLPParser('plate.qlp')
wells = parser.parse_to_dataframe(include_filtered=True)

# ddPCR
parser = BioRadDdpcrParser('plate.ddpcr')
wells = parser.parse_to_dataframe(include_filtered=True)

for well_id, df in wells.items():
    ch = df.attrs['channel_names']       # e.g. ['FAM', 'HEX']
    meta = df.attrs['well_metadata']     # thresholds, concentrations, …
    print(f"{well_id}: {len(df)} droplets, channels={ch}")
```

Both parsers return an identical `dict[str, pd.DataFrame]` — one entry per well
that contained data.

---

## Output

```
plots/
├── well_csvs/
│   ├── A01.csv
│   ├── A02.csv
│   └── …
├── A01_1D_plot.png
├── A01_2D_plot.png
├── …
├── <stem>_96well_plate_<timestamp>.png
└── <stem>_band_statistics_<timestamp>.csv
```

Each per-well CSV contains: `Droplet_Index`, `Timestamp`, `Ch1_Amplitude`,
`Ch2_Amplitude`, `Width`, `Quality`, `Cluster`, `Cluster_Label`, plus all
scalar well metadata broadcast across every row.

---

## DataFrame Schema

### Droplet columns

| Column | Type | Source | Notes |
|--------|------|--------|-------|
| `Droplet_Index` | int | Both | 0-based index within the well |
| `Timestamp` | uint32 | Both | Instrument acquisition timestamp |
| `Ch1_Amplitude` | float | Both | Channel 1 fluorescence amplitude |
| `Ch2_Amplitude` | float | Both | Channel 2 fluorescence amplitude |
| `Ch1_Width` | float | QLP | Droplet width (same value as `Ch2_Width`) |
| `Ch1_Quality` | float | QLP | Droplet quality score (same value as `Ch2_Quality`) |
| `Ch2_Width` | float | QLP | Alias of `Ch1_Width` |
| `Ch2_Quality` | float | QLP | Alias of `Ch1_Quality` |
| `Width` | float | ddPCR | Droplet width |
| `Quality` | float | ddPCR | Droplet quality score |
| `Gating_Flag` | uint32 | Both | Per-droplet gating flag (0 = accepted) |
| `Cluster` | int | Both | Population index 0–4 (see below) |
| `Cluster_Label` | str | Both | Human-readable label (see below) |
| `<Dye>_Amplitude` | float | Both | Named dye columns when dye names are available, e.g. `FAM_Amplitude` |

**Cluster values**

| `Cluster` | Byte (QLP) | Label (QLP) | Label (ddPCR) |
|-----------|------------|-------------|---------------|
| 0 | `0x00` | `Filtered` | `Gated` |
| 1 | `0x11` | `NN` | `FAM-/HEX-` (example) |
| 2 | `0x22` | `Ch1+` | `FAM+/HEX-` (example) |
| 3 | `0x33` | `Ch1+Ch2+` | `FAM+/HEX+` (example) |
| 4 | `0x44` | `Ch2+` | `FAM-/HEX+` (example) |

ddPCR labels are constructed from BioRad's per-target Positive/Negative calls
and actual dye names, so they reflect the specific assay configuration.

### `df.attrs`

| Key | Description |
|-----|-------------|
| `channel_names` | Ordered dye names, e.g. `['FAM', 'HEX']` |
| `channel_map` | Maps amplitude column to dye, e.g. `{'Ch1_Amplitude': 'FAM'}` |
| `well_metadata` | Full per-well metadata dict (superset of broadcast columns) |
| `file_metadata` | Plate-level metadata (software version, instrument, run date) |

### Well metadata keys

| Key | QLP | ddPCR | Description |
|-----|:---:|:-----:|-------------|
| `well_id` | ✓ | ✓ | e.g. `"A01"` |
| `droplet_volume_nL` | ✓ | ✓ | Droplet volume in nanolitres |
| `rejected_droplets` | ✓ | ✓ | Count of rejected droplets |
| `saturated_droplets` | ✓ | ✓ | Count of saturated droplets |
| `flow_rate` | ✓ | — | Droplet flow rate (µL/hr) |
| `sample_name` | ✓ | — | Combined sample ID string (multiple IDs joined with `-`) |
| `supermix` | ✓ | — | Supermix name |
| `target_names` | ✓ | — | Comma-separated target names, e.g. `"Syt1,GAPDH"` |
| `target_ch1` / `target_ch2` | ✓ | — | Split from `target_names` |
| `experiment_type` | ✓ | — | e.g. `"Absolute Quantification"` |
| `experiment_name` | ✓ | — | Experiment name string |
| `threshold_ch1` / `threshold_ch2` | ✓ | — | Auto-threshold per channel |
| `threshold_confidence_ch1/ch2` | ✓ | — | Threshold quality score per channel |
| `manual_threshold_ch1/ch2` | ✓ | — | Manual threshold value, if set |
| `concentration_ch1` / `concentration_ch2` | ✓ | — | Copies/µL (BioRad computed) |
| `conc_ci_lower_ch1/ch2` | ✓ | — | 95% CI lower bound |
| `conc_ci_upper_ch1/ch2` | ✓ | — | 95% CI upper bound |
| `positives_ch1` / `positives_ch2` | ✓ | — | Positive droplet count |
| `negatives_ch1` / `negatives_ch2` | ✓ | — | Negative droplet count |
| `total_accepted_droplets` | ✓ | — | Accepted droplet count |
| `was_thresholded` | ✓ | — | Bool: thresholding was applied |
| `system_version` | ✓ | — | QuantaSoft system version (uint16) |
| `threshold_<DYE>` | — | ✓ | BioRad threshold for each dye |
| `threshold_<DYE>_manual` | — | ✓ | Whether that threshold was manually set |
| `was_thresholded_ch1v2` | — | ✓ | Bool: Ch1 vs Ch2 threshold applied |
| `was_thresholded_ch3v4` | — | ✓ | Bool: Ch3 vs Ch4 threshold applied |
| `was_thresholded_ch5v6` | — | ✓ | Bool: Ch5 vs Ch6 threshold applied |
| `channel_count` | — | ✓ | Number of acquisition channels |

---

## File Format Reference

### QLP

QLP is a **TIFF-derived binary format** written by QuantaSoft 1.7.

#### File header (8 bytes)

```
Offset  Size  Content
0       2     Byte-order marker: b'II' = little-endian, b'MM' = big-endian
2       2     TIFF magic: 0x002A
4       4     Absolute offset of the first IFD
```

Files from QX200/QX600 instruments are almost always little-endian (`II`).

#### IFD (Image File Directory)

Each IFD represents one well; the first IFD is a plate-level header. IFDs are
linked in a singly-linked list.

```
Offset from IFD start  Size      Content
0                      2         Number of tag entries N
2                      N × 12    Tag entries (see below)
2 + N×12               4         Offset to next IFD; 0 = end of chain
```

Each 12-byte tag entry:

```
Bytes 0–1   Tag ID (uint16)
Bytes 2–3   Data type (uint16)
Bytes 4–7   Element count (uint32)
Bytes 8–11  Value inline (if total bytes ≤ 4) or absolute offset to data
```

TIFF data types: `1`=BYTE(1), `2`=ASCII(1), `3`=SHORT(2), `4`=LONG(4),
`5`=RATIONAL(8), `7`=UNDEFINED(1), `9`=SLONG(4), `11`=FLOAT(4), `12`=DOUBLE(8).

#### Tag registry

Standard TIFF tags:

| Tag | Type | Description |
|-----|------|-------------|
| 270 | ASCII | ImageDescription — plate-level XML or JSON |
| 305 | ASCII | Software version string |

Bio-Rad custom tags:

| Tag | Hex | Type | Count | Description |
|-----|-----|------|-------|-------------|
| 65004 | 0xFE04 | ASCII | 1 | Channel names, comma-separated, e.g. `"FAM,HEX"` |
| 65005 | 0xFE05 | FLOAT | 1 | Droplet flow rate (µL/hr) |
| 65006 | 0xFE06 | FLOAT | 2 | Concentration per channel (copies/µL) |
| 65007 | 0xFE07 | FLOAT | 2 | 95% CI lower bound per channel |
| 65008 | 0xFE08 | FLOAT | 2 | 95% CI upper bound per channel |
| 65009 | 0xFE09 | LONG | 2 | Positive droplet count per channel |
| 65010 | 0xFE0A | LONG | 2 | Negative droplet count per channel |
| 65011 | 0xFE0B | LONG | 1 | Total accepted droplet count |
| 65012 | 0xFE0C | FLOAT | 2 | Auto-threshold per channel |
| 65013 | 0xFE0D | FLOAT | 2 | Threshold quality score per channel |
| 65014 | 0xFE0E | FLOAT | 2 | Manual threshold per channel |
| 65015 | 0xFE0F | ASCII | 1 | Sample name (multiple IDs joined with `-`) |
| 65016 | 0xFE10 | ASCII | 1 | Supermix name |
| 65017 | 0xFE11 | ASCII | 1 | Target names, comma-separated |
| 65018 | 0xFE12 | ASCII | 1 | Experiment type string |
| 65019 | 0xFE13 | ASCII | 1 | Well name / ID, e.g. `"A01"` |
| 65020 | 0xFE14 | ASCII | 1 | Experiment name |
| 65021 | 0xFE15 | LONG | 1 | Absolute file offset to first droplet record |
| 65022 | 0xFE16 | LONG | 1 | Rejected droplet count |
| 65023 | 0xFE17 | LONG | 1 | Saturated droplet count |
| 65030 | 0xFE1E | ASCII | 1 | Plate ID (plate-header IFD) |
| 65031 | 0xFE1F | ASCII | 1 | Run date (plate-header IFD) |
| 65032 | 0xFE20 | ASCII | 1 | Instrument serial number (plate-header IFD) |
| 65033 | 0xFE21 | ASCII | 1 | Instrument make (plate-header IFD) |
| 65034 | 0xFE22 | ASCII | 1 | Instrument model (plate-header IFD) |
| 65054 | 0xFE36 | UNDEFINED | NbChannels×N | Per-channel per-droplet gating flags (uint32 each) |
| 65057 | 0xFE39 | BYTE | N | Per-droplet cluster assignment byte array |
| 65058 | 0xFE3A | UNDEFINED | NbChannels | Cluster algorithm mode flags (bool per channel) |
| 65065 | 0xFE41 | FLOAT | 1 | Single-well cluster events confidence |
| 65066 | 0xFE42 | FLOAT | 1 | Multi-well cluster events confidence |
| 65067 | 0xFE43 | BYTE | 1 | Was-thresholded flag (bool) |
| 65074 | 0xFE52 | SHORT | 1 | System version (uint16) |
| 65075 | 0xFE53 | LONG | 1 | Color compensation matrix group index |
| 65078 | 0xFE56 | FLOAT | 1 | Per-well droplet volume (nanolitres) |
| 65079 | 0xFE57 | SHORT | 1 | Well compensation matrix selected mask |
| 65081 | 0xFE59 | LONG | 1 | Production gating mask |

Tags 65022 and 65023 (rejected/saturated counts) and 65031–65034 (run date,
instrument info) are written to plate-header and setup IFDs by QuantaSoft's
native writer. Per-well data IFDs contain structured blobs at some of the same
tag IDs; the parser reads these tags from whichever IFD it finds them in.

#### Droplet record (28 bytes)

Confirmed from the `ManagedQLEvent` struct layout
(`[StructLayout(LayoutKind.Sequential, Pack = 1)]`):

```
Offset  Size  Type    Field
0       4     uint32  Timestamp
4       4     float   Amplitude0  → Ch1_Amplitude
8       4     float   Width0      → Ch1_Width
12      4     float   Quality0    → Ch1_Quality
16      4     float   Amplitude1  → Ch2_Amplitude
20      4     float   Width1      → Ch2_Width
24      4     float   Quality1    → Ch2_Quality
```

`Width0 == Width1` and `Quality0 == Quality1` in every record. They are single
per-droplet scalars that BioRad writes into both channel slots.

Record count equals the length of the cluster byte array (tag 65057).

#### Cluster byte encoding

Confirmed from `ddPcrToQlpConverter.GetEventData`:

| Byte | Decimal | Population |
|------|---------|------------|
| `0x00` | 0 | Gated / Filtered |
| `0x11` | 17 | Ch1−  Ch2−  (NN, double-negative) |
| `0x22` | 34 | Ch1+  Ch2− |
| `0x33` | 51 | Ch1+  Ch2+ (double-positive) |
| `0x44` | 68 | Ch1−  Ch2+ |

The byte array is zero-initialised; droplets not assigned to any cluster default
to `0x00`. Gated clusters (`cluster.IsGated == true`) are also assigned `0x00`.

---

### ddPCR

#### Encryption

A `.ddpcr` file is a **7-Zip archive encrypted with AES-256**. The password is
a static plain-ASCII GUID, confirmed by calling BioRad's
`EncryptDecryptMsgHandler.Decrypt` on the constant
`ComponentIdentityInfo.s_FileIDForRUOFile` with a zero IV:

```
1b53402e-503a-4303-bf86-71af1f3178dd
```

Extract manually with:

```bash
7za x -p"1b53402e-503a-4303-bf86-71af1f3178dd" -o./extracted plate.ddpcr
```

This password is shared across all `.ddpcr` files from the Research edition of
QuantaSoft / Droplet Analysis Software.

#### Archive structure

The archive contains two directories of per-well JSON files — one with raw
amplitude data and one with BioRad's analysis results. Some archives nest these
inside a subdirectory; the parser handles both layouts.

**Amplitude JSON** — one file per well, confirmed keys:

```jsonc
{
  "PeakInfo": {
    "PeakCount":   <int>,
    "Amplitudes":  [[...], [...]],   // [channel_index][droplet_index], floats
    "Widths":      [...],            // float per droplet
    "Timestamps":  [...],            // float per droplet
    "Qualities":   [...],            // float per droplet
    "GatingFlags": [...]             // int per droplet; 0 = accepted
  },
  "DataAcquisitionInfo": {
    "DropletVolume": <float>,
    "ChannelCount":  <int>,
    "ChannelMap": [
      { "Dye": "FAM", "Channel": 1 },
      { "Dye": "HEX", "Channel": 2 }
    ]
  },
  "RejectedInfo": {
    "RejectedDropletCount":  <int>,
    "SaturatedDropletCount": <int>
  }
}
```

**Analysis JSON** — one file per well, confirmed keys:

```jsonc
{
  "WellIndex": <int>,                // 0-based row-major (0=A01, 1=A02, …)
  "Clusters": [
    {
      "Cluster":  <int>,             // population index 0–4 (same as QLP encoding)
      "Targets":  [ { "Name": "Syt1", "Dyes": [{ "Name": "FAM", "Channel": 1 }] } ],
      "Results":  ["Positive", "Negative"],  // one per target
      "Droplets": [0, 5, 12, ...]    // droplet indices for this cluster
    }
  ],
  "ThresholdKeys":   ["FAM", "HEX"],
  "ThresholdValues": [
    [ { "ThresholdValue": 3500.0, "IsManuallySet": false } ],
    [ { "ThresholdValue": 2800.0, "IsManuallySet": false } ]
  ],
  "WasThreshed1v2": <bool>,
  "WasThreshed3v4": <bool>,
  "WasThreshed5v6": <bool>
}
```

`Results` entries are `"Positive"`, `"Negative"`, or `"Gated"`. The parser
builds `Cluster_Label` strings from dye names and results, e.g. `"FAM+/HEX-"`.

---

## License

MIT