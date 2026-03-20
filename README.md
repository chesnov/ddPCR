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
    ch    = df.attrs['channel_names']    # e.g. ['FAM', 'HEX']
    meta  = df.attrs['well_metadata']    # thresholds, concentrations, targets, …
    fmeta = df.attrs['file_metadata']    # instrument info, run dates, …
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
| `<Dye>_Amplitude` | float | Both | Named dye columns when dye names are available |

**Cluster values**

| `Cluster` | Byte (QLP) | Label (QLP) | Label (ddPCR) |
|-----------|------------|-------------|---------------|
| 0 | `0x00` | `Filtered` | `Gated` |
| 1 | `0x11` | `NN` | `FAM-/HEX-` (example) |
| 2 | `0x22` | `Ch1+` | `FAM+/HEX-` (example) |
| 3 | `0x33` | `Ch1+Ch2+` | `FAM+/HEX+` (example) |
| 4 | `0x44` | `Ch2+` | `FAM-/HEX+` (example) |

### `df.attrs`

| Key | Description |
|-----|-------------|
| `channel_names` | Ordered dye names, e.g. `['FAM', 'HEX']` |
| `channel_map` | Maps amplitude column to dye |
| `well_metadata` | Full per-well metadata dict (superset of broadcast columns) |
| `file_metadata` | Plate-level metadata |

### `file_metadata` keys

| Key | QLP | ddPCR | Source | Description |
|-----|:---:|:-----:|--------|-------------|
| `software` | ✓ | — | TAG 305 | QuantaSoft version string |
| `equipment_make` | ✓ | — | TAG 271 | Instrument make |
| `equipment_model` | ✓ | — | TAG 272 | Instrument model |
| `equipment_serial` | ✓ | — | TAG 65033 | Instrument serial number |
| `equipment_droplet_volume_nL` | ✓ | — | TAG 65049 | Equipment-level droplet volume |
| `created_by_user` | — | ✓ | `PersistableHeader.json` → `CreatedByUser` | Operator who created the file |
| `created_date` | — | ✓ | `PersistableHeader.json` → `CreatedDate` | File creation date |
| `created_app_version` | — | ✓ | `PersistableHeader.json` → `CreatedByAppVersion` | Software version |
| `modified_by_user` | — | ✓ | `PersistableHeader.json` → `ModifiedByUserName` | Operator who last saved |
| `modified_date` | — | ✓ | `PersistableHeader.json` → `ModifiedDate` | Last save date |
| `run_start_date` | — | ✓ | `RunInfo.json` → `RunStartDate` | Run start date/time |
| `run_end_date` | — | ✓ | `RunInfo.json` → `RunEndDate` | Run end date/time |
| `plate_file_name` | — | ✓ | `RunInfo.json` → `PlateFileName` | Plate file name |
| `data_file_name` | — | ✓ | `RunInfo.json` → `DataFileName` | Data file name |
| `run_type` | — | ✓ | `RunInfo.json` → `RunType` | Run type string |
| `software_edition` | — | ✓ | `RunInfo.json` → `SoftwareEdition` | e.g. `"Research"` |
| `plate_setup` | — | ✓ | `*.ddplt` | Full raw plate setup JSON |
| `header` | — | ✓ | `PersistableHeader.json` | Full raw header JSON |
| `plate_info` | — | ✓ | `PlateInfo.json` | Full raw instrument/acquisition JSON |
| `run_info` | — | ✓ | `RunInfo.json` | Full raw run info JSON |

### Well metadata keys

| Key | QLP | ddPCR | Source | Description |
|-----|:---:|:-----:|--------|-------------|
| `well_id` | ✓ | ✓ | Both | e.g. `"A01"` |
| `droplet_volume_nL` | ✓ | ✓ | TAG 65078 / `DataAcquisitionInfo` | Per-well droplet volume |
| `flow_rate` | ✓ | — | TAG 65005 | Droplet flow rate (µL/hr) |
| `was_thresholded` | ✓ | — | TAG 65067 | Bool: thresholding was applied |
| `system_version` | ✓ | — | TAG 65074 | QuantaSoft system version (uint16) |
| `concentration_ch1` / `concentration_ch2` | ✓ | — | TAG 65031 `QLBFileQuantitationData.Concentration` | Copies/µL (BioRad computed) |
| `conc_ci_lower_ch1/ch2` | ✓ | — | TAG 65031 `ConfidenceLowerBound` | 95% CI lower bound |
| `conc_ci_upper_ch1/ch2` | ✓ | — | TAG 65031 `ConfidenceUpperBound` | 95% CI upper bound |
| `rejected_droplets` | ✓ | ✓ | TAG 65031 offset−8 / `RejectedInfo` | Rejected droplet count |
| `saturated_droplets` | ✓ | ✓ | TAG 65031 offset−8 / `RejectedInfo` | Saturated droplet count |
| `threshold_ch1` / `threshold_ch2` | ✓ | — | TAG 65032 `QLBFileQuantitationProcessingDetail.Threshold` | Auto-threshold per channel |
| `threshold_confidence_ch1/ch2` | ✓ | — | TAG 65032 `ThresholdConfidence` | Threshold quality score |
| `manual_threshold_ch1/ch2` | ✓ | — | TAG 65032 `ManualThreshold` | Manual threshold value |
| `multi_well_threshold_ch1/ch2` | ✓ | — | TAG 65032 `MultiWellThreshold` | Multi-well threshold |
| `multi_well_manual_threshold_ch1/ch2` | ✓ | — | TAG 65032 `MultiWellManualThreshold` | Multi-well manual threshold |
| `min_width_gate_ch1/ch2` | ✓ | — | TAG 65032 `MinWidthGate` | Minimum width gate |
| `max_width_gate_ch1/ch2` | ✓ | — | TAG 65032 `MaxWidthGate` | Maximum width gate |
| `min_quality_gate_ch1/ch2` | ✓ | — | TAG 65032 `MinQualityGate` | Minimum quality gate |
| `use_auto_threshold_ch1/ch2` | ✓ | — | TAG 65032 `UseAutoThreshold` | Bool: auto threshold in use |
| `use_single_well_ch1/ch2` | ✓ | — | TAG 65032 `UseSingleWell` | Bool: single-well analysis |
| `sample_name` | ✓ | ✓ | Setup blob `SampleName` / `.ddplt` `SampleIds` | Sample ID(s) joined with `-` |
| `experiment_type` | ✓ | ✓ | Setup blob `ExperimentType` / `.ddplt` `OrcaExperimentType` | e.g. `"Absolute Quantification"` |
| `experiment_name` | ✓ | ✓ | Setup blob `ExperimentName` / `.ddplt` `ExperimentName` | Experiment name |
| `experiment_comment` | ✓ | — | Setup blob `ExperimentComment` | Experiment comment |
| `target_ch1` | ✓ | ✓ | Setup blob `Target0` / `.ddplt` targets | Target name for channel 1 |
| `target_ch1_type` | ✓ | ✓ | Setup blob `Type0` / `.ddplt` `TargetType` | `"Reference"` or `"Unknown"` |
| `target_ch2` | ✓ | ✓ | Setup blob `Target1` / `.ddplt` targets | Target name for channel 2 |
| `target_ch2_type` | ✓ | ✓ | Setup blob `Type1` / `.ddplt` `TargetType` | `"Reference"` or `"Unknown"` |
| `target_ch<N>` … | — | ✓ | `.ddplt` targets | Additional targets for 3–6 channel assays |
| `target_ch<N>_is_reference` | — | ✓ | `.ddplt` `TargetType == "Reference"` | Bool |
| `target_ch<N>_dye` | — | ✓ | `.ddplt` `Dye.DyeName` | Dye name, e.g. `"FAM"` |
| `supermix` | ✓ | — | TAG 65077 | Supermix name |
| `droplet_generator_cartridge` | ✓ | — | TAG 65076 | Cartridge name |
| `reaction_volume` | ✓ | — | TAG 65063 | Reaction volume (float) |
| `dilution_factor` | ✓ | — | TAG 65064 | Dilution factor (float) |
| `plex_mode` | — | ✓ | `.ddplt` `PlexMode` | e.g. `"Singleplex"`, `"Duplex"` |
| `assay_name` | — | ✓ | `.ddplt` `Panel.Name` | Panel/assay name |
| `channel_count` | — | ✓ | `DataAcquisitionInfo.ChannelCount` | Number of channels |
| `threshold_<DYE>` | — | ✓ | `.ddmetajson` `ThresholdValues` | User threshold per dye |
| `threshold_<DYE>_manual` | — | ✓ | `.ddmetajson` `ThresholdValues.IsManuallySet` | Whether manually set |
| `auto_threshold_<DYE>` | — | ✓ | `.ddmetajson` `AutoThresholdValues` | BioRad auto-computed threshold per dye |
| `was_thresholded_ch1v2` | — | ✓ | `.ddmetajson` `WasThreshed1v2` | Bool |
| `was_thresholded_ch3v4` | — | ✓ | `.ddmetajson` `WasThreshed3v4` | Bool |
| `was_thresholded_ch5v6` | — | ✓ | `.ddmetajson` `WasThreshed5v6` | Bool |
| `cluster_quality_ch1v2` | — | ✓ | `.ddmetajson` `ClusterQual1v2` | Cluster quality score |
| `cluster_quality_ch3v4` | — | ✓ | `.ddmetajson` `ClusterQual3v4` | Cluster quality score |

---

## File Format Reference

### QLP

QLP is a **TIFF-derived binary format** written by QuantaSoft 1.7's
`ProcessedWriter` class.

#### File header (8 bytes)

```
Offset  Size  Content
0       2     Byte-order marker: b'II' = little-endian, b'MM' = big-endian
2       2     TIFF magic: 0x002A
4       4     Absolute offset to first IFD
```

#### IFD chain structure

The IFD chain contains five distinct IFD types, written in this order:

```
IFD 1        File header IFD    — TAG 305 (software), TAG 65004 (channel names)
IFD 2        Setup IFD          — TAG 65018 blob (QLBFile2ChannelSetupData × N_wells)
IFDs 3…N+2   Additional setup   — one per well: TAGs 65019, 65063, 65064, 65076, 65077
IFD N+3      Equipment IFD      — TAGs 271, 272, 65033, 65049
IFDs N+4…    Well data IFDs     — one per well: TAGs 65021, 65031, 65032, 65054, 65057…
```

Each 12-byte tag entry:
```
Bytes 0–1   Tag ID (uint16)
Bytes 2–3   Data type (uint16)
Bytes 4–7   Element count (uint32)
Bytes 8–11  Value inline (if total bytes ≤ 4) or absolute file offset to data
```

#### Tag registry

**File header / equipment IFD:**

| Tag | Type | Description |
|-----|------|-------------|
| 270 | ASCII | ImageDescription |
| 271 | ASCII | Equipment make (max 15 chars) |
| 272 | ASCII | Equipment model (max 15 chars) |
| 305 | ASCII | Software version |
| 65004 | ASCII | Channel names, comma-separated |
| 65033 | ASCII | Equipment serial number (max 15 chars) |
| 65049 | FLOAT | Equipment droplet volume in nanolitres (inline) |

**Setup IFD:**

| Tag | Type | Count | Description |
|-----|------|-------|-------------|
| 65018 | UNDEFINED | N_wells | `QLBFile2ChannelSetupData` struct array |

**Additional setup IFDs (one per well):**

| Tag | Type | Description |
|-----|------|-------------|
| 65019 | ASCII | Well ID, e.g. `"A01"` (inline, ≤4 bytes) |
| 65063 | FLOAT | Reaction volume (inline) |
| 65064 | FLOAT | Dilution factor (inline) |
| 65076 | ASCII | Droplet generator cartridge name |
| 65077 | ASCII | Supermix name |

**Well data IFDs (one per well):**

| Tag | Type | Description |
|-----|------|-------------|
| 65005 | FLOAT | Droplet flow rate µL/hr (inline) |
| 65019 | ASCII | Well ID (inline) |
| 65021 | UNDEFINED | Absolute offset to first droplet record |
| 65031 | UNDEFINED | Inline ptr to `QLBFileQuantitationData[NbChannels]`; count = NbChannels |
| 65032 | UNDEFINED | Inline ptr to `QLBFileQuantitationProcessingDetail[NbChannels]`; count = NbChannels |
| 65054 | UNDEFINED | Per-channel per-droplet gating flags (uint32 each) |
| 65057 | BYTE | Per-droplet cluster assignment byte array |
| 65058 | UNDEFINED | Cluster algorithm mode flags (bool per channel) |
| 65065 | FLOAT | Single-well cluster events confidence (inline) |
| 65066 | FLOAT | Multi-well cluster events confidence (inline) |
| 65067 | BYTE | Was-thresholded flag (bool, inline) |
| 65074 | SHORT | System version uint16 (inline) |
| 65075 | LONG | Color compensation matrix group index (inline) |
| 65078 | FLOAT | Per-well droplet volume in nanolitres (inline) |
| 65079 | SHORT | Well compensation matrix selected mask (inline) |
| 65081 | LONG | Production gating mask (inline) |

#### Setup blob — `QLBFile2ChannelSetupData`

`[StructLayout(LayoutKind.Sequential, Pack = 1)]`

```
Offset  Size  Field
0       4     NextSetupDataOffset (uint32)
4       1     Format (byte)
5       1     SaveRawData (bool)
6       4     Well (char[4])        — e.g. "A01\0"
10      256   SampleName (char[256])
266     32    ExperimentType (char[32])
298     256   ExperimentName (char[256])
554     256   ExperimentComment (char[256])
810     32    Type0 (char[32])      — "Reference" or "Unknown"
842     256   Target0 (char[256])   — target name for channel 1
1098    32    Type1 (char[32])      \  2-channel only
1130    256   Target1 (char[256])   /
```
Total: 1386 bytes (2-channel), 1098 bytes (1-channel).

#### Quantitation data blob — `QLBFileQuantitationData`

`[StructLayout(LayoutKind.Sequential, Pack = 1)]` — 12 bytes per channel.
Stored at the absolute file offset held in tag 65031's inline value.
`NumberOfRejectedPeaks` and `NumberOfSaturatedPeaks` (uint32 each) are written
immediately before this array (i.e. at the file offset stored in tag 65031 minus 8).

```
Offset  Size  Field
0       4     Concentration (float, copies/µL)
4       4     ConfidenceLowerBound (float)
8       4     ConfidenceUpperBound (float)
```

#### Quantitation processing detail blob — `QLBFileQuantitationProcessingDetail`

`[StructLayout(LayoutKind.Sequential, Pack = 1)]` — 54 bytes per channel.
Stored at the absolute file offset held in tag 65032's inline value.

```
Offset  Size  Field
0       4     WidthGateSigmaMultiplier (float)
4       4     MinWidthGate (float)
8       4     MinWidthGateConfidence (float)
12      4     MaxWidthGate (float)
16      4     MaxWidthGateConfidence (float)
20      4     MinQualityGate (float)
24      4     MinQualityGateConfidence (float)
28      4     Threshold (float)
32      4     ThresholdConfidence (float)
36      4     ManualThreshold (float)
40      4     MultiWellThreshold (float)
44      4     MultiWellThresholdConfidence (float)
48      4     MultiWellManualThreshold (float)
52      1     UseAutoThreshold (bool)
53      1     UseSingleWell (bool)
```

#### Droplet record (28 bytes)

`ManagedQLEvent` — `[StructLayout(LayoutKind.Sequential, Pack = 1)]`:

```
Offset  Type    Field
0       uint32  Timestamp
4       float   Amplitude0  → Ch1_Amplitude
8       float   Width0      → Ch1_Width
12      float   Quality0    → Ch1_Quality
16      float   Amplitude1  → Ch2_Amplitude
20      float   Width1      → Ch2_Width
24      float   Quality1    → Ch2_Quality
```

`Width0 == Width1` and `Quality0 == Quality1` — single per-droplet scalars
duplicated into both channel slots.

For 1-channel instruments the record is 16 bytes (last 12 bytes absent).

#### Cluster byte encoding

| Byte | Decimal | Population |
|------|---------|------------|
| `0x00` | 0 | Gated / Filtered |
| `0x11` | 17 | Ch1−  Ch2−  (NN) |
| `0x22` | 34 | Ch1+  Ch2− |
| `0x33` | 51 | Ch1+  Ch2+ |
| `0x44` | 68 | Ch1−  Ch2+ |

---

### ddPCR

#### Encryption

Password confirmed via `EncryptDecryptMsgHandler.Decrypt` on
`ComponentIdentityInfo.s_FileIDForRUOFile` with key `"DbCdNa2OrCaDx56#4"` and
zero IV:

```
1b53402e-503a-4303-bf86-71af1f3178dd
```

```bash
7za x -p"1b53402e-503a-4303-bf86-71af1f3178dd" -o./extracted plate.ddpcr
```

#### Archive structure

```
<stem>.ddplt                      PlateSetup — sample names, targets, experiment config
PersistableHeader.json            File header — operator, timestamps
PlateInfo.json                    Instrument and acquisition info
RunInfo.json                      Run-level metadata
<stem>.csv                        BioRad's own analysis export
PeakData/<WellID>.ddpeakjson      Per-droplet amplitude data
PeakMetaData/<WellID>.ddmetajson  Cluster assignments and thresholds
```

#### Important note on concentrations

Concentrations, CI bounds, and positive/negative counts are **not stored** in
the ddPCR archive. They are computed by BioRad's `CalculateResults()` at load
time and are not serialised to the ddmetajson.  These values are available only
in the BioRad-generated `.csv` file inside the archive (which is also extracted
to `file_metadata['plate_setup']` indirectly) or by computing them from the
per-droplet cluster assignments.

#### Plate setup file (`*.ddplt`)

```jsonc
{
  "WellSamples": [
    {
      "WellIndex": <int>,               // 0-based row-major
      "SampleIds": ["<id>", ...],       // joined with "-" → sample_name
      "ExperimentName": "<string>",
      "OrcaExperimentType": "<string>",
      "PlexMode": "<string>",
      "Panel": {
        "Name": "<string>",
        "Targets": [
          {
            "TargetName": "<string>",
            "TargetType": "Reference" | "Unknown",
            "IsUIReferenceTarget": <bool>,
            "Dye":  { "DyeName": "<string>", "Channel": <int> },  // 0-based
            "DyeAmount": <float>,
            "Dye2": { "DyeName": "<string>", "Channel": <int> }   // multiplex
          }
        ]
      }
    }
  ]
}
```

#### File header (`PersistableHeader.json`)

```jsonc
{
  "CreatedByUser": "<operator>",    "CreatedDate": "<ISO datetime>",
  "CreatedByAppVersion": "<ver>",   "ModifiedByUserName": "<operator>",
  "ModifiedDate": "<ISO datetime>"
}
```

#### Run info (`RunInfo.json`)

```jsonc
{
  "RunStartDate": "<ISO>",  "RunEndDate": "<ISO>",
  "PlateFileName": "<str>", "DataFileName": "<str>",
  "RunType": "<str>",       "SoftwareEdition": "<str>",
  "QX600InstrumentInfo": { "SerialNumber": "<str>", "FirmwareVersion": "<str>", ... }
}
```

#### Amplitude JSON (`*.ddpeakjson`)

```jsonc
{
  "PeakInfo": {
    "PeakCount": <int>, "Amplitudes": [[...], [...]],
    "Widths": [...], "Timestamps": [...], "Qualities": [...], "GatingFlags": [...]
  },
  "DataAcquisitionInfo": {
    "DropletVolume": <float>, "ChannelCount": <int>,
    "ChannelMap": [ { "Dye": "FAM", "Channel": 1 }, ... ]
  },
  "RejectedInfo": { "RejectedDropletCount": <int>, "SaturatedDropletCount": <int> }
}
```

#### Analysis JSON (`*.ddmetajson`)

Confirmed from `PersistedWellAnalysisResults.ToJObject()`:

```jsonc
{
  "WellIndex": <int>,
  "Clusters": [ { "Cluster": <int>, "Targets": [...], "Results": [...], "Droplets": [...] } ],
  "ThresholdKeys": ["FAM", "HEX"],
  "ThresholdValues": [ [ { "ThresholdValue": 3500.0, "IsManuallySet": false } ], ... ],
  "AutoThresholdKeys": ["FAM", "HEX"],
  "AutoThresholdValues": [ [ { "ThresholdValue": 3480.0, "IsManuallySet": false } ], ... ],
  "WasThreshed1v2": <bool>, "WasThreshed3v4": <bool>, "WasThreshed5v6": <bool>,
  "ClusterQual1v2": <float>, "ClusterQual3v4": <float>,
  "UseClusteringAlg1v2": <bool>, "UseClusteringAlg3v4": <bool>
}
```

---

## License

MIT