# ddPCRvis - Bio-Rad QLP File Parser & Visualization Tool

Parse and visualize droplet amplitude data from Bio-Rad ddPCR `.qlp` files.

## Quick Start

```bash
pip install pandas numpy matplotlib PyQt5
python ddPCRvis.py
```

---

## Features

- Parse `.qlp` binary files from Bio-Rad QX200/QX600
- Extract well IDs automatically (e.g., C03, D05)
- Extract channel names (e.g., FAM, HEX, ZEB2)
- CSV exports with actual channel names
- 1D/2D amplitude plots  
- 96-well plate overview
- Band statistics

---

## Usage

### GUI
1. Run `python ddPCRvis.py`
2. Select QLP files
3. Assign conditions (optional)
4. Process

### Programmatic
```python
from ddPCRvis import BioRadQLPParser

parser = BioRadQLPParser('file.qlp')
wells = parser.parse_to_dataframe()

for well_id, df in wells.items():
    print(f"{well_id}: {df.attrs['channel_names']}")
```

---

## Output

```
plots/
├── well_csvs/
│   ├── C03.csv           # FAM_Amplitude, HEX_Amplitude, Cluster, Quality_Flag
│   └── ...
├── C03_1D_plot.pdf
├── C03_2D_plot.pdf
├── 96well_overview.pdf
└── band_statistics.csv
```

---

## QLP File Format (What We Know)

> ⚠️ **WARNING**: This is from reverse-engineering. It's **INCOMPLETE** and **MAY BE WRONG**. Bio-Rad hasn't documented this format.

> ✅ **UPDATE**: Well ID extraction is now **reliable** - we discovered they're at a fixed offset (+10) in each IFD header, not scattered through the file. This eliminates previous ordering assumptions.

### Structure
```
QLP File = TIFF-like structure with custom tags
├── TIFF Header (8 bytes: "II" + magic + IFD offset)
└── IFD Chain (linked list of Image File Directories)
```

### What Works ✅

#### 1. Well IDs: Direct Read from IFD Structure
```
IFD Header Pattern (each well block):
Offset: +0  +2      +6          +10        +14
        1C 00 FB FD 02 00 01 00 00 00 [WELL_ID]\x00 ED FD
        ││││  │││││                    └─────┬─────┘
        ││││  │││││                          │
        ││││  │││││                     "C03", "D04", etc.
        ││││  │││││
        ││││  │││││ (This is Tag 65531/0xFDFB entry)
        ││││  ││││└─ Tag 65533 (0xFDFD) follows
        ││││  │││└── Count: 1
        ││││  ││└─── Type: 2 (ASCII)
        ││││  │└──── Tag ID: FB FD (little-endian)
        ││││  └───── (part of tag structure)
        │││└──────── Number of tag entries (0x1C = 28)
        └─────────── IFD marker
```

**How it works:**
- Well ID is **embedded in IFD header** at fixed offset +10
- Read 4 bytes: `[Row][Col1][Col2]\x00` (e.g., "C03\x00")
- No pattern matching needed - just read directly!
- Row = A-H, Col = 01-12

**Why this is reliable:**
- Fixed position in structured data (not floating in binary)
- No false positives from data sections
- No dependency on file ordering
- Works even with empty wells (cluster_size = 0)

**Implementation:**
```python
well_id_offset = ifd_offset + 10
well_id_bytes = data[well_id_offset:well_id_offset+4]
well_id = well_id_bytes.split(b'\x00')[0].decode('ascii')
# Returns: "C03", "D04", etc.
```

#### 2. Channel Names: `[zeros] "Unknown\x00" [+32] [Name]`
```
...00 00 55 6E 6B 6E 6F 77 6E 00...46 41 4D 00
      └──────┬─────┘              └───┬───┘
         "Unknown"                  "FAM"
```
- Search for "Unknown\x00" with zero padding
- Read +32 bytes ahead for channel name
- Take first 2 unique names

**Unknown**: What is "Unknown"? Why +32 offset?

#### 3. Droplet Records (28 bytes)
```
[ID:4][Ch1:4][???:4][???:4][Ch2:4][???:4][???:4]
      └float        └float
```
- Byte 4-7: Ch1 amplitude ✅
- Byte 16-19: Ch2 amplitude ✅
- Rest: **UNKNOWN**

#### 4. Clusters (1 byte/droplet)
```
0x00 → 0 (filtered)
0x11 → 1
0x22 → 2
0x33 → 3
0x44 → 4
```
**Unknown**: Why these values? More clusters?

### Known Tags

| Tag | Hex | Content | Status |
|-----|-----|---------|--------|
| 305 | 0x0131 | Software version | ✅ |
| 65004 | 0xFDF4 | "Ch1,Ch2" (generic) | ⚠️ Not useful |
| 65019 | 0xFDFB | Row letter (A-H) | ✅ |
| 65021 | 0xFDFD | Droplet data pointer | ✅ |
| 65057 | 0xFE41 | Cluster array | ✅ |

### Major Unknowns ❓

1. **Why 2 channel name locations?**
   - Tag 65004: Always "Ch1,Ch2"
   - Unknown pattern: Actual names

2. **Droplet record fields**
   - 16 bytes unknown per droplet (out of 28 total)

3. **Quality metrics**
   - Multiple tags, meanings unclear

4. **Empty wells**
   - Files contain IFDs for ALL 96 wells
   - Most have cluster_size = 0 (no data)
   - Why include them?

---

## Limitations

1. **Only 2 channels** (no multiplex support yet)
2. **Channel name extraction is heuristic** (Unknown+32 pattern)
3. **No validation** (no checksums or error detection)
4. **Partial droplet record** (only Ch1/Ch2 amplitudes decoded)

---

## Contributing

Please test on:
- Different instruments (QX200/QX600)
- Different software versions
- Different assays

Report issues with:
- Software version & instrument model

---

## License

MIT

---

**REMEMBER**: This is **REVERSE-ENGINEERED**. It works on our files but **may not be universal**. QLP is proprietary and undocumented.