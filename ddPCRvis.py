#!/usr/bin/env python3
"""
ddPCRvis - ddPCR Visualization Tool (Enhanced Version v2 - FIXED)
A GUI application for visualizing ddPCR amplitude data from Bio-Rad CSV exports and QLP binary files

ENHANCED FEATURES (v2 - Over-segmentation Fix):
- 2D Gaussian Mixture Model clustering optimized for standard ddPCR
- FIXED: Now correctly identifies 2 bands per channel (negative + positive) for standard assays
- Shared thresholds across wells in the same experimental condition
- Percentile-based bounds (robust to outliers)
- Improved threshold placement using largest gap detection

KEY FIX IN V2:
- Standard ddPCR has 4 populations in 2D space: negative, FAM+, HEX+, double+
- This means 2 bands per channel (negative and positive), NOT 4
- Threshold placement based on largest gap between population means
- Secondary thresholds only added if there's another significant gap (>50% of primary)

CLUSTERING METHODS AVAILABLE:
1. detect_bands_1d_gmm() - 1D GMM, defaults to 2 bands (recommended for ddPCR)
2. detect_populations_2d_gmm() - 2D GMM with smart threshold selection (best for complex samples)
3. detect_bands_kmeans() - Legacy k-means (kept for backwards compatibility)

THRESHOLD MODES:
- Individual: Each well gets its own thresholds (default for CSV files)
- Shared: Wells in same condition share thresholds (default for QLP with assignments)
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import re
from datetime import datetime
import struct
from collections import defaultdict
import tempfile

try:
    from sklearn.mixture import GaussianMixture
    from scipy import stats
    from scipy.signal import find_peaks
    from scipy.stats import gaussian_kde
except ImportError:
    print("Installing scikit-learn and scipy...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-learn", "scipy", "--break-system-packages"])
    from sklearn.mixture import GaussianMixture
    from scipy import stats

try:
    from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                                 QHBoxLayout, QLabel, QPushButton, QTextEdit, 
                                 QFileDialog, QProgressBar, QDialog, QTableWidget,
                                 QTableWidgetItem, QHeaderView, QComboBox, QMessageBox,
                                 QScrollArea, QGridLayout, QGroupBox, QCheckBox)
    from PyQt5.QtCore import Qt, QThread, pyqtSignal
    from PyQt5.QtGui import QFont, QPalette, QColor, QDragEnterEvent, QDropEvent
except ImportError:
    print("PyQt5 not found. Installing...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "PyQt5", "--break-system-packages"])
    from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                                 QHBoxLayout, QLabel, QPushButton, QTextEdit, 
                                 QFileDialog, QProgressBar, QDialog, QTableWidget,
                                 QTableWidgetItem, QHeaderView, QComboBox, QMessageBox,
                                 QScrollArea, QGridLayout, QGroupBox, QCheckBox)
    from PyQt5.QtCore import Qt, QThread, pyqtSignal
    from PyQt5.QtGui import QFont, QPalette, QColor, QDragEnterEvent, QDropEvent


class BioRadQLPParser:
    """Minimal robust parser - reads well IDs directly from IFD data"""
    
    def __init__(self, filepath, debug=False):
        self.filepath = filepath
        self.debug = debug
        
        with open(filepath, 'rb') as f:
            self.data = f.read()
        
        self.endian = '<' if self.data[:2] == b'II' else '>'
        
        # Constants
        self.RECORD_SIZE = 28
        self.RECORD_FMT = "I fff fff"
        
        # Tag IDs
        self.TAG_WELL_NAME = 65019
        self.TAG_DATA_START = 65021
        self.TAG_CLUSTER_ARRAY = 65057
        self.TAG_QUALITY_ARRAY = 65054
        self.TAG_WELL_QUALITY = 65005
        self.TAG_CHANNEL_NAMES = 65004
        self.TAG_SOFTWARE = 305
        
        self.channel_names = None
        
        self.CLUSTER_BYTE_MAP = {
            0x00: 0, 0x11: 1, 0x22: 2, 0x33: 3, 0x44: 4
        }
        
        self.metadata = {}
        self.well_metadata = {}

    def _get_val(self, offset, fmt):
        try:
            return struct.unpack(f"{self.endian}{fmt}", 
                               self.data[offset:offset+struct.calcsize(fmt)])[0]
        except:
            return None
    
    def _extract_tags(self, ifd_offset):
        """Extract all tags from an IFD"""
        num_entries = self._get_val(ifd_offset, "H")
        tags = {}
        
        if num_entries is None:
            return tags
            
        for i in range(num_entries):
            entry_ptr = ifd_offset + 2 + (i * 12)
            tid = self._get_val(entry_ptr, "H")
            ttype = self._get_val(entry_ptr + 2, "H")
            count = self._get_val(entry_ptr + 4, "I")
            
            type_sizes = {1:1, 2:1, 3:2, 4:4, 5:8, 7:1, 9:4, 11:4, 12:8}
            item_size = type_sizes.get(ttype, 1)
            total_size = item_size * count
            
            val_or_offset = self._get_val(entry_ptr + 8, "I")
            ptr = entry_ptr + 8 if total_size <= 4 else val_or_offset
            
            tags[tid] = (ptr, total_size, ttype, count)
        
        return tags

    def _read_well_id_from_ifd(self, ifd_offset):
        """
        Read well ID directly from IFD data.
        Pattern: IFD starts with: 1c 00 fb fd 02 00 01 00 00 00 [WELL_ID]\x00
        Well ID is at offset +10 from IFD start (e.g., "C03", "D04")
        """
        # Read 4 bytes starting at offset +10
        well_id_offset = ifd_offset + 10
        well_id_bytes = self.data[well_id_offset:well_id_offset+4]
        
        # Parse as null-terminated string
        try:
            well_id = well_id_bytes.split(b'\x00')[0].decode('ascii')
            # Validate format: [A-H][0-9][0-9]
            if len(well_id) == 3 and well_id[0] in 'ABCDEFGH' and well_id[1:].isdigit():
                return well_id
        except:
            pass
        
        return None

    def _extract_metadata(self):
        """Extract file metadata"""
        ifd_offset = self._get_val(4, "I")
        if not ifd_offset:
            return

        tags = self._extract_tags(ifd_offset)
        
        if self.TAG_SOFTWARE in tags:
            ptr, size, _, _ = tags[self.TAG_SOFTWARE]
            if ptr and size:
                self.metadata['software'] = self.data[ptr:ptr+size].decode('ascii', errors='ignore').split('\x00')[0]
        
        if self.TAG_CHANNEL_NAMES in tags:
            ptr, size, _, _ = tags[self.TAG_CHANNEL_NAMES]
            if ptr and size:
                channel_str = self.data[ptr:ptr+size].decode('ascii', errors='ignore').split('\x00')[0]
                self.channel_names = [ch.strip() for ch in channel_str.split(',')]
        
        if not self.channel_names:
            self.channel_names = ['Ch1', 'Ch2']

    def parse_to_dataframe(self, include_filtered=True):
        """Parse QLP file and return dict of well_id -> DataFrame"""
        self._extract_metadata()
        
        ifd_offset = self._get_val(4, "I")
        well_data = defaultdict(list)

        while ifd_offset and ifd_offset < len(self.data):
            tags = self._extract_tags(ifd_offset)
            num_entries = self._get_val(ifd_offset, "H")
            
            if not num_entries:
                break

            # Check if this is a well block
            if (self.TAG_WELL_NAME in tags and 
                self.TAG_DATA_START in tags and 
                self.TAG_CLUSTER_ARRAY in tags):
                
                cluster_size = tags[self.TAG_CLUSTER_ARRAY][1]
                
                # Skip empty wells
                if cluster_size == 0:
                    ifd_offset = self._get_val(ifd_offset + 2 + (num_entries * 12), "I")
                    continue
                
                # Read well ID directly from IFD
                well_id = self._read_well_id_from_ifd(ifd_offset)
                
                if not well_id:
                    if self.debug:
                        print(f"Warning: Could not read well ID at IFD offset 0x{ifd_offset:08x}")
                    ifd_offset = self._get_val(ifd_offset + 2 + (num_entries * 12), "I")
                    continue
                
                if self.debug:
                    print(f"Well {well_id}: {cluster_size} droplets")
                
                # Extract well metadata
                well_meta = {}
                if self.TAG_WELL_QUALITY in tags:
                    ptr = tags[self.TAG_WELL_QUALITY][0]
                    well_meta['well_quality_flag'] = self._get_val(ptr, 'f')
                
                self.well_metadata[well_id] = well_meta
                
                # Read droplet data
                cluster_array = self.data[tags[self.TAG_CLUSTER_ARRAY][0]:
                                         tags[self.TAG_CLUSTER_ARRAY][0]+cluster_size]
                
                quality_array = None
                if self.TAG_QUALITY_ARRAY in tags:
                    q_ptr, q_size = tags[self.TAG_QUALITY_ARRAY][0], tags[self.TAG_QUALITY_ARRAY][1]
                    if q_size > 0:
                        quality_array = self.data[q_ptr:q_ptr+q_size]
                
                cursor = tags[self.TAG_DATA_START][0]
                droplet_idx = 0
                
                while cursor < len(self.data) - self.RECORD_SIZE and droplet_idx < len(cluster_array):
                    r = struct.unpack(f"{self.endian}{self.RECORD_FMT}", 
                                     self.data[cursor:cursor+self.RECORD_SIZE])
                    
                    cluster = self.CLUSTER_BYTE_MAP.get(cluster_array[droplet_idx], 0)
                    quality_flag = quality_array[droplet_idx] if quality_array and droplet_idx < len(quality_array) else 0
                    
                    if include_filtered or cluster > 0:
                        droplet = {
                            'DropletID': r[0],
                            'Ch1_Amplitude': r[1],
                            'Ch2_Amplitude': r[4],
                            'Cluster': cluster,
                            'Quality_Flag': quality_flag,
                        }
                        droplet.update(well_meta)
                        well_data[well_id].append(droplet)
                    
                    droplet_idx += 1
                    cursor += self.RECORD_SIZE

            # Move to next IFD
            ifd_offset = self._get_val(ifd_offset + 2 + (num_entries * 12), "I")

        # Convert to DataFrames
        well_dataframes = {}
        for well_id, data in well_data.items():
            if data:
                df = pd.DataFrame(data)
                df.attrs['channel_names'] = self.channel_names
                df.attrs['channel_map'] = {
                    'Ch1_Amplitude': self.channel_names[0],
                    'Ch2_Amplitude': self.channel_names[1]
                }
                well_dataframes[well_id] = df
        
        return well_dataframes

    def get_channel_names(self):
        if self.channel_names is None:
            self._extract_metadata()
        return self.channel_names or ['Ch1', 'Ch2']


class WellAssignmentDialog(QDialog):
    """Dialog for assigning wells to experimental conditions"""
    def __init__(self, well_ids, parent=None):
        super().__init__(parent)
        self.well_ids = sorted(well_ids)
        self.assignments = {}
        self.conditions = set()
        self.include_filtered = True  # Default: include cluster 0
        
        self.setWindowTitle("Assign Wells to Conditions")
        self.setModal(True)
        self.resize(900, 650)
        
        self.setup_ui()
        
    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Instructions
        instructions = QLabel(
            "Assign each well to an experimental condition.\n"
            "Wells with the same condition will be plotted together."
        )
        instructions.setFont(QFont("Arial", 10))
        instructions.setStyleSheet("color: #2c3e50; margin: 10px; padding: 10px; background-color: #ecf0f1; border-radius: 5px;")
        layout.addWidget(instructions)
        
        # Filter options
        filter_group = QGroupBox("Data Filtering Options")
        filter_layout = QVBoxLayout()
        
        self.include_filtered_checkbox = QCheckBox("Include filtered droplets (Cluster 0)")
        self.include_filtered_checkbox.setChecked(True)
        self.include_filtered_checkbox.setToolTip(
            "Cluster 0 droplets are filtered out by BioRad's software.\n"
            "Check this to include ALL droplets, uncheck to match BioRad's CSV export behavior."
        )
        self.include_filtered_checkbox.stateChanged.connect(self.on_filter_changed)
        
        filter_info = QLabel(
            "ℹ️ Cluster 0 droplets are excluded by BioRad when exporting to CSV.\n"
            "   Unchecking this option will match BioRad's default behavior."
        )
        filter_info.setStyleSheet("color: #7f8c8d; font-size: 9pt; margin-left: 20px;")
        
        filter_layout.addWidget(self.include_filtered_checkbox)
        filter_layout.addWidget(filter_info)
        filter_group.setLayout(filter_layout)
        layout.addWidget(filter_group)
        
        # Quick assign section
        quick_group = QGroupBox("Quick Assign")
        quick_layout = QVBoxLayout()
        
        quick_row = QHBoxLayout()
        quick_row.addWidget(QLabel("Condition name:"))
        self.quick_condition_input = QComboBox()
        self.quick_condition_input.setEditable(True)
        self.quick_condition_input.addItems(["Control", "Treatment", "Sample"])
        quick_row.addWidget(self.quick_condition_input)
        
        quick_row.addWidget(QLabel("Wells:"))
        self.quick_wells_input = QComboBox()
        self.quick_wells_input.setEditable(True)
        self.quick_wells_input.addItem("e.g., A01-A12, B01, C05-C08")
        quick_row.addWidget(self.quick_wells_input)
        
        assign_btn = QPushButton("Assign")
        assign_btn.clicked.connect(self.quick_assign)
        quick_row.addWidget(assign_btn)
        
        quick_layout.addLayout(quick_row)
        quick_group.setLayout(quick_layout)
        layout.addWidget(quick_group)
        
        # Scroll area for well grid
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll_widget = QWidget()
        
        # Create grid layout for wells (8 rows x 12 columns)
        grid = QGridLayout(scroll_widget)
        grid.setSpacing(5)
        
        # Add column headers
        for col in range(12):
            header = QLabel(f"{col+1:02d}")
            header.setAlignment(Qt.AlignCenter)
            header.setFont(QFont("Arial", 9, QFont.Bold))
            grid.addWidget(header, 0, col + 1)
        
        # Add row headers and well widgets
        self.well_widgets = {}
        for row in range(8):
            row_letter = chr(ord('A') + row)
            header = QLabel(row_letter)
            header.setAlignment(Qt.AlignCenter)
            header.setFont(QFont("Arial", 9, QFont.Bold))
            grid.addWidget(header, row + 1, 0)
            
            for col in range(12):
                well_id = f"{row_letter}{col+1:02d}"
                
                if well_id in self.well_ids:
                    well_widget = QWidget()
                    well_layout = QVBoxLayout(well_widget)
                    well_layout.setContentsMargins(2, 2, 2, 2)
                    
                    well_label = QLabel(well_id)
                    well_label.setAlignment(Qt.AlignCenter)
                    well_label.setFont(QFont("Arial", 8))
                    
                    condition_combo = QComboBox()
                    condition_combo.setEditable(True)
                    condition_combo.addItem("(unassigned)")
                    condition_combo.addItems(["Control", "Treatment", "Sample"])
                    condition_combo.currentTextChanged.connect(
                        lambda text, wid=well_id: self.update_assignment(wid, text)
                    )
                    
                    well_layout.addWidget(well_label)
                    well_layout.addWidget(condition_combo)
                    
                    well_widget.setStyleSheet("""
                        QWidget {
                            background-color: #e8f4f8;
                            border: 1px solid #3498db;
                            border-radius: 3px;
                        }
                    """)
                    
                    self.well_widgets[well_id] = condition_combo
                    grid.addWidget(well_widget, row + 1, col + 1)
                else:
                    # Empty well
                    empty = QLabel("")
                    empty.setStyleSheet("background-color: #ecf0f1;")
                    grid.addWidget(empty, row + 1, col + 1)
        
        scroll.setWidget(scroll_widget)
        layout.addWidget(scroll)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        clear_btn = QPushButton("Clear All")
        clear_btn.clicked.connect(self.clear_all)
        button_layout.addWidget(clear_btn)
        
        button_layout.addStretch()
        
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(cancel_btn)
        
        ok_btn = QPushButton("OK")
        ok_btn.setStyleSheet("""
            QPushButton {
                background-color: #27ae60;
                color: white;
                padding: 8px 20px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #229954;
            }
        """)
        ok_btn.clicked.connect(self.validate_and_accept)
        button_layout.addWidget(ok_btn)
        
        layout.addLayout(button_layout)
    
    def on_filter_changed(self, state):
        """Handle filter checkbox state change"""
        self.include_filtered = (state == Qt.Checked)
    
    def quick_assign(self):
        """Quick assign wells based on range input"""
        condition = self.quick_condition_input.currentText().strip()
        wells_text = self.quick_wells_input.currentText().strip()
        
        if not condition or condition == "(unassigned)":
            QMessageBox.warning(self, "Warning", "Please enter a condition name")
            return
        
        # Parse wells text
        wells_to_assign = self.parse_wells_range(wells_text)
        
        if not wells_to_assign:
            QMessageBox.warning(self, "Warning", "No valid wells found in range")
            return
        
        # Assign wells
        for well_id in wells_to_assign:
            if well_id in self.well_widgets:
                self.well_widgets[well_id].setCurrentText(condition)
        
        QMessageBox.information(self, "Success", f"Assigned {len(wells_to_assign)} wells to '{condition}'")
    
    def parse_wells_range(self, text):
        """Parse well range string like 'A01-A12, B01, C05-C08'"""
        wells = set()
        parts = [p.strip() for p in text.split(',')]
        
        for part in parts:
            if '-' in part:
                # Range
                try:
                    start, end = part.split('-')
                    start = start.strip()
                    end = end.strip()
                    
                    start_row = ord(start[0]) - ord('A')
                    start_col = int(start[1:])
                    end_row = ord(end[0]) - ord('A')
                    end_col = int(end[1:])
                    
                    if start_row == end_row:
                        # Same row
                        row_letter = start[0]
                        for col in range(start_col, end_col + 1):
                            wells.add(f"{row_letter}{col:02d}")
                    elif start_col == end_col:
                        # Same column
                        col = start_col
                        for row in range(start_row, end_row + 1):
                            wells.add(f"{chr(ord('A') + row)}{col:02d}")
                except:
                    pass
            else:
                # Single well
                if len(part) >= 2:
                    wells.add(part.upper())
        
        # Filter to only valid wells
        return [w for w in wells if w in self.well_ids]
    
    def update_assignment(self, well_id, condition):
        """Update assignment for a well"""
        if condition and condition != "(unassigned)":
            self.assignments[well_id] = condition
            self.conditions.add(condition)
        elif well_id in self.assignments:
            del self.assignments[well_id]
    
    def clear_all(self):
        """Clear all assignments"""
        for combo in self.well_widgets.values():
            combo.setCurrentIndex(0)
        self.assignments.clear()
        self.conditions.clear()
    
    def validate_and_accept(self):
        """Validate assignments before accepting"""
        if not self.assignments:
            reply = QMessageBox.question(
                self, "No Assignments",
                "No wells have been assigned to conditions. Continue anyway?",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
        
        self.accept()
    
    def get_assignments(self):
        """Return the assignments dict: {well_id: condition}"""
        return self.assignments
    
    def get_include_filtered(self):
        """Return whether to include cluster 0 droplets"""
        return self.include_filtered


class ProcessingThread(QThread):
    """Thread for processing ddPCR files without blocking the GUI"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)
    
    def __init__(self, files, file_type='csv', well_assignments=None, include_filtered=True, output_format='png'):
        super().__init__()
        self.files = files
        self.file_type = file_type
        self.well_assignments = well_assignments
        self.include_filtered = include_filtered
        self.output_format = output_format # Store format
        
    def run(self):
        try:
            if self.file_type == 'qlp':
                output_dir, stats_file = process_qlp_file(
                    self.files[0], 
                    self.well_assignments,
                    self.include_filtered,
                    self.progress.emit,
                    self.output_format # Pass to function
                )
            else:
                output_dir, stats_file = process_ddpcr_files(
                    self.files, 
                    self.progress.emit,
                    self.output_format # Pass to function
                )
            
            msg = f"Processing complete!\n\nOutputs saved to:\n{output_dir}"
            if stats_file:
                msg += f"\n\nStats saved to:\n{stats_file}"
            self.finished.emit(msg)
        except Exception as e:
            import traceback
            self.error.emit(f"Error processing files: {str(e)}\n\n{traceback.format_exc()}")


class DropZone(QLabel):
    """Custom label widget that accepts drag and drop"""
    files_dropped = pyqtSignal(list, str)  # files, file_type
    
    def __init__(self, text):
        super().__init__(text)
        self.setAcceptDrops(True)
        self.setAlignment(Qt.AlignCenter)
        self.setStyleSheet("""
            QLabel {
                border: 3px dashed #3498db;
                border-radius: 10px;
                background-color: #ecf0f1;
                padding: 40px;
                font-size: 16px;
                color: #2c3e50;
            }
            QLabel:hover {
                background-color: #d5dbdb;
                border-color: #2980b9;
            }
        """)
        
    def dragEnterEvent(self, event: QDragEnterEvent):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
            self.setStyleSheet("""
                QLabel {
                    border: 3px dashed #27ae60;
                    border-radius: 10px;
                    background-color: #d5f4e6;
                    padding: 40px;
                    font-size: 16px;
                    color: #2c3e50;
                }
            """)
        
    def dragLeaveEvent(self, event):
        self.setStyleSheet("""
            QLabel {
                border: 3px dashed #3498db;
                border-radius: 10px;
                background-color: #ecf0f1;
                padding: 40px;
                font-size: 16px;
                color: #2c3e50;
            }
        """)
        
    def dropEvent(self, event: QDropEvent):
        csv_files = []
        qlp_files = []
        
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            if file_path.endswith('.csv'):
                csv_files.append(file_path)
            elif file_path.endswith('.qlp'):
                qlp_files.append(file_path)
        
        if qlp_files:
            if len(qlp_files) > 1:
                QMessageBox.warning(None, "Warning", "Only one QLP file can be processed at a time. Using the first file.")
            self.files_dropped.emit([qlp_files[0]], 'qlp')
        elif csv_files:
            self.files_dropped.emit(csv_files, 'csv')
        
        self.setStyleSheet("""
            QLabel {
                border: 3px dashed #3498db;
                border-radius: 10px;
                background-color: #ecf0f1;
                padding: 40px;
                font-size: 16px;
                color: #2c3e50;
            }
        """)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ddPCRvis - ddPCR Visualization Tool")
        self.setGeometry(100, 100, 800, 600)
        
        # Create central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        
        # Title
        title = QLabel("ddPCR Visualization Tool")
        title.setFont(QFont("Arial", 24, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("color: #2c3e50; margin: 20px;")
        layout.addWidget(title)
        
        # Subtitle
        subtitle = QLabel("Drag and drop Bio-Rad CSV files or QLP binary files")
        subtitle.setFont(QFont("Arial", 12))
        subtitle.setAlignment(Qt.AlignCenter)
        subtitle.setStyleSheet("color: #7f8c8d; margin-bottom: 20px;")
        layout.addWidget(subtitle)
        
        # Drop zone
        self.drop_zone = DropZone("📁 Drop CSV or QLP files here\n\nor click 'Browse Files' below")
        self.drop_zone.files_dropped.connect(self.process_files)
        layout.addWidget(self.drop_zone)

        # Output Format Selection
        format_layout = QHBoxLayout()
        format_label = QLabel("Output Format:")
        format_label.setFont(QFont("Arial", 10, QFont.Bold))
        
        self.format_combo = QComboBox()
        self.format_combo.addItems(["PNG", "PDF (Vector)"])
        self.format_combo.setToolTip("PNG is recommended for large datasets (faster to save and view).\nPDF preserves vector data but files can be huge.")
        self.format_combo.setStyleSheet("padding: 5px;")
        
        format_layout.addStretch()
        format_layout.addWidget(format_label)
        format_layout.addWidget(self.format_combo)
        format_layout.addStretch()
        layout.addLayout(format_layout)
        
        # Browse buttons
        browse_layout = QHBoxLayout()
        
        browse_csv_btn = QPushButton("📂 Browse CSV Files")
        browse_csv_btn.setFont(QFont("Arial", 11))
        browse_csv_btn.setStyleSheet("""
            QPushButton {
                background-color: #3498db;
                color: white;
                border: none;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #2980b9;
            }
        """)
        browse_csv_btn.clicked.connect(self.browse_csv_files)
        
        browse_qlp_btn = QPushButton("📂 Browse QLP File")
        browse_qlp_btn.setFont(QFont("Arial", 11))
        browse_qlp_btn.setStyleSheet("""
            QPushButton {
                background-color: #9b59b6;
                color: white;
                border: none;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #8e44ad;
            }
        """)
        browse_qlp_btn.clicked.connect(self.browse_qlp_file)
        
        browse_layout.addStretch()
        browse_layout.addWidget(browse_csv_btn)
        browse_layout.addWidget(browse_qlp_btn)
        browse_layout.addStretch()
        layout.addLayout(browse_layout)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 2px solid #bdc3c7;
                border-radius: 5px;
                text-align: center;
                height: 25px;
            }
            QProgressBar::chunk {
                background-color: #27ae60;
            }
        """)
        layout.addWidget(self.progress_bar)
        
        # Status/log area
        self.log_area = QTextEdit()
        self.log_area.setReadOnly(True)
        self.log_area.setFont(QFont("Courier", 9))
        self.log_area.setStyleSheet("""
            QTextEdit {
                background-color: #2c3e50;
                color: #ecf0f1;
                border: 1px solid #34495e;
                border-radius: 5px;
                padding: 10px;
            }
        """)
        layout.addWidget(self.log_area)
        
        # Initial log message
        self.log("Welcome to ddPCRvis!")
        self.log("Drop CSV files, QLP files, or click 'Browse Files' to get started.")
        self.log("- CSV files: Will be plotted individually and combined")
        self.log("- QLP files: Will prompt for well-to-condition assignments")
        
    def log(self, message):
        """Add a message to the log area"""
        self.log_area.append(message)
        
    def browse_csv_files(self):
        """Open file dialog to select CSV files"""
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Select CSV Files",
            "",
            "CSV Files (*.csv)"
        )
        if files:
            self.process_files(files, 'csv')
    
    def browse_qlp_file(self):
        """Open file dialog to select a QLP file"""
        file, _ = QFileDialog.getOpenFileName(
            self,
            "Select QLP File",
            "",
            "QLP Files (*.qlp)"
        )
        if file:
            self.process_files([file], 'qlp')
    
    def process_files(self, files, file_type):
        """Process the selected files"""
        if not files:
            return
        
        # Determine format
        is_png = self.format_combo.currentIndex() == 0
        output_format = 'png' if is_png else 'pdf'
        
        self.log(f"\n{'='*60}")
        self.log(f"Processing {len(files)} {file_type.upper()} file(s)...")
        
        well_assignments = None
        include_filtered = True  # Default for CSV files
        
        if file_type == 'qlp':
            # Parse QLP file to get well IDs (quick parse to get well names only)
            try:
                self.log("Parsing QLP file...")
                parser = BioRadQLPParser(files[0])
                # Quick parse just to get well IDs
                well_dataframes = parser.parse_to_dataframe(include_filtered=True)
                well_ids = list(well_dataframes.keys())
                
                self.log(f"Found {len(well_ids)} wells: {', '.join(sorted(well_ids))}")
                
                # Show well assignment dialog
                dialog = WellAssignmentDialog(well_ids, self)
                if dialog.exec_() == QDialog.Accepted:
                    well_assignments = dialog.get_assignments()
                    include_filtered = dialog.get_include_filtered()
                    self.log(f"Well assignments: {well_assignments}")
                    self.log(f"Include filtered droplets (Cluster 0): {include_filtered}")
                else:
                    self.log("Processing cancelled by user")
                    return
                    
            except Exception as e:
                self.log(f"Error parsing QLP file: {str(e)}")
                QMessageBox.critical(self, "Error", f"Failed to parse QLP file:\n{str(e)}")
                return
        
        # Start processing in background thread
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        # UPDATE: Pass output_format to Thread
        self.thread = ProcessingThread(files, file_type, well_assignments, include_filtered, output_format)
        self.thread.progress.connect(self.log)
        self.thread.finished.connect(self.on_processing_finished)
        self.thread.error.connect(self.on_processing_error)
        self.thread.start()
    
    def on_processing_finished(self, message):
        """Handle successful processing completion"""
        self.progress_bar.setVisible(False)
        self.log("\n" + message)
        QMessageBox.information(self, "Success", message)
    
    def on_processing_error(self, error_message):
        """Handle processing error"""
        self.progress_bar.setVisible(False)
        self.log(f"\n❌ ERROR: {error_message}")
        QMessageBox.critical(self, "Error", error_message)


def format_population_stats(count, total):
    """
    Format count and percentage with scientific notation for rare events.
    Format: "Count (Percentage%)"
    Example: "1250 (12.5%)" or "3 (2.5e-3%)"
    """
    if total == 0:
        return "0 (0%)"
    
    pct = (count / total) * 100
    
    if count == 0:
        pct_str = "0%"
    elif pct < 0.01:
        # Scientific notation for very small percentages
        # e.g. 1.5e-04%
        pct_str = f"{pct:.1e}%"
    else:
        # Standard fixed point for normal percentages
        pct_str = f"{pct:.1f}%"
        
    return f"{count} ({pct_str})"

def get_density_based_thresholds(values, min_separation=500, rain_sensitivity=0.1):
    """
    Revised V3: Log-Space Peak Detection.
    
    Fixes the issue where dense negative populations drown out distinct 
    but smaller positive populations (middle/top bands).
    
    1. Uses Log-Density for peak finding (compresses dynamic range).
    2. Maintains strict safety buffers from the previous version.
    3. Detects multiple bands even if counts vary by 100x.
    """
    values = np.array(values)
    if len(values) < 10:
        return []
    
    val_range = values.max() - values.min()
    if val_range < min_separation:
        return [] 

    # 1. Estimate Density
    grid_points = np.linspace(values.min(), values.max(), 500)
    try:
        kde = gaussian_kde(values, bw_method='scott')
        density = kde(grid_points)
    except:
        return []

    # 2. Log-Transformation (The Fix)
    # This makes small populations visible even if the negative peak is massive
    # epsilon prevents log(0)
    epsilon = 1e-10
    log_density = np.log(density + epsilon)
    
    # 3. Find Peaks in Log-Space
    # Prominence is now relative to the log-dynamic range. 
    # 0.05 (5%) of the log range is usually robust enough to ignore noise 
    # but catch real populations.
    log_range = log_density.max() - log_density.min()
    peaks, properties = find_peaks(log_density, prominence=0.05 * log_range)
    
    # Calculate robust stats for safety checks (based on main negative cluster)
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    sigma_est = 1.4826 * mad
    
    thresholds = []
    
    # Case A: Multiple populations (Valley detection)
    if len(peaks) >= 2:
        for i in range(len(peaks) - 1):
            p1_idx = peaks[i]
            p2_idx = peaks[i+1]
            
            # Find valley in the ORIGINAL linear density (better for finding true minimum)
            # We use the indices found from log-space, but search the linear density for the low point
            valley_slice = density[p1_idx:p2_idx]
            local_min_idx = np.argmin(valley_slice) + p1_idx
            raw_threshold = grid_points[local_min_idx]
            
            # SAFETY CHECK: 
            # Ensure threshold is at least 4 sigmas away from the left peak
            peak_val = grid_points[p1_idx]
            min_safe_dist = peak_val + (4 * sigma_est)
            
            final_threshold = max(raw_threshold, min_safe_dist)
            thresholds.append(final_threshold)

    # Case B: Single Peak (Sparse detection)
    elif len(peaks) == 1:
        peak_val = grid_points[peaks[0]]
        
        # Strict sparse detection (12 sigma)
        cutoff = peak_val + max(12 * sigma_est, min_separation)
        
        if values.max() > cutoff:
            thresholds.append(cutoff)

    # 4. Filter thresholds based on minimum separation
    valid_thresholds = []
    thresholds = sorted(thresholds)
    
    if thresholds:
        peak_values = grid_points[peaks]
        for thresh in thresholds:
            # Check distance to nearest peak
            dist_to_peak = min(abs(thresh - p) for p in peak_values)
            # Check distance to existing valid thresholds
            dist_to_thresh = min([abs(thresh - t) for t in valid_thresholds]) if valid_thresholds else float('inf')
            
            if dist_to_peak > (min_separation / 3) and dist_to_thresh > min_separation:
                valid_thresholds.append(thresh)

    return valid_thresholds


def plot_well_2d_scatter(df, well_id, ch1_bands, ch2_bands):
    """
    Create 2D scatter plot for Ch1 vs Ch2 with detailed stats in legend.
    """
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    
    channel_map = df.attrs.get('channel_map', {'Ch1_Amplitude': 'Ch1', 'Ch2_Amplitude': 'Ch2'})
    ch1_name = channel_map.get('Ch1_Amplitude', 'Ch1')
    ch2_name = channel_map.get('Ch2_Amplitude', 'Ch2')
    
    ch1_colors = create_gradient_colors('#0000FF', len(ch1_bands))
    ch2_colors = create_gradient_colors('#00FF00', len(ch2_bands))
    
    total_droplets = len(df)
    
    # Plot each combination of Ch1 and Ch2 bands
    for ch1_idx, ch1_band in enumerate(ch1_bands):
        ch1_mask = (df['Ch1_Amplitude'] >= ch1_band['min']) & (df['Ch1_Amplitude'] <= ch1_band['max'])
        
        for ch2_idx, ch2_band in enumerate(ch2_bands):
            ch2_mask = (df['Ch2_Amplitude'] >= ch2_band['min']) & (df['Ch2_Amplitude'] <= ch2_band['max'])
            
            combined_mask = ch1_mask & ch2_mask
            band_data = df[combined_mask]
            
            count = len(band_data)
            # Always plot/label if valid combination, even if count is 0 (though loop logic usually prevents 0)
            if count >= 0: 
                # Determine Color
                color = ch1_colors[ch1_idx]
                alpha = 0.3 + (0.4 * ch2_idx / max(len(ch2_bands) - 1, 1))
                
                # Format Label with Stats
                stats_str = format_population_stats(count, total_droplets)
                
                # Determine Population Name for Legend
                # Logic: 0=Neg, >0=Pos
                if ch1_idx == 0 and ch2_idx == 0:
                    pop_name = "Neg/Neg"
                elif ch1_idx > 0 and ch2_idx == 0:
                    pop_name = f"{ch1_name}+"
                elif ch1_idx == 0 and ch2_idx > 0:
                    pop_name = f"{ch2_name}+"
                elif ch1_idx > 0 and ch2_idx > 0:
                    pop_name = "Double+"
                else:
                    pop_name = f"B{ch1_idx}/B{ch2_idx}"
                
                label_text = f"{pop_name}: {stats_str}"

                if count > 0:
                    ax.scatter(band_data['Ch1_Amplitude'], 
                              band_data['Ch2_Amplitude'],
                              c=color, alpha=alpha, s=5, edgecolors='none',
                              label=label_text)
                else:
                    # Plot ghost point for legend entry
                    ax.scatter([], [], c=color, label=label_text)

    # Draw separators
    for ch1_idx in range(len(ch1_bands) - 1):
        # Use existing thresholds if available in band dict, else midpoint
        if ch1_bands[ch1_idx].get('threshold_upper'):
            sep = ch1_bands[ch1_idx]['threshold_upper']
        else:
            sep = (ch1_bands[ch1_idx]['max'] + ch1_bands[ch1_idx + 1]['min']) / 2
        ax.axvline(x=sep, color='red', linestyle='--', linewidth=1, alpha=0.5)
    
    for ch2_idx in range(len(ch2_bands) - 1):
        if ch2_bands[ch2_idx].get('threshold_upper'):
            sep = ch2_bands[ch2_idx]['threshold_upper']
        else:
            sep = (ch2_bands[ch2_idx]['max'] + ch2_bands[ch2_idx + 1]['min']) / 2
        ax.axhline(y=sep, color='red', linestyle='--', linewidth=1, alpha=0.5)
    
    ax.set_xlabel(f'{ch1_name} Amplitude', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'{ch2_name} Amplitude', fontsize=12, fontweight='bold')
    ax.set_title(f'Well {well_id} - 2D Scatter ({ch1_name} vs {ch2_name})', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.2)
    
    # Place Legend outside if too many items, or best location
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=10)
    
    plt.tight_layout()
    return fig


def plot_well_with_bands(df, well_id, ch1_bands, ch2_bands, limits=None):
    """
    Create band-colored 1D amplitude plot with proper axis limits and threshold lines.
    UPDATED: Legends now use scientific notation for percentages via format_population_stats.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Get actual channel names
    channel_map = df.attrs.get('channel_map', {'Ch1_Amplitude': 'Ch1', 'Ch2_Amplitude': 'Ch2'})
    ch1_name = channel_map.get('Ch1_Amplitude', 'Ch1')
    ch2_name = channel_map.get('Ch2_Amplitude', 'Ch2')
    
    # Generate gradient colors
    ch1_colors = create_gradient_colors('#0000FF', len(ch1_bands))
    ch2_colors = create_gradient_colors('#00FF00', len(ch2_bands))
    
    total_droplets = len(df)
    
    # Set Limits
    if limits and 'ch1_limits' in limits:
        ylim_ch1 = limits['ch1_limits']
        ylim_ch2 = limits['ch2_limits']
    else:
        # Fallback to local scaling
        ch1_all = df['Ch1_Amplitude'].values
        ch2_all = df['Ch2_Amplitude'].values
        ylim_ch1 = (ch1_all.min() * 0.95, ch1_all.max() * 1.05) if len(ch1_all) > 0 else (0, 10000)
        ylim_ch2 = (ch2_all.min() * 0.95, ch2_all.max() * 1.05) if len(ch2_all) > 0 else (0, 10000)
    
    # --- Plot Ch1 with bands ---
    for band_idx, band in enumerate(ch1_bands):
        mask = (df['Ch1_Amplitude'] >= band['min']) & (df['Ch1_Amplitude'] <= band['max'])
        band_data = df[mask]
        
        if len(band_data) > 0:
            x_jitter = np.random.uniform(-0.3, 0.3, len(band_data))
            
            # Format label using the helper function
            stats_str = format_population_stats(band["count"], total_droplets)
            label_text = f'Pop {band_idx+1}: {stats_str}'
            
            ax1.scatter(x_jitter, band_data['Ch1_Amplitude'], 
                       c=ch1_colors[band_idx], alpha=0.6, s=10,
                       label=label_text)
            
            # Draw threshold line
            if band['threshold_upper'] is not None:
                ax1.axhline(y=band['threshold_upper'], color='red', linestyle='--', 
                           linewidth=2, alpha=0.8)
    
    ax1.set_ylabel(f'{ch1_name} Amplitude', fontsize=12, fontweight='bold')
    ax1.set_xlim(-0.5, 0.5)
    ax1.set_ylim(ylim_ch1)
    ax1.set_xticks([])
    ax1.set_title(f'{ch1_name} ({len(ch1_bands)} populations)', fontsize=14, fontweight='bold')
    if ch1_bands:
        ax1.legend(loc='best', fontsize=9, framealpha=0.9)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # --- Plot Ch2 with bands ---
    for band_idx, band in enumerate(ch2_bands):
        mask = (df['Ch2_Amplitude'] >= band['min']) & (df['Ch2_Amplitude'] <= band['max'])
        band_data = df[mask]
        
        if len(band_data) > 0:
            x_jitter = np.random.uniform(-0.3, 0.3, len(band_data))
            
            # Format label using the helper function
            stats_str = format_population_stats(band["count"], total_droplets)
            label_text = f'Pop {band_idx+1}: {stats_str}'
            
            ax2.scatter(x_jitter, band_data['Ch2_Amplitude'],
                       c=ch2_colors[band_idx], alpha=0.6, s=10,
                       label=label_text)
            
            # Draw threshold line
            if band['threshold_upper'] is not None:
                ax2.axhline(y=band['threshold_upper'], color='red', linestyle='--', 
                           linewidth=2, alpha=0.8)
    
    ax2.set_ylabel(f'{ch2_name} Amplitude', fontsize=12, fontweight='bold')
    ax2.set_xlim(-0.5, 0.5)
    ax2.set_ylim(ylim_ch2)
    ax2.set_xticks([])
    ax2.set_title(f'{ch2_name} ({len(ch2_bands)} populations)', fontsize=14, fontweight='bold')
    if ch2_bands:
        ax2.legend(loc='best', fontsize=9, framealpha=0.9)
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle(f'Well {well_id} - Population Analysis\n' + 
                 f'Total droplets: {total_droplets:,}', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    return fig


def plot_1d_amplitude(df, title="ddPCR Amplitude Plot"):
    """Create 1D amplitude plot for a single sample (legacy, kept for compatibility)"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Cluster colors
    cluster_colors = {1: '#808080', 2: '#0000FF', 3: '#00FF00', 4: '#FF0000'}
    cluster_names = {1: 'Negative', 2: 'Ch1+', 3: 'Ch2+', 4: 'Double+'}
    
    # Plot Ch1
    for cluster in [1, 2, 3, 4]:
        mask = df['Cluster'] == cluster
        if mask.any():
            cluster_data = df[mask]
            x_jitter = np.random.uniform(-0.3, 0.3, len(cluster_data))
            ax1.scatter(x_jitter, cluster_data['Ch1_Amplitude'], 
                       c=cluster_colors[cluster], alpha=0.5, s=10,
                       label=cluster_names[cluster])
    
    ax1.set_ylabel('Ch1 Amplitude (FAM)', fontsize=12)
    ax1.set_xlim(-0.5, 0.5)
    ax1.set_xticks([])
    ax1.set_title('Channel 1 (FAM)', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    # Plot Ch2
    for cluster in [1, 2, 3, 4]:
        mask = df['Cluster'] == cluster
        if mask.any():
            cluster_data = df[mask]
            x_jitter = np.random.uniform(-0.3, 0.3, len(cluster_data))
            ax2.scatter(x_jitter, cluster_data['Ch2_Amplitude'],
                       c=cluster_colors[cluster], alpha=0.5, s=10,
                       label=cluster_names[cluster])
    
    ax2.set_ylabel('Ch2 Amplitude (HEX)', fontsize=12)
    ax2.set_xlim(-0.5, 0.5)
    ax2.set_xticks([])
    ax2.set_title('Channel 2 (HEX)', fontsize=14, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    plt.suptitle(title, fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    return fig


def create_gradient_colors(base_color_hex, n_bands):
    """Create gradient colors from dark (bottom/low amplitude) to light (top/high amplitude)
    
    Since bands are sorted by center value (ascending), band index 0 is the lowest/darkest,
    and the last band index is the highest/lightest.
    """
    base_r = int(base_color_hex[1:3], 16)
    base_g = int(base_color_hex[3:5], 16)
    base_b = int(base_color_hex[5:7], 16)
    
    colors = []
    for i in range(n_bands):
        # i=0 is darkest (bottom/lowest amplitude), i=n_bands-1 is lightest (top/highest amplitude)
        factor = 0.3 + (0.7 * i / max(n_bands - 1, 1))
        r = int(base_r * factor)
        g = int(base_g * factor)
        b = int(base_b * factor)
        colors.append(f'#{r:02x}{g:02x}{b:02x}')
    
    return colors




def detect_populations_2d_gmm_robust(ch1_values, ch2_values, 
                                      max_populations=4,
                                      min_population_size=10,
                                      min_population_fraction=0.001):
    """
    Refined logic: Uses 1D density thresholding on Ch1 and Ch2 independently,
    then combines them to form the 2D populations (Neg, Ch1+, Ch2+, Double+).
    This handles 'rain' and sparse positives much better than GMM.
    """
    n_droplets = len(ch1_values)
    
    # Get thresholds using the new robust method
    ch1_thresholds = get_density_based_thresholds(ch1_values, min_separation=500)
    ch2_thresholds = get_density_based_thresholds(ch2_values, min_separation=500)
    
    # Assign labels based on grid logic (Standard ddPCR approach)
    # 0=Neg, 1=Ch1+, 2=Ch2+, 3=Double+ (Simplified logic for 1 threshold per channel)
    # If multiple thresholds exist (multiplexing), logic adapts.
    
    labels = np.zeros(n_droplets, dtype=int)
    
    # Simple case: 1 threshold per channel (Standard)
    t1 = ch1_thresholds[0] if ch1_thresholds else float('inf')
    t2 = ch2_thresholds[0] if ch2_thresholds else float('inf')
    
    # Create mask for standard 4-cluster output
    is_ch1_pos = ch1_values > t1
    is_ch2_pos = ch2_values > t2
    
    # Vectorized label assignment
    # 0: Neg/Neg
    # 1: Pos/Neg (Ch1 only)
    # 2: Neg/Pos (Ch2 only)
    # 3: Pos/Pos (Double)
    labels = (is_ch1_pos.astype(int) * 1) + (is_ch2_pos.astype(int) * 2)
    
    # Remap to match specific request of user visualization if needed, 
    # but standardizing on: 0=Neg, 1=Ch1+, 2=Ch2+, 3=Dbl+ is safest.
    
    # Calculate stats for the populations found
    unique_labels = np.unique(labels)
    populations = []
    
    for label in unique_labels:
        mask = labels == label
        count = np.sum(mask)
        if count > 0:
            populations.append({
                'id': int(label),
                'center_ch1': float(np.mean(ch1_values[mask])),
                'center_ch2': float(np.mean(ch2_values[mask])),
                'count': int(count),
                'proportion': float(count / n_droplets * 100),
                'label': int(label)
            })

    return {
        'n_populations': len(unique_labels),
        'populations': populations,
        'ch1_thresholds': ch1_thresholds,
        'ch2_thresholds': ch2_thresholds,
        'labels': labels
    }




# Compatibility alias for old function name
def detect_populations_2d_gmm(ch1_values, ch2_values, max_populations=4, min_population_size=50):
    """Compatibility wrapper for old function name - calls new robust version"""
    return detect_populations_2d_gmm_robust(
        ch1_values, ch2_values,
        max_populations=max_populations,
        min_population_size=min_population_size,
        min_population_fraction=0.001
    )


def detect_bands_1d_gmm(values, max_bands=4, percentile_bounds=True, force_two_bands=False):
    """
    Detect bands in 1D - DEPRECATED in favor of 2D approach.
    Kept for backwards compatibility only.
    
    For ddPCR, use detect_populations_2d_gmm_robust() instead, then apply_thresholds_to_get_bands().
    """
    # This function is now a thin wrapper - real work done in 2D
    if len(values) < 10:
        return []
    
    # Simple percentile-based approach for 1D-only case
    values_array = np.array(values).reshape(-1, 1)
    
    try:
        # Try 2-component model
        gmm = GaussianMixture(n_components=2, covariance_type='full', random_state=42, n_init=10)
        gmm.fit(values_array)
        
        labels = gmm.predict(values_array)
        centers = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())
        
        # Check if well separated
        sorted_centers = np.sort(centers)
        gap = sorted_centers[1] - sorted_centers[0]
        avg_std = np.mean(stds)
        
        if gap < 2 * avg_std:
            # Fall back to 1 component
            gmm = GaussianMixture(n_components=1, covariance_type='full', random_state=42, n_init=10)
            gmm.fit(values_array)
            labels = gmm.predict(values_array)
            centers = gmm.means_.flatten()
            stds = np.sqrt(gmm.covariances_.flatten())
    
    except Exception:
        gmm = GaussianMixture(n_components=1, covariance_type='full', random_state=42, n_init=10)
        gmm.fit(values_array)
        labels = gmm.predict(values_array)
        centers = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())
    
    # Sort by center and create bands
    indices = np.argsort(centers)
    
    bands = []
    for idx, i in enumerate(indices):
        mask = labels == i
        band_values = values[mask]
        
        if len(band_values) > 0:
            if percentile_bounds:
                band_min = np.percentile(band_values, 0.5)
                band_max = np.percentile(band_values, 99.5)
            else:
                band_min = float(band_values.min())
                band_max = float(band_values.max())
            
            band = {
                'min': band_min,
                'max': band_max,
                'center': float(centers[i]),
                'mean': float(np.mean(band_values)),
                'std': float(stds[i]),
                'count': len(band_values),
                'proportion': len(band_values) / len(values) * 100,
                'threshold_lower': None,
                'threshold_upper': None
            }
            
            if idx > 0:
                prev_center = bands[-1]['center']
                band['threshold_lower'] = (prev_center + band['center']) / 2
                bands[-1]['threshold_upper'] = band['threshold_lower']
            
            bands.append(band)
    
    return bands


def calculate_shared_thresholds(well_data_dict, condition_assignments=None, method='density'):
    """
    Calculate shared thresholds AND Axis Limits across wells in the same condition.
    """
    if condition_assignments is None:
        condition_assignments = {well: 'All' for well in well_data_dict.keys()}
    
    # Group wells by condition
    condition_wells = defaultdict(list)
    for well, condition in condition_assignments.items():
        condition_wells[condition].append(well)
    
    shared_data = {}
    
    for condition, wells in condition_wells.items():
        # Pool all droplets from wells in this condition
        all_ch1 = []
        all_ch2 = []
        
        for well in wells:
            if well in well_data_dict:
                df = well_data_dict[well]
                all_ch1.extend(df['Ch1_Amplitude'].values)
                all_ch2.extend(df['Ch2_Amplitude'].values)
        
        if len(all_ch1) == 0:
            continue
            
        all_ch1 = np.array(all_ch1)
        all_ch2 = np.array(all_ch2)

        # 1. Calculate Thresholds using new Robust Density method
        ch1_thresh = get_density_based_thresholds(all_ch1, min_separation=500)
        ch2_thresh = get_density_based_thresholds(all_ch2, min_separation=500)
        
        # 2. Calculate Shared Axis Limits (Max/Min with padding)
        # This ensures all plots in the condition look identical in scale
        c1_min, c1_max = all_ch1.min(), all_ch1.max()
        c2_min, c2_max = all_ch2.min(), all_ch2.max()
        
        pad_c1 = (c1_max - c1_min) * 0.05 if (c1_max - c1_min) > 0 else 100
        pad_c2 = (c2_max - c2_min) * 0.05 if (c2_max - c2_min) > 0 else 100
        
        shared_data[condition] = {
            'ch1_thresholds': ch1_thresh,
            'ch2_thresholds': ch2_thresh,
            'ch1_limits': (c1_min - pad_c1, c1_max + pad_c1),
            'ch2_limits': (c2_min - pad_c2, c2_max + pad_c2)
        }
    
    return shared_data


def apply_thresholds_to_get_bands(values, thresholds, percentile_bounds=False):
    """
    Updated: Adds 'band_index' to allow consistent coloring across wells.
    """
    if len(values) == 0:
        return []
    
    values = np.array(values)
    thresholds = sorted([t for t in thresholds if np.isfinite(t)])
    bins = [-np.inf] + thresholds + [np.inf]
    
    bands = []
    # Loop through ALL defined bins (Logical Bands)
    for i in range(len(bins) - 1):
        lower = bins[i]
        upper = bins[i + 1]
        
        mask = (values > lower) & (values <= upper)
        band_values = values[mask]
        
        if len(band_values) > 0:
            band_min = float(band_values.min())
            band_max = float(band_values.max())
            
            bands.append({
                'band_index': i, # CRITICAL: Stores which tier this band belongs to (0=Lowest)
                'min': band_min,
                'max': band_max,
                'center': float(np.median(band_values)),
                'mean': float(np.mean(band_values)),
                'std': float(np.std(band_values)),
                'count': len(band_values),
                'proportion': len(band_values) / len(values) * 100,
                'threshold_lower': float(lower) if not np.isinf(lower) else None,
                'threshold_upper': float(upper) if not np.isinf(upper) else None
            })
    
    return bands


# Keep old function for backwards compatibility
def detect_bands_kmeans(values, max_bands=4, max_sample_size=10000):
    """Legacy k-means function - DEPRECATED. Use detect_populations_2d_gmm_robust instead."""
    from sklearn.cluster import KMeans
    
    if len(values) < 10:
        return []
    
    # Simplified version - just use GMM instead
    return detect_bands_1d_gmm(values, max_bands=max_bands)



def plot_96well_layout(data_dict, output_base_path, log_func=print, 
                       condition_assignments=None, use_shared_thresholds=True, 
                       output_format='png'):
    """
    Create 96-well plate visualization.
    FIXED: Uses robust density thresholding for individual wells (replacing legacy GMM)
    to ensure 'band_index' exists and strict thresholding is applied everywhere.
    """
    try:
        from sklearn.cluster import KMeans
    except ImportError:
        pass # imports handled at top of file
    
    log_func("\nCreating 96-well plate visualizations...")
    
    # Calculate shared thresholds if requested
    shared_thresh = None
    if use_shared_thresholds and condition_assignments:
        log_func("  Calculating shared thresholds per condition...")
        shared_thresh = calculate_shared_thresholds(
            data_dict,
            condition_assignments,
            method='density' # Use the new method
        )
        
        for condition, thresholds in shared_thresh.items():
            ch1_str = ', '.join([f'{t:.1f}' for t in thresholds['ch1_thresholds']])
            ch2_str = ', '.join([f'{t:.1f}' for t in thresholds['ch2_thresholds']])
            log_func(f"    {condition}: Ch1=[{ch1_str}], Ch2=[{ch2_str}]")
    
    well_info = {}
    for filename, df in data_dict.items():
        well_match = re.search(r'([A-H])(\d{1,2})', filename)
        if well_match:
            row_letter = well_match.group(1)
            col_num = int(well_match.group(2))
            well_id = f'{row_letter}{col_num:02d}'
            
            current_ch1_thresh = []
            current_ch2_thresh = []
            limits = None
            
            # LOGIC 1: Use Shared Thresholds (Grouped)
            if shared_thresh and condition_assignments:
                condition = condition_assignments.get(well_id, 'All')
                if condition in shared_thresh:
                    current_ch1_thresh = shared_thresh[condition]['ch1_thresholds']
                    current_ch2_thresh = shared_thresh[condition]['ch2_thresholds']
                    # Apply shared limits
                    limits = {
                        'ch1_limits': shared_thresh[condition]['ch1_limits'],
                        'ch2_limits': shared_thresh[condition]['ch2_limits']
                    }
            
            # LOGIC 2: Individual Analysis (Fallback)
            # Use new robust functions to ensure consistency (band_index, strictness)
            if not current_ch1_thresh and not current_ch2_thresh:
                current_ch1_thresh = get_density_based_thresholds(df['Ch1_Amplitude'].values, min_separation=500)
                current_ch2_thresh = get_density_based_thresholds(df['Ch2_Amplitude'].values, min_separation=500)
            
            # Calculate Bands using the thresholds
            ch1_bands = apply_thresholds_to_get_bands(
                df['Ch1_Amplitude'].values,
                current_ch1_thresh,
                percentile_bounds=False # Use False to avoid truncation
            )
            ch2_bands = apply_thresholds_to_get_bands(
                df['Ch2_Amplitude'].values,
                current_ch2_thresh,
                percentile_bounds=False
            )
            
            well_info[well_id] = {
                'df': df,
                'ch1_bands': ch1_bands,
                'ch2_bands': ch2_bands,
                'ch1_thresholds': current_ch1_thresh,
                'ch2_thresholds': current_ch2_thresh,
                'condition': condition_assignments.get(well_id, 'All') if condition_assignments else 'All',
                'limits': limits
            }
    
    if not well_info:
        log_func("  ⚠️  No valid well identifiers found in filenames")
        return
    
    # PLOTTING LOGIC
    output_base_path = Path(output_base_path)
    
    if output_format == 'pdf':
        # PDF Mode: Use PdfPages for a single file with multiple pages
        pdf_path = output_base_path.with_suffix('.pdf')
        with PdfPages(pdf_path) as pdf:
            create_96well_1d_page(well_info, pdf, log_func)
            create_96well_2d_page(well_info, pdf, log_func)
        log_func(f"  ✓ 96-well plate saved to: {pdf_path.name}")
        
    else:
        # PNG Mode: Save separate high-res images
        path_1d = output_base_path.parent / f"{output_base_path.stem}_1D_view.png"
        path_2d = output_base_path.parent / f"{output_base_path.stem}_2D_view.png"
        
        # Pass the PATH string instead of a PdfPages object
        create_96well_1d_page(well_info, path_1d, log_func)
        create_96well_2d_page(well_info, path_2d, log_func)
        
        log_func(f"  ✓ 96-well plate saved to: {path_1d.name} and {path_2d.name}")


def create_96well_1d_page(well_info, destination, log_func):
    """
    Updated: Lists Count + % for all bands in the 1D view.
    """
    all_ch1 = [v for info in well_info.values() for v in info['df']['Ch1_Amplitude']]
    all_ch2 = [v for info in well_info.values() for v in info['df']['Ch2_Amplitude']]
    ch1_min_global, ch1_max_global = min(all_ch1), max(all_ch1)
    ch2_min_global, ch2_max_global = min(all_ch2), max(all_ch2)
    
    fig = plt.figure(figsize=(20, 12))
    
    for row in range(8):
        for col in range(12):
            well_id = f'{chr(ord("A") + row)}{col+1:02d}'
            
            if well_id not in well_info:
                continue
            
            info = well_info[well_id]
            df = info['df']
            total_droplets = len(df)
            
            pos_ch1 = row * 24 + col * 2 + 1
            pos_ch2 = row * 24 + col * 2 + 2
            
            ax1 = plt.subplot(8, 24, pos_ch1)
            ax2 = plt.subplot(8, 24, pos_ch2)
            
            # --- Colors & Plotting (Same as before) ---
            n_bands_ch1 = len(info.get('ch1_thresholds', [])) + 1
            n_bands_ch2 = len(info.get('ch2_thresholds', [])) + 1
            grad_ch1 = create_gradient_colors('#0000FF', n_bands_ch1)
            grad_ch2 = create_gradient_colors('#00FF00', n_bands_ch2)
            
            np.random.seed(42)
            
            # Plot Ch1
            for band in info['ch1_bands']:
                mask = (df['Ch1_Amplitude'] >= band['min']) & (df['Ch1_Amplitude'] <= band['max'])
                band_data = df[mask]
                if len(band_data) > 0:
                    x_positions = np.random.uniform(-0.4, 0.4, len(band_data))
                    c_idx = min(band.get('band_index', 0), len(grad_ch1)-1)
                    ax1.scatter(x_positions, band_data['Ch1_Amplitude'],
                              c=grad_ch1[c_idx], alpha=0.6, s=1.5, edgecolors='none')
            
            # Lines Ch1
            if 'ch1_thresholds' in info:
                for threshold in info['ch1_thresholds']:
                    ax1.axhline(y=threshold, color='red', linestyle='--', linewidth=0.5, alpha=0.7)

            # Plot Ch2
            for band in info['ch2_bands']:
                mask = (df['Ch2_Amplitude'] >= band['min']) & (df['Ch2_Amplitude'] <= band['max'])
                band_data = df[mask]
                if len(band_data) > 0:
                    x_positions = np.random.uniform(-0.4, 0.4, len(band_data))
                    c_idx = min(band.get('band_index', 0), len(grad_ch2)-1)
                    ax2.scatter(x_positions, band_data['Ch2_Amplitude'],
                              c=grad_ch2[c_idx], alpha=0.6, s=1.5, edgecolors='none')

            # Lines Ch2
            if 'ch2_thresholds' in info:
                for threshold in info['ch2_thresholds']:
                    ax2.axhline(y=threshold, color='red', linestyle='--', linewidth=0.5, alpha=0.7)
            
            # --- Limits ---
            ax1.set_xlim(-0.5, 0.5)
            ax1.set_xticks([])
            ax1.set_yticks([])
            for spine in ax1.spines.values(): spine.set_visible(False)
            if info.get('limits'): ax1.set_ylim(info['limits']['ch1_limits'])
            else: ax1.set_ylim(ch1_min_global, ch1_max_global)

            ax2.set_xlim(-0.5, 0.5)
            ax2.set_xticks([])
            ax2.set_yticks([])
            for spine in ax2.spines.values(): spine.set_visible(False)
            if info.get('limits'): ax2.set_ylim(info['limits']['ch2_limits'])
            else: ax2.set_ylim(ch2_min_global, ch2_max_global)
            
            # --- Text Labels (Updated) ---
            ax1.text(0, 1.05, well_id, transform=ax1.transAxes, 
                    fontsize=6, fontweight='bold', ha='center', va='bottom')
            
            # Create stacked string for stats: "100 (10%)"
            # We use a very small font and \n to stack them if needed, or | separator
            # Given narrow width, we'll try a condensed vertical stack or tight list
            
            # Helper to make concise string: "1.2k (50%)"
            def concise_stats(band):
                return format_population_stats(band['count'], total_droplets)

            ch1_stats_list = [concise_stats(b) for b in info['ch1_bands']]
            # Reverse order so top band is on top textually
            ch1_text = '\n'.join(reversed(ch1_stats_list))
            
            ax1.text(0.5, 0.02, ch1_text, transform=ax1.transAxes,
                    fontsize=3.5, ha='center', va='bottom', color='#000080', 
                    family='monospace', fontweight='bold',
                    bbox=dict(facecolor='white', alpha=0.6, linewidth=0, pad=0.5))
            
            ch2_stats_list = [concise_stats(b) for b in info['ch2_bands']]
            ch2_text = '\n'.join(reversed(ch2_stats_list))
            
            ax2.text(0.5, 0.02, ch2_text, transform=ax2.transAxes,
                    fontsize=3.5, ha='center', va='bottom', color='#008000', 
                    family='monospace', fontweight='bold',
                    bbox=dict(facecolor='white', alpha=0.6, linewidth=0, pad=0.5))

    # Clean up grid labels
    for col in range(12):
        ax = plt.subplot(8, 24, col * 2 + 1)
        ax.text(1, 1.15, f'{col+1:02d}', transform=ax.transAxes,
               fontsize=7, ha='center', va='bottom', fontweight='bold')
    
    for row in range(8):
        ax = plt.subplot(8, 24, row * 24 + 1)
        row_letter = chr(ord('A') + row)
        ax.text(-0.3, 0.5, row_letter, transform=ax.transAxes,
               fontsize=7, ha='right', va='center', fontweight='bold')
               
    for row in range(8):
        for col in range(12):
            well_id = f'{chr(ord("A") + row)}{col+1:02d}'
            if well_id not in well_info:
                plt.subplot(8, 24, row * 24 + col * 2 + 1).axis('off')
                plt.subplot(8, 24, row * 24 + col * 2 + 2).axis('off')

    plt.suptitle('96-Well Plate - 1D Band Detection (Stats: Count + %)', 
                fontsize=13, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0.02, 0, 1, 0.99])
    
    if isinstance(destination, (str, Path)):
        fig.savefig(destination, dpi=300, bbox_inches='tight')
    else:
        destination.savefig(fig, dpi=200, bbox_inches='tight')
    plt.close(fig)


def create_96well_2d_page(well_info, destination, log_func):
    """
    Updated: Uses a floating legend list to show stats for ALL populations (no corners).
    """
    all_ch1 = [v for info in well_info.values() for v in info['df']['Ch1_Amplitude']]
    all_ch2 = [v for info in well_info.values() for v in info['df']['Ch2_Amplitude']]
    ch1_min_global, ch1_max_global = min(all_ch1), max(all_ch1)
    ch2_min_global, ch2_max_global = min(all_ch2), max(all_ch2)
    
    fig = plt.figure(figsize=(20, 12))
    
    for row in range(8):
        for col in range(12):
            well_id = f'{chr(ord("A") + row)}{col+1:02d}'
            
            if well_id not in well_info:
                continue
            
            info = well_info[well_id]
            df = info['df']
            total_droplets = len(df)
            
            pos = row * 12 + col + 1
            ax = plt.subplot(8, 12, pos)
            
            n_bands_ch1 = len(info.get('ch1_thresholds', [])) + 1
            grad_ch1 = create_gradient_colors('#0000FF', n_bands_ch1)
            
            ch1_bands = info['ch1_bands']
            ch2_bands = info['ch2_bands']
            
            # Store stats for the legend
            stats_entries = []
            
            for ch1_idx, ch1_band in enumerate(ch1_bands):
                ch1_mask = (df['Ch1_Amplitude'] >= ch1_band['min']) & (df['Ch1_Amplitude'] <= ch1_band['max'])
                
                for ch2_idx, ch2_band in enumerate(ch2_bands):
                    ch2_mask = (df['Ch2_Amplitude'] >= ch2_band['min']) & (df['Ch2_Amplitude'] <= ch2_band['max'])
                    combined_mask = ch1_mask & ch2_mask
                    band_data = df[combined_mask]
                    count = len(band_data)
                    
                    if count > 0:
                        c_idx = min(ch1_band.get('band_index', 0), len(grad_ch1)-1)
                        color = grad_ch1[c_idx]
                        alpha = 0.3 + (0.4 * ch2_idx / max(len(ch2_bands) - 1, 1))
                        
                        ax.scatter(band_data['Ch1_Amplitude'], 
                                  band_data['Ch2_Amplitude'],
                                  c=color, alpha=alpha, s=0.5, edgecolors='none')
                        
                        # Name logic for stats list
                        if ch1_idx == 0 and ch2_idx == 0: name = "Neg"
                        elif ch1_idx > 0 and ch2_idx == 0: name = "Fam+"
                        elif ch1_idx == 0 and ch2_idx > 0: name = "Hex+"
                        elif ch1_idx > 0 and ch2_idx > 0: name = "Dbl+"
                        else: name = f"B{ch1_idx}/{ch2_idx}" # Fallback for complex multiplex
                        
                        # Add to list
                        stats_str = format_population_stats(count, total_droplets)
                        stats_entries.append(f"{name}: {stats_str}")

            # Draw Lines
            if 'ch1_thresholds' in info:
                for threshold in info['ch1_thresholds']:
                    ax.axvline(x=threshold, color='red', linestyle='--', linewidth=0.3, alpha=0.5)
            if 'ch2_thresholds' in info:
                for threshold in info['ch2_thresholds']:
                    ax.axhline(y=threshold, color='red', linestyle='--', linewidth=0.3, alpha=0.5)
            
            # Limits
            ax.set_xticks([])
            ax.set_yticks([])
            if info.get('limits'):
                ax.set_xlim(info['limits']['ch1_limits'])
                ax.set_ylim(info['limits']['ch2_limits'])
            else:
                ax.set_xlim(ch1_min_global, ch1_max_global)
                ax.set_ylim(ch2_min_global, ch2_max_global)
            
            ax.text(0.5, 1.02, well_id, transform=ax.transAxes, 
                    fontsize=6, fontweight='bold', ha='center', va='bottom')
            
            # --- CREATE LEGEND BLOCK ---
            # Sort: Neg first, then others
            stats_entries.sort(key=lambda x: 0 if "Neg" in x else 1)
            
            full_legend_text = '\n'.join(stats_entries)
            
            # Place in top-right corner with background
            ax.text(0.96, 0.96, full_legend_text, transform=ax.transAxes,
                   fontsize=3.5, ha='right', va='top', family='monospace',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.75, linewidth=0, pad=0.3))

    # Headers and cleanup
    for col in range(12):
        ax = plt.subplot(8, 12, col + 1)
        ax.text(0.5, 1.12, f'{col+1:02d}', transform=ax.transAxes,
               fontsize=7, ha='center', va='bottom', fontweight='bold')
    for row in range(8):
        ax = plt.subplot(8, 12, row * 12 + 1)
        row_letter = chr(ord('A') + row)
        ax.text(-0.15, 0.5, row_letter, transform=ax.transAxes,
               fontsize=7, ha='right', va='center', fontweight='bold')
    for row in range(8):
        for col in range(12):
            well_id = f'{chr(ord("A") + row)}{col+1:02d}'
            if well_id not in well_info:
                plt.subplot(8, 12, row * 12 + col + 1).axis('off')

    plt.suptitle('96-Well Plate - 2D Scatter (Full Stats in Legend)', 
                fontsize=13, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0.02, 0, 1, 0.99])
    
    if isinstance(destination, (str, Path)):
        fig.savefig(destination, dpi=300, bbox_inches='tight')
    else:
        destination.savefig(fig, dpi=200, bbox_inches='tight')
    plt.close(fig)


def calculate_well_band_statistics(well_id, df, ch1_bands=None, ch2_bands=None, use_gmm=True):
    """
    Calculate band statistics for a well.
    ALWAYS uses 2D GMM first if bands not provided (2D-first approach).
    """
    if ch1_bands is None or ch2_bands is None:
        # Use 2D GMM to detect populations
        result = detect_populations_2d_gmm_robust(
            df['Ch1_Amplitude'].values,
            df['Ch2_Amplitude'].values,
            max_populations=4
        )
        
        # Apply thresholds to get bands
        ch1_bands = apply_thresholds_to_get_bands(
            df['Ch1_Amplitude'].values,
            result['ch1_thresholds'],
            percentile_bounds=True
        )
        
        ch2_bands = apply_thresholds_to_get_bands(
            df['Ch2_Amplitude'].values,
            result['ch2_thresholds'],
            percentile_bounds=True
        )
    
    # Calculate statistics for each band combination
    stats_rows = []
    
    # Get channel names from DataFrame attributes if available
    channel_map = df.attrs.get('channel_map', {'Ch1_Amplitude': 'Ch1', 'Ch2_Amplitude': 'Ch2'})
    ch1_name = channel_map.get('Ch1_Amplitude', 'Ch1')
    ch2_name = channel_map.get('Ch2_Amplitude', 'Ch2')
    
    for ch1_idx, ch1_band in enumerate(ch1_bands):
        ch1_mask = (df['Ch1_Amplitude'] >= ch1_band['min']) & (df['Ch1_Amplitude'] <= ch1_band['max'])
        
        for ch2_idx, ch2_band in enumerate(ch2_bands):
            ch2_mask = (df['Ch2_Amplitude'] >= ch2_band['min']) & (df['Ch2_Amplitude'] <= ch2_band['max'])
            
            combined_mask = ch1_mask & ch2_mask
            count = combined_mask.sum()
            
            if count > 0:
                stats_rows.append({
                    'Well': well_id,
                    f'{ch1_name}_Band': ch1_idx + 1,
                    f'{ch2_name}_Band': ch2_idx + 1,
                    'Count': int(count),
                    'Percentage': float(count / len(df) * 100),
                    f'{ch1_name}_Center': float(ch1_band['center']),
                    f'{ch2_name}_Center': float(ch2_band['center']),
                })
    
    return stats_rows, ch1_bands, ch2_bands


def process_ddpcr_files(csv_files, log_func=print, output_format='png'):
    """Process multiple ddPCR CSV files and generate plots and statistics"""
    
    if not csv_files:
        raise ValueError("No CSV files provided")
    
    first_file = Path(csv_files[0])
    base_dir = first_file.parent
    output_dir = base_dir / "plots"
    output_dir.mkdir(exist_ok=True)
    
    log_func(f"\nOutput directory: {output_dir}")
    
    all_stats_rows = []
    all_data = {}
    
    for csv_file in csv_files:
        csv_path = Path(csv_file)
        log_func(f"\nProcessing: {csv_path.name}")
        
        try:
            df = pd.read_csv(csv_path)
            log_func(f"  Loaded {len(df)} droplets")
            
            # Detect and normalize column names
            # Handle both old format (Ch1_Amplitude) and new format (FAM_Amplitude, HEX_Amplitude, etc.)
            if 'Ch1 Amplitude' in df.columns:
                df = df.rename(columns={'Ch1 Amplitude': 'Ch1_Amplitude', 'Ch2 Amplitude': 'Ch2_Amplitude'})
            
            # Extract actual channel names from column headers
            amplitude_cols = [col for col in df.columns if col.endswith('_Amplitude')]
            
            if len(amplitude_cols) >= 2:
                # New format with actual channel names
                ch1_col = amplitude_cols[0]
                ch2_col = amplitude_cols[1]
                ch1_name = ch1_col.replace('_Amplitude', '')
                ch2_name = ch2_col.replace('_Amplitude', '')
                
                # Rename to standard Ch1/Ch2 for internal processing
                df = df.rename(columns={ch1_col: 'Ch1_Amplitude', ch2_col: 'Ch2_Amplitude'})
                
                # Store channel names in DataFrame attributes
                df.attrs['channel_names'] = [ch1_name, ch2_name]
                df.attrs['channel_map'] = {
                    'Ch1_Amplitude': ch1_name,
                    'Ch2_Amplitude': ch2_name
                }
            else:
                # Old format or default
                df.attrs['channel_names'] = ['Ch1', 'Ch2']
                df.attrs['channel_map'] = {
                    'Ch1_Amplitude': 'Ch1',
                    'Ch2_Amplitude': 'Ch2'
                }
            
            # Extract well ID from filename
            well_match = re.search(r'([A-H])(\d{1,2})', csv_path.stem)
            if well_match:
                row_letter = well_match.group(1)
                col_num = int(well_match.group(2))
                well_id = f'{row_letter}{col_num:02d}'
            else:
                well_id = csv_path.stem
            
            all_data[well_id] = df
            
            # Calculate band statistics for this well
            stats_rows, ch1_bands, ch2_bands = calculate_well_band_statistics(well_id, df)
            all_stats_rows.extend(stats_rows)
            
            # Log summary
            log_func(f"  Ch1: {len(ch1_bands)} bands, Ch2: {len(ch2_bands)} bands")
            
            # Create individual well plot with band coloring
            fig = plot_well_with_bands(df, well_id, ch1_bands, ch2_bands)
            
            # DYNAMIC EXTENSION
            output_plot = output_dir / f"{well_id}_1D_plot.{output_format}"
            fig.savefig(output_plot, format=output_format, dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            # Create 2D scatter plot
            fig_2d = plot_well_2d_scatter(df, well_id, ch1_bands, ch2_bands)
            output_2d = output_dir / f"{well_id}_2D_plot.{output_format}"
            fig_2d.savefig(output_2d, format=output_format, dpi=300, bbox_inches='tight')
            plt.close(fig_2d)
            
            log_func(f"  Saved: {output_plot.name} and {output_2d.name}")
            
        except Exception as e:
            log_func(f"  ⚠️  Error processing {csv_path.name}: {str(e)}")
            import traceback
            log_func(traceback.format_exc())
            continue
    
    # Create 96-well plate overview
    if all_data:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        # Note: We pass the base name, the function will handle extension for PNG vs PDF
        plate_output = output_dir / f"96well_plate_overview_{timestamp}" 
        plot_96well_layout(all_data, plate_output, log_func, None, True, output_format)

    stats_file = None
    if all_stats_rows:
        stats_df = pd.DataFrame(all_stats_rows)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        stats_file = output_dir / f"band_statistics_{timestamp}.csv"
        stats_df.to_csv(stats_file, index=False)
        log_func(f"\nBand statistics saved to: {stats_file.name}")
    
    return output_dir, stats_file


def process_qlp_file(qlp_file, well_assignments, include_filtered, log_func=print, output_format='png'):
    """Process a QLP file with well assignments to conditions"""
    
    qlp_path = Path(qlp_file)
    base_dir = qlp_path.parent
    output_dir = base_dir / "plots"
    output_dir.mkdir(exist_ok=True)
    
    # Create subdirectory for individual well CSVs
    csv_dir = output_dir / "well_csvs"
    csv_dir.mkdir(exist_ok=True)
    
    log_func(f"\nOutput directory: {output_dir}")
    log_func(f"Well CSV directory: {csv_dir}")
    log_func(f"Parsing QLP file: {qlp_path.name}")
    log_func(f"Include filtered droplets (Cluster 0): {include_filtered}")
    
    # Parse QLP file
    parser = BioRadQLPParser(qlp_file)
    well_dataframes = parser.parse_to_dataframe(include_filtered=include_filtered)
    
    log_func(f"Extracted {len(well_dataframes)} wells")
    
    # Calculate shared thresholds if we have condition assignments
    use_shared_thresholds = well_assignments and len(set(well_assignments.values())) > 0
    shared_thresh = None
    
    if use_shared_thresholds:
        log_func("\nCalculating shared thresholds per condition...")
        shared_thresh = calculate_shared_thresholds(
            well_dataframes,
            well_assignments,
            method='2d_gmm'
        )
        log_func(f"  Calculated thresholds for {len(shared_thresh)} condition(s)")
        
        # Log thresholds for each condition
        for condition, thresholds in shared_thresh.items():
            ch1_str = ', '.join([f'{t:.1f}' for t in thresholds['ch1_thresholds']])
            ch2_str = ', '.join([f'{t:.1f}' for t in thresholds['ch2_thresholds']])
            log_func(f"    {condition}: Ch1=[{ch1_str}], Ch2=[{ch2_str}]")
    
    # Export individual well CSVs and create plots
    all_stats_rows = []
    
    for well_id, df in sorted(well_dataframes.items()):
        # Get channel names for this well
        channel_map = df.attrs.get('channel_map', {})
        ch1_name = channel_map.get('Ch1_Amplitude', 'Ch1')
        ch2_name = channel_map.get('Ch2_Amplitude', 'Ch2')
        
        # Rename columns to use actual channel names for export
        df_export = df.copy()
        df_export = df_export.rename(columns={
            'Ch1_Amplitude': f'{ch1_name}_Amplitude',
            'Ch2_Amplitude': f'{ch2_name}_Amplitude'
        })
        
        # Export CSV for this well
        csv_path = csv_dir / f"{well_id}.csv"
        df_export.to_csv(csv_path, index=False)
        
        # Count cluster distribution
        cluster_counts = df['Cluster'].value_counts().to_dict()
        cluster_str = ', '.join([f"C{c}:{cluster_counts.get(c, 0)}" for c in sorted(cluster_counts.keys())])
        log_func(f"  Exported: {well_id}.csv ({len(df)} droplets - {cluster_str}) [{ch1_name}, {ch2_name}]")
        
        # Calculate band statistics using shared thresholds if available
        if use_shared_thresholds:
            condition = well_assignments.get(well_id, 'Unassigned')
            thresholds = shared_thresh.get(condition, {'ch1_thresholds': [], 'ch2_thresholds': []})
            
            ch1_bands = apply_thresholds_to_get_bands(
                df['Ch1_Amplitude'].values,
                thresholds['ch1_thresholds'],
                percentile_bounds=True
            )
            
            ch2_bands = apply_thresholds_to_get_bands(
                df['Ch2_Amplitude'].values,
                thresholds['ch2_thresholds'],
                percentile_bounds=True
            )
            
            stats_rows, _, _ = calculate_well_band_statistics(well_id, df, ch1_bands, ch2_bands)
        else:
            # Individual well analysis
            stats_rows, ch1_bands, ch2_bands = calculate_well_band_statistics(well_id, df, use_gmm=True)
        
        all_stats_rows.extend(stats_rows)
        
        # Add condition info to stats if assigned
        condition = well_assignments.get(well_id, 'Unassigned')
        for row in stats_rows:
            row['Condition'] = condition
        
        # Create individual well 1D plot with band coloring
        fig = plot_well_with_bands(df, well_id, ch1_bands, ch2_bands)
        
        # Add condition to title if assigned
        if condition != 'Unassigned':
            fig.suptitle(f'Well {well_id} ({condition}) - Band Detection', fontsize=16, fontweight='bold')
        
        plot_path = output_dir / f"{well_id}_1D_plot.{output_format}"
        fig.savefig(plot_path, format=output_format, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        # Create individual well 2D scatter plot
        fig_2d = plot_well_2d_scatter(df, well_id, ch1_bands, ch2_bands)
        
        # Add condition to title if assigned
        if condition != 'Unassigned':
            fig_2d.axes[0].set_title(f'Well {well_id} ({condition}) - 2D Scatter (Ch1 vs Ch2)', 
                                     fontsize=14, fontweight='bold')
        
        plot_2d_path = output_dir / f"{well_id}_2D_plot.{output_format}"
        fig_2d.savefig(plot_2d_path, format=output_format, dpi=300, bbox_inches='tight')
        plt.close(fig_2d)
    
    log_func(f"\nGenerated {len(well_dataframes)} individual well plots (1D and 2D)")
    
    # Create 96-well plate overview
    log_func(f"\nCreating 96-well plate overview...")
    plate_data = {}
    for well_id, df in well_dataframes.items():
        plate_data[well_id] = df
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # Pass base path without extension
    plate_output = output_dir / f"{qlp_path.stem}_96well_plate_{timestamp}"
    plot_96well_layout(plate_data, plate_output, log_func,
                       condition_assignments=well_assignments,
                       use_shared_thresholds=use_shared_thresholds,
                       output_format=output_format)
    
    # Save band statistics
    stats_file = None
    if all_stats_rows:
        stats_df = pd.DataFrame(all_stats_rows)
        filtered_suffix = "_all" if include_filtered else "_accepted"
        stats_file = output_dir / f"{qlp_path.stem}_band_statistics{filtered_suffix}_{timestamp}.csv"
        stats_df.to_csv(stats_file, index=False)
        log_func(f"\nBand statistics saved to: {stats_file.name}")
        log_func(f"Individual well CSVs saved to: {csv_dir.name}/")
    
    return output_dir, stats_file


def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()