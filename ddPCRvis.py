#!/usr/bin/env python3
"""
ddPCRvis - ddPCR Visualization Tool
A GUI application for visualizing ddPCR amplitude data from Bio-Rad CSV exports and QLP binary files
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
    """Comprehensive parser for BioRad .QLP binary files - extracts ALL droplet data"""
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
        self.TAG_WELL_NAME = 65019  # Contains only row letter (A-H)
        self.TAG_WELL_ID = 65018     # Contains complete well ID (A01-H12)
        self.TAG_CHANNEL_NAMES = 65004  # Contains channel names (Ch1,Ch2)
        self.TAG_DATA_START = 65021
        self.TAG_CLUSTER_ARRAY = 65057
        self.TAG_QUALITY_ARRAY = 65054  # Quality flags per droplet
        self.TAG_WELL_QUALITY = 65005   # Per-well quality (0.0 or 600.0)
        self.TAG_QUALITY_SCORE = 65065  # Per-well quality score (0.0-1.0)
        self.TAG_QUALITY_THRESHOLD = 65078  # Quality threshold (0.85)
        self.TAG_SOFTWARE = 305
        
        # Channel names (parsed from file)
        self.channel_names = None
        
        # Cluster mapping
        self.CLUSTER_BYTE_MAP = {
            0x00: 0,  # Filtered/Invalid
            0x11: 1,  # Cluster 1
            0x22: 2,  # Cluster 2
            0x33: 3,  # Cluster 3
            0x44: 4,  # Cluster 4
        }
        
        self.metadata = {}
        self.well_metadata = {}  # Per-well metadata
        
        # Find all valid well IDs at initialization
        self._find_valid_well_ids()
    
    def _find_valid_well_ids(self):
        """Find all valid well IDs in the file using the 02 01 pattern"""
        pattern = re.compile(rb'\x02\x01([A-H])(0[1-9]|1[0-2])\x00')
        
        self.valid_well_ids = []
        for match in pattern.finditer(self.data):
            row = match.group(1).decode('ascii')
            col = match.group(2).decode('ascii')
            well_id = f"{row}{col}"
            self.valid_well_ids.append(well_id)
        
        if self.debug:
            print(f"\nFound {len(self.valid_well_ids)} valid well IDs: {self.valid_well_ids}")

    def _get_val(self, offset, fmt):
        try:
            return struct.unpack(f"{self.endian}{fmt}", self.data[offset:offset+struct.calcsize(fmt)])[0]
        except: 
            return None
    
    def _extract_tags(self, ifd_offset):
        """Extract all tags from an IFD"""
        num_entries = self._get_val(ifd_offset, "H")
        tags = {}
        
        for i in range(num_entries):
            entry_ptr = ifd_offset + 2 + (i * 12)
            tid = self._get_val(entry_ptr, "H")
            ttype = self._get_val(entry_ptr + 2, "H")
            count = self._get_val(entry_ptr + 4, "I")
            size = {1:1, 2:1, 3:2, 4:4, 5:8, 7:1, 9:4, 11:4, 12:8}.get(ttype, 1) * count
            ptr = entry_ptr + 8 if size <= 4 else self._get_val(entry_ptr + 8, "I")
            tags[tid] = (ptr, size, ttype, count)
        
        return tags
    
    def _extract_metadata(self):
        """Extract file metadata"""
        ifd_offset = self._get_val(4, "I")
        tags = self._extract_tags(ifd_offset)
        
        if self.TAG_SOFTWARE in tags:
            ptr, size, _, _ = tags[self.TAG_SOFTWARE]
            try:
                value = self.data[ptr:ptr+size].decode('ascii', errors='ignore').split('\x00')[0].strip()
                self.metadata['software'] = value
            except:
                pass
        
        # Try to extract channel names from Tag 65004 first
        if self.TAG_CHANNEL_NAMES in tags:
            ptr, size, _, _ = tags[self.TAG_CHANNEL_NAMES]
            try:
                channel_str = self.data[ptr:ptr+size].decode('ascii', errors='ignore').rstrip('\x00')
                self.channel_names = [ch.strip() for ch in channel_str.split(',')]
                self.metadata['channel_names'] = self.channel_names
                if self.debug:
                    print(f"Found channel names from Tag 65004: {self.channel_names}")
            except:
                pass
        
        # If Tag 65004 didn't work or gave generic names, try the "Unknown" pattern
        if self.channel_names is None or self.channel_names == ['Ch1', 'Ch2']:
            self._extract_channel_names_from_unknown_pattern()
        
        # Fallback to default channel names if not found
        if self.channel_names is None:
            self.channel_names = ['Ch1', 'Ch2']
            if self.debug:
                print(f"Using default channel names: {self.channel_names}")
    
    def _extract_channel_names_from_unknown_pattern(self):
        """
        Extract channel names using the Unknown + 32 bytes pattern.
        Structure: [padding] "Unknown" [padding] [channel_name]
        """
        channels_found = []
        idx = 0
        
        while idx < len(self.data):
            idx = self.data.find(b'Unknown\x00', idx)
            if idx == -1:
                break
            
            # Check if preceded by zeros (padding)
            if idx >= 32:
                preceding = self.data[idx-32:idx]
                if preceding.count(b'\x00') >= 28:  # Mostly zeros
                    # Channel name is 32 bytes after "Unknown\x00" (8 bytes)
                    channel_name_offset = idx + 32
                    channel_data = self.data[channel_name_offset:channel_name_offset+64]
                    
                    # Find null terminator
                    null_pos = channel_data.find(b'\x00')
                    if null_pos > 0:
                        channel_name = channel_data[:null_pos].decode('ascii', errors='ignore')
                        
                        # Validate it's a reasonable channel name
                        if len(channel_name) >= 2 and channel_name.isprintable() and channel_name != 'Unknown':
                            channels_found.append(channel_name)
                            if self.debug:
                                print(f"Found channel '{channel_name}' at offset {channel_name_offset}")
            
            idx += 1
        
        # Use first two unique channels
        unique_channels = []
        seen = set()
        for name in channels_found:
            if name not in seen:
                unique_channels.append(name)
                seen.add(name)
                if len(unique_channels) >= 2:
                    break
        
        if len(unique_channels) >= 1:
            if self.channel_names is None:
                self.channel_names = ['Ch1', 'Ch2']
            self.channel_names[0] = unique_channels[0]
            self.metadata['channel_1_name'] = unique_channels[0]
            
        if len(unique_channels) >= 2:
            self.channel_names[1] = unique_channels[1]
            self.metadata['channel_2_name'] = unique_channels[1]
        
        if self.debug and unique_channels:
            print(f"Extracted channel names from Unknown pattern: {self.channel_names}")

    def parse_to_dataframe(self, include_filtered=True):
        """
        Parse QLP file and return dict of well_id -> DataFrame
        
        Args:
            include_filtered: If True, include droplets with cluster=0 (default: True)
        """
        self._extract_metadata()
        
        ifd_offset = self._get_val(4, "I")
        well_data = defaultdict(list)
        well_index = 0  # Index into self.valid_well_ids list

        while ifd_offset != 0 and ifd_offset < len(self.data):
            tags = self._extract_tags(ifd_offset)
            num_entries = self._get_val(ifd_offset, "H")
            
            # Only process IFDs with well name, data, and cluster arrays
            if (self.TAG_WELL_NAME in tags and 
                self.TAG_DATA_START in tags and 
                self.TAG_CLUSTER_ARRAY in tags):
                
                # Assign the next valid well ID in sequential order
                if well_index < len(self.valid_well_ids):
                    well_id = self.valid_well_ids[well_index]
                    well_index += 1
                else:
                    # No more valid well IDs - skip this IFD
                    if self.debug:
                        print(f"  Skipping IFD at {ifd_offset} - no more well IDs")
                    ifd_offset = self._get_val(ifd_offset + 2 + (num_entries * 12), "I")
                    continue
                
                if self.debug:
                    print(f"Processing well: {well_id}")
                
                # Extract per-well metadata
                well_meta = {}
                
                if self.TAG_WELL_QUALITY in tags:
                    ptr, _, _, _ = tags[self.TAG_WELL_QUALITY]
                    well_meta['well_quality_flag'] = self._get_val(ptr, 'f')
                
                if self.TAG_QUALITY_SCORE in tags:
                    ptr, _, _, _ = tags[self.TAG_QUALITY_SCORE]
                    well_meta['well_quality_score'] = self._get_val(ptr, 'f')
                
                if self.TAG_QUALITY_THRESHOLD in tags:
                    ptr, _, _, _ = tags[self.TAG_QUALITY_THRESHOLD]
                    well_meta['quality_threshold'] = self._get_val(ptr, 'f')
                
                self.well_metadata[well_id] = well_meta
                
                # Read cluster array
                tag_info = tags[self.TAG_CLUSTER_ARRAY]
                cluster_array_size = tag_info[1]
                
                # Skip wells with empty cluster arrays
                if cluster_array_size == 0:
                    ifd_offset = self._get_val(ifd_offset + 2 + (num_entries * 12), "I")
                    continue
                
                cluster_array = self.data[tag_info[0]:tag_info[0]+cluster_array_size]
                
                # Read quality flag array (Tag 65054) if present
                quality_array = None
                if self.TAG_QUALITY_ARRAY in tags:
                    tag_info = tags[self.TAG_QUALITY_ARRAY]
                    if tag_info[1] > 0:
                        quality_array = self.data[tag_info[0]:tag_info[0]+tag_info[1]]
                
                # Read droplet data
                cursor = tags[self.TAG_DATA_START][0]
                droplet_idx = 0
                
                while cursor < len(self.data) - self.RECORD_SIZE and droplet_idx < len(cluster_array):
                    # Parse complete 28-byte record
                    r = struct.unpack(f"{self.endian}{self.RECORD_FMT}", 
                                     self.data[cursor:cursor+self.RECORD_SIZE])
                    
                    # Get cluster assignment
                    cluster_byte = cluster_array[droplet_idx]
                    cluster = self.CLUSTER_BYTE_MAP.get(cluster_byte, 0)
                    
                    # Get quality flag if available
                    quality_flag = None
                    if quality_array and droplet_idx < len(quality_array):
                        quality_flag = quality_array[droplet_idx]
                    
                    # Include ALL droplets or filter based on parameter
                    if include_filtered or cluster > 0:
                        droplet = {
                            'DropletID': r[0],
                            'Ch1_Amplitude': r[1],
                            'Ch2_Amplitude': r[4],
                            'Cluster': cluster,
                            'Cluster_Byte': cluster_byte,
                            'Quality_Flag': quality_flag if quality_flag is not None else 0,
                            'Field_2': r[2],
                            'Field_3': r[3],
                            'Field_5': r[5],
                            'Field_6': r[6],
                        }
                        
                        # Add well metadata to each droplet
                        for key, value in well_meta.items():
                            droplet[key] = value
                        
                        well_data[well_id].append(droplet)
                    
                    droplet_idx += 1
                    cursor += self.RECORD_SIZE

            ifd_offset = self._get_val(ifd_offset + 2 + (num_entries * 12), "I")

        # Convert to DataFrames
        well_dataframes = {}
        for well_id, data in well_data.items():
            if data:
                df = pd.DataFrame(data)
                # Store channel names as DataFrame attribute for reference
                df.attrs['channel_names'] = self.channel_names if self.channel_names else ['Ch1', 'Ch2']
                df.attrs['channel_map'] = {
                    'Ch1_Amplitude': self.channel_names[0] if self.channel_names and len(self.channel_names) > 0 else 'Ch1',
                    'Ch2_Amplitude': self.channel_names[1] if self.channel_names and len(self.channel_names) > 1 else 'Ch2'
                }
                well_dataframes[well_id] = df
        
        return well_dataframes

    def get_channel_names(self):
        """Return the channel names parsed from the file"""
        if self.channel_names is None:
            self._extract_metadata()
        return self.channel_names if self.channel_names else ['Ch1', 'Ch2']


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
    
    def __init__(self, files, file_type='csv', well_assignments=None, include_filtered=True):
        super().__init__()
        self.files = files
        self.file_type = file_type
        self.well_assignments = well_assignments
        self.include_filtered = include_filtered
        
    def run(self):
        try:
            if self.file_type == 'qlp':
                output_dir, stats_file = process_qlp_file(
                    self.files[0], 
                    self.well_assignments,
                    self.include_filtered,
                    self.progress.emit
                )
            else:
                output_dir, stats_file = process_ddpcr_files(
                    self.files, 
                    self.progress.emit
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
        self.progress_bar.setRange(0, 0)  # Indeterminate
        
        self.thread = ProcessingThread(files, file_type, well_assignments, include_filtered)
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


def plot_well_2d_scatter(df, well_id, ch1_bands, ch2_bands):
    """Create 2D scatter plot for Ch1 vs Ch2 with band coloring"""
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    
    # Get actual channel names from DataFrame attributes
    channel_map = df.attrs.get('channel_map', {'Ch1_Amplitude': 'Ch1', 'Ch2_Amplitude': 'Ch2'})
    ch1_name = channel_map.get('Ch1_Amplitude', 'Ch1')
    ch2_name = channel_map.get('Ch2_Amplitude', 'Ch2')
    
    # Generate gradient colors for bands
    ch1_colors = create_gradient_colors('#0000FF', len(ch1_bands))
    ch2_colors = create_gradient_colors('#00FF00', len(ch2_bands))
    
    # Create a color map based on which band each droplet belongs to
    # We'll color by Ch1 band primarily, with some transparency
    np.random.seed(42)
    
    # Plot each combination of Ch1 and Ch2 bands
    for ch1_idx, ch1_band in enumerate(ch1_bands):
        ch1_mask = (df['Ch1_Amplitude'] >= ch1_band['min']) & (df['Ch1_Amplitude'] <= ch1_band['max'])
        
        for ch2_idx, ch2_band in enumerate(ch2_bands):
            ch2_mask = (df['Ch2_Amplitude'] >= ch2_band['min']) & (df['Ch2_Amplitude'] <= ch2_band['max'])
            
            # Get droplets in both bands
            combined_mask = ch1_mask & ch2_mask
            band_data = df[combined_mask]
            
            if len(band_data) > 0:
                # Mix colors based on both channels
                # Use Ch1 color with varying intensity based on Ch2
                color = ch1_colors[ch1_idx]
                alpha = 0.3 + (0.4 * ch2_idx / max(len(ch2_bands) - 1, 1))
                
                ax.scatter(band_data['Ch1_Amplitude'], 
                          band_data['Ch2_Amplitude'],
                          c=color, alpha=alpha, s=5, edgecolors='none',
                          label=f'{ch1_name}-B{ch1_idx+1} × {ch2_name}-B{ch2_idx+1}' if ch1_idx < 2 and ch2_idx < 2 else '')
    
    # Draw separator lines for bands
    for ch1_idx in range(len(ch1_bands) - 1):
        separator_x = (ch1_bands[ch1_idx]['max'] + ch1_bands[ch1_idx + 1]['min']) / 2
        ax.axvline(x=separator_x, color='red', linestyle='--', linewidth=1, alpha=0.5)
    
    for ch2_idx in range(len(ch2_bands) - 1):
        separator_y = (ch2_bands[ch2_idx]['max'] + ch2_bands[ch2_idx + 1]['min']) / 2
        ax.axhline(y=separator_y, color='red', linestyle='--', linewidth=1, alpha=0.5)
    
    ax.set_xlabel(f'{ch1_name} Amplitude', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'{ch2_name} Amplitude', fontsize=12, fontweight='bold')
    ax.set_title(f'Well {well_id} - 2D Scatter ({ch1_name} vs {ch2_name})', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.2)
    
    # Add band count annotation
    ax.text(0.02, 0.98, f'{ch1_name}: {len(ch1_bands)} bands\n{ch2_name}: {len(ch2_bands)} bands', 
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    return fig


def plot_well_with_bands(df, well_id, ch1_bands, ch2_bands):
    """Create band-colored amplitude plot for a single well"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Get actual channel names from DataFrame attributes
    channel_map = df.attrs.get('channel_map', {'Ch1_Amplitude': 'Ch1', 'Ch2_Amplitude': 'Ch2'})
    ch1_name = channel_map.get('Ch1_Amplitude', 'Ch1')
    ch2_name = channel_map.get('Ch2_Amplitude', 'Ch2')
    
    # Generate gradient colors for bands
    ch1_colors = create_gradient_colors('#0000FF', len(ch1_bands))
    ch2_colors = create_gradient_colors('#00FF00', len(ch2_bands))
    
    np.random.seed(42)
    
    # Plot Ch1 with bands
    for band_idx, band in enumerate(ch1_bands):
        mask = (df['Ch1_Amplitude'] >= band['min']) & (df['Ch1_Amplitude'] <= band['max'])
        band_data = df[mask]
        
        if len(band_data) > 0:
            x_jitter = np.random.uniform(-0.3, 0.3, len(band_data))
            ax1.scatter(x_jitter, band_data['Ch1_Amplitude'], 
                       c=ch1_colors[band_idx], alpha=0.6, s=10,
                       label=f'Band {band_idx+1} ({band["proportion"]:.1f}%)')
            
            # Draw separator line between bands
            if band_idx < len(ch1_bands) - 1:
                next_band = ch1_bands[band_idx + 1]
                separator_y = (band['max'] + next_band['min']) / 2
                ax1.axhline(y=separator_y, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    
    ax1.set_ylabel(f'{ch1_name} Amplitude', fontsize=12)
    ax1.set_xlim(-0.5, 0.5)
    ax1.set_xticks([])
    ax1.set_title(f'Channel 1 ({ch1_name})', fontsize=14, fontweight='bold')
    if ch1_bands:
        ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot Ch2 with bands
    for band_idx, band in enumerate(ch2_bands):
        mask = (df['Ch2_Amplitude'] >= band['min']) & (df['Ch2_Amplitude'] <= band['max'])
        band_data = df[mask]
        
        if len(band_data) > 0:
            x_jitter = np.random.uniform(-0.3, 0.3, len(band_data))
            ax2.scatter(x_jitter, band_data['Ch2_Amplitude'],
                       c=ch2_colors[band_idx], alpha=0.6, s=10,
                       label=f'Band {band_idx+1} ({band["proportion"]:.1f}%)')
            
            # Draw separator line
            if band_idx < len(ch2_bands) - 1:
                next_band = ch2_bands[band_idx + 1]
                separator_y = (band['max'] + next_band['min']) / 2
                ax2.axhline(y=separator_y, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    
    ax2.set_ylabel(f'{ch2_name} Amplitude', fontsize=12)
    ax2.set_xlim(-0.5, 0.5)
    ax2.set_xticks([])
    ax2.set_title(f'Channel 2 ({ch2_name})', fontsize=14, fontweight='bold')
    if ch2_bands:
        ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle(f'Well {well_id} - Band Detection', fontsize=16, fontweight='bold')
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


def detect_bands_kmeans(values, max_bands=4, max_sample_size=10000):
    """Detect bands using k-means clustering (optimized for large datasets)"""
    from sklearn.cluster import KMeans
    
    if len(values) < 10:
        return []
    
    # For large datasets, sample for clustering but apply to all data
    use_sampling = len(values) > max_sample_size
    if use_sampling:
        sample_indices = np.random.choice(len(values), max_sample_size, replace=False)
        values_for_clustering = values.iloc[sample_indices]
    else:
        values_for_clustering = values
    
    best_bands = []
    best_score = -np.inf
    
    values_array = values_for_clustering.values.reshape(-1, 1)
    
    for n_clusters in range(1, min(max_bands + 1, len(values_for_clustering) + 1)):
        try:
            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            labels = kmeans.fit_predict(values_array)
            
            if n_clusters > 1:
                from sklearn.metrics import silhouette_score
                score = silhouette_score(values_array, labels)
            else:
                score = 0
            
            if score > best_score:
                best_score = score
                centers = sorted(kmeans.cluster_centers_.flatten())
                
                # Now apply the learned centers to ALL data
                all_values_array = values.values.reshape(-1, 1)
                
                # Predict cluster for all points
                if use_sampling:
                    # Re-fit on all data with the best n_clusters
                    kmeans_full = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
                    all_labels = kmeans_full.fit_predict(all_values_array)
                    centers = sorted(kmeans_full.cluster_centers_.flatten())
                else:
                    all_labels = labels
                
                # Get the original centers before sorting
                if use_sampling:
                    original_centers = kmeans_full.cluster_centers_.flatten()
                else:
                    original_centers = kmeans.cluster_centers_.flatten()
                
                bands = []
                for i, center in enumerate(centers):
                    # Find which ORIGINAL cluster index corresponds to this sorted center
                    # This is the actual label that kmeans assigned
                    original_cluster_idx = np.argmin(np.abs(original_centers - center))
                    
                    # Get all points assigned to this cluster using the ORIGINAL label
                    if use_sampling:
                        mask = all_labels == original_cluster_idx
                    else:
                        mask = labels == original_cluster_idx
                    
                    band_values = values[mask] if use_sampling else values.iloc[mask]
                    
                    if len(band_values) > 0:
                        bands.append({
                            'min': float(band_values.min()),
                            'max': float(band_values.max()),
                            'center': float(center),
                            'count': len(band_values),
                            'proportion': len(band_values) / len(values) * 100
                        })
                
                best_bands = bands
                
        except Exception as e:
            continue
    
    return best_bands


def plot_96well_layout(data_dict, output_path, log_func=print):
    """Create 96-well plate layout visualization with both 1D and 2D views"""
    try:
        from sklearn.cluster import KMeans
    except ImportError:
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-learn", "--break-system-packages"])
        from sklearn.cluster import KMeans
    
    log_func("\nCreating 96-well plate visualizations...")
    
    well_info = {}
    for filename, df in data_dict.items():
        well_match = re.search(r'([A-H])(\d{1,2})', filename)
        if well_match:
            row_letter = well_match.group(1)
            col_num = int(well_match.group(2))
            well_id = f'{row_letter}{col_num:02d}'
            
            ch1_bands = detect_bands_kmeans(df['Ch1_Amplitude'])
            ch2_bands = detect_bands_kmeans(df['Ch2_Amplitude'])
            
            well_info[well_id] = {
                'df': df,
                'ch1_bands': ch1_bands,
                'ch2_bands': ch2_bands
            }
    
    if not well_info:
        log_func("  ⚠️  No valid well identifiers found in filenames")
        return
    
    # Create PDF with multiple pages
    with PdfPages(output_path) as pdf:
        # Page 1: 1D band view (existing)
        create_96well_1d_page(well_info, pdf, log_func)
        
        # Page 2: 2D scatter view
        create_96well_2d_page(well_info, pdf, log_func)
    
    log_func(f"  ✓ 96-well plate plots saved (1D and 2D views)")


def create_96well_1d_page(well_info, pdf, log_func):
    """Create the 1D band detection page for 96-well plate"""
    ch1_band_colors = {}
    ch2_band_colors = {}
    
    for well_id, info in well_info.items():
        ch1_colors = create_gradient_colors('#0000FF', len(info['ch1_bands']))
        ch2_colors = create_gradient_colors('#00FF00', len(info['ch2_bands']))
        
        for idx, color in enumerate(ch1_colors):
            ch1_band_colors[(well_id, idx)] = color
        for idx, color in enumerate(ch2_colors):
            ch2_band_colors[(well_id, idx)] = color
    
    all_ch1 = [v for info in well_info.values() for v in info['df']['Ch1_Amplitude']]
    all_ch2 = [v for info in well_info.values() for v in info['df']['Ch2_Amplitude']]
    ch1_min, ch1_max = min(all_ch1), max(all_ch1)
    ch2_min, ch2_max = min(all_ch2), max(all_ch2)
    
    fig = plt.figure(figsize=(20, 12))
    
    for row in range(8):
        for col in range(12):
            well_id = f'{chr(ord("A") + row)}{col+1:02d}'
            
            if well_id not in well_info:
                continue
            
            info = well_info[well_id]
            
            pos_ch1 = row * 24 + col * 2 + 1
            pos_ch2 = row * 24 + col * 2 + 2
            
            ax1 = plt.subplot(8, 24, pos_ch1)
            ax2 = plt.subplot(8, 24, pos_ch2)
            
            np.random.seed(42)
            df = info['df']
            
            for band_idx, band in enumerate(info['ch1_bands']):
                mask = (df['Ch1_Amplitude'] >= band['min']) & (df['Ch1_Amplitude'] <= band['max'])
                band_data = df[mask]
                
                if len(band_data) > 0:
                    x_positions = np.random.uniform(-0.4, 0.4, len(band_data))
                    color = ch1_band_colors.get((well_id, band_idx), '#0000FF')
                    
                    ax1.scatter(x_positions, band_data['Ch1_Amplitude'],
                              c=color, alpha=0.6, s=1.5, edgecolors='none')
                    
                    if band_idx < len(info['ch1_bands']) - 1:
                        next_band = info['ch1_bands'][band_idx + 1]
                        separator_y = (band['max'] + next_band['min']) / 2
                        ax1.axhline(y=separator_y, color='red', linestyle='--', linewidth=0.5, alpha=0.7)
            
            for band_idx, band in enumerate(info['ch2_bands']):
                mask = (df['Ch2_Amplitude'] >= band['min']) & (df['Ch2_Amplitude'] <= band['max'])
                band_data = df[mask]
                
                if len(band_data) > 0:
                    x_positions = np.random.uniform(-0.4, 0.4, len(band_data))
                    color = ch2_band_colors.get((well_id, band_idx), '#00FF00')
                    
                    ax2.scatter(x_positions, band_data['Ch2_Amplitude'],
                              c=color, alpha=0.6, s=1.5, edgecolors='none')
                    
                    if band_idx < len(info['ch2_bands']) - 1:
                        next_band = info['ch2_bands'][band_idx + 1]
                        separator_y = (band['max'] + next_band['min']) / 2
                        ax2.axhline(y=separator_y, color='red', linestyle='--', linewidth=0.5, alpha=0.7)
            
            # Set axis properties for Ch1
            ax1.set_xlim(-0.5, 0.5)
            ax1.set_ylim(ch1_min, ch1_max)
            ax1.set_xticks([])
            ax1.set_yticks([])
            for spine in ax1.spines.values():
                spine.set_visible(False)

            # Set axis properties for Ch2
            ax2.set_xlim(-0.5, 0.5)
            ax2.set_ylim(ch2_min, ch2_max)
            ax2.set_xticks([])
            ax2.set_yticks([])
            for spine in ax2.spines.values():
                spine.set_visible(False)
            
            ax1.text(0, 1.05, well_id, transform=ax1.transAxes, 
                    fontsize=6, fontweight='bold', ha='center', va='bottom')
            
            ch1_stats = ' | '.join([f'{b["proportion"]:.0f}%' for b in info['ch1_bands']])
            ax1.text(0, -0.15, ch1_stats if ch1_stats else '0%', transform=ax1.transAxes,
                    fontsize=4.5, ha='center', va='top', color='#000080', family='monospace')
            
            ch2_stats = ' | '.join([f'{b["proportion"]:.0f}%' for b in info['ch2_bands']])
            ax2.text(0, -0.15, ch2_stats if ch2_stats else '0%', transform=ax2.transAxes,
                    fontsize=4.5, ha='center', va='top', color='#008000', family='monospace')
    
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
                pos_ch1 = row * 24 + col * 2 + 1
                pos_ch2 = row * 24 + col * 2 + 2
                ax1 = plt.subplot(8, 24, pos_ch1)
                ax2 = plt.subplot(8, 24, pos_ch2)
                ax1.axis('off')
                ax2.axis('off')

    plt.suptitle('96-Well Plate - 1D Band Detection\nLeft: FAM (Blue) | Right: HEX (Green) | Numbers: Band proportions | Red lines: Separators', 
                fontsize=13, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0.02, 0, 1, 0.99])
    pdf.savefig(fig, dpi=200, bbox_inches='tight')
    plt.close(fig)


def create_96well_2d_page(well_info, pdf, log_func):
    """Create the 2D scatter plot page for 96-well plate"""
    # Get global min/max for consistent axes
    all_ch1 = [v for info in well_info.values() for v in info['df']['Ch1_Amplitude']]
    all_ch2 = [v for info in well_info.values() for v in info['df']['Ch2_Amplitude']]
    ch1_min, ch1_max = min(all_ch1), max(all_ch1)
    ch2_min, ch2_max = min(all_ch2), max(all_ch2)
    
    fig = plt.figure(figsize=(20, 12))
    
    for row in range(8):
        for col in range(12):
            well_id = f'{chr(ord("A") + row)}{col+1:02d}'
            
            if well_id not in well_info:
                continue
            
            info = well_info[well_id]
            df = info['df']
            ch1_bands = info['ch1_bands']
            ch2_bands = info['ch2_bands']
            
            # Position in grid (one subplot per well)
            pos = row * 12 + col + 1
            ax = plt.subplot(8, 12, pos)
            
            # Generate colors
            ch1_colors = create_gradient_colors('#0000FF', len(ch1_bands))
            
            # Plot each combination of bands
            for ch1_idx, ch1_band in enumerate(ch1_bands):
                ch1_mask = (df['Ch1_Amplitude'] >= ch1_band['min']) & (df['Ch1_Amplitude'] <= ch1_band['max'])
                
                for ch2_idx, ch2_band in enumerate(ch2_bands):
                    ch2_mask = (df['Ch2_Amplitude'] >= ch2_band['min']) & (df['Ch2_Amplitude'] <= ch2_band['max'])
                    
                    combined_mask = ch1_mask & ch2_mask
                    band_data = df[combined_mask]
                    
                    if len(band_data) > 0:
                        color = ch1_colors[ch1_idx]
                        alpha = 0.3 + (0.4 * ch2_idx / max(len(ch2_bands) - 1, 1))
                        
                        ax.scatter(band_data['Ch1_Amplitude'], 
                                  band_data['Ch2_Amplitude'],
                                  c=color, alpha=alpha, s=0.5, edgecolors='none')
            
            # Draw separator lines
            for ch1_idx in range(len(ch1_bands) - 1):
                separator_x = (ch1_bands[ch1_idx]['max'] + ch1_bands[ch1_idx + 1]['min']) / 2
                ax.axvline(x=separator_x, color='red', linestyle='--', linewidth=0.3, alpha=0.5)
            
            for ch2_idx in range(len(ch2_bands) - 1):
                separator_y = (ch2_bands[ch2_idx]['max'] + ch2_bands[ch2_idx + 1]['min']) / 2
                ax.axhline(y=separator_y, color='red', linestyle='--', linewidth=0.3, alpha=0.5)
            
            # Format axes
            ax.set_xlim(ch1_min, ch1_max)
            ax.set_ylim(ch2_min, ch2_max)
            ax.set_xticks([])
            ax.set_yticks([])
            
            # Well label
            ax.text(0.5, 1.02, well_id, transform=ax.transAxes, 
                    fontsize=6, fontweight='bold', ha='center', va='bottom')
    
    # Add column headers
    for col in range(12):
        ax = plt.subplot(8, 12, col + 1)
        ax.text(0.5, 1.12, f'{col+1:02d}', transform=ax.transAxes,
               fontsize=7, ha='center', va='bottom', fontweight='bold')
    
    # Add row labels
    for row in range(8):
        ax = plt.subplot(8, 12, row * 12 + 1)
        row_letter = chr(ord('A') + row)
        ax.text(-0.15, 0.5, row_letter, transform=ax.transAxes,
               fontsize=7, ha='right', va='center', fontweight='bold')
    
    # Hide empty wells
    for row in range(8):
        for col in range(12):
            well_id = f'{chr(ord("A") + row)}{col+1:02d}'
            if well_id not in well_info:
                pos = row * 12 + col + 1
                ax = plt.subplot(8, 12, pos)
                ax.axis('off')

    plt.suptitle('96-Well Plate - 2D Scatter (Ch1 vs Ch2)\nRed lines: Band separators', 
                fontsize=13, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0.02, 0, 1, 0.99])
    pdf.savefig(fig, dpi=200, bbox_inches='tight')
    plt.close(fig)


def calculate_well_band_statistics(well_id, df):
    """Calculate band-based statistics for a single well"""
    total_droplets = len(df)
    
    # Detect bands for both channels
    ch1_bands = detect_bands_kmeans(df['Ch1_Amplitude'])
    ch2_bands = detect_bands_kmeans(df['Ch2_Amplitude'])
    
    stats_rows = []
    
    # Ch1 bands
    for band_idx, band in enumerate(ch1_bands):
        stats_rows.append({
            'Well': well_id,
            'Channel': 'Ch1_FAM',
            'Band_Number': band_idx + 1,
            'Band_Min': band['min'],
            'Band_Max': band['max'],
            'Band_Center': band['center'],
            'Droplet_Count': band['count'],
            'Proportion_Percent': band['proportion'],
            'Total_Droplets_In_Well': total_droplets
        })
    
    # Ch2 bands
    for band_idx, band in enumerate(ch2_bands):
        stats_rows.append({
            'Well': well_id,
            'Channel': 'Ch2_HEX',
            'Band_Number': band_idx + 1,
            'Band_Min': band['min'],
            'Band_Max': band['max'],
            'Band_Center': band['center'],
            'Droplet_Count': band['count'],
            'Proportion_Percent': band['proportion'],
            'Total_Droplets_In_Well': total_droplets
        })
    
    return stats_rows, ch1_bands, ch2_bands


def process_ddpcr_files(csv_files, log_func=print):
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
            output_pdf = output_dir / f"{well_id}_1D_plot.pdf"
            fig.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            # Create 2D scatter plot
            fig_2d = plot_well_2d_scatter(df, well_id, ch1_bands, ch2_bands)
            output_2d_pdf = output_dir / f"{well_id}_2D_plot.pdf"
            fig_2d.savefig(output_2d_pdf, format='pdf', dpi=300, bbox_inches='tight')
            plt.close(fig_2d)
            
            log_func(f"  Saved: {output_pdf.name} and {output_2d_pdf.name}")
            
        except Exception as e:
            log_func(f"  ⚠️  Error processing {csv_path.name}: {str(e)}")
            import traceback
            log_func(traceback.format_exc())
            continue
    
    # Create 96-well plate overview
    if all_data:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        plate_output = output_dir / f"96well_plate_overview_{timestamp}.pdf"
        plot_96well_layout(all_data, plate_output, log_func)

    stats_file = None
    if all_stats_rows:
        stats_df = pd.DataFrame(all_stats_rows)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        stats_file = output_dir / f"band_statistics_{timestamp}.csv"
        stats_df.to_csv(stats_file, index=False)
        log_func(f"\nBand statistics saved to: {stats_file.name}")
    
    return output_dir, stats_file


def process_qlp_file(qlp_file, well_assignments, include_filtered, log_func=print):
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
        
        # Calculate band statistics
        stats_rows, ch1_bands, ch2_bands = calculate_well_band_statistics(well_id, df)
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
        
        plot_path = output_dir / f"{well_id}_1D_plot.pdf"
        fig.savefig(plot_path, format='pdf', dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        # Create individual well 2D scatter plot
        fig_2d = plot_well_2d_scatter(df, well_id, ch1_bands, ch2_bands)
        
        # Add condition to title if assigned
        if condition != 'Unassigned':
            fig_2d.axes[0].set_title(f'Well {well_id} ({condition}) - 2D Scatter (Ch1 vs Ch2)', 
                                     fontsize=14, fontweight='bold')
        
        plot_2d_path = output_dir / f"{well_id}_2D_plot.pdf"
        fig_2d.savefig(plot_2d_path, format='pdf', dpi=300, bbox_inches='tight')
        plt.close(fig_2d)
    
    log_func(f"\nGenerated {len(well_dataframes)} individual well plots (1D and 2D)")
    
    # Create 96-well plate overview
    log_func(f"\nCreating 96-well plate overview...")
    plate_data = {}
    for well_id, df in well_dataframes.items():
        plate_data[well_id] = df
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    plate_output = output_dir / f"{qlp_path.stem}_96well_plate_{timestamp}.pdf"
    plot_96well_layout(plate_data, plate_output, log_func)
    
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
