#!/usr/bin/env python3
"""
ddPCRvis — ddPCR Visualization Tool
A GUI application for visualizing Bio-Rad ddPCR data from CSV exports,
QLP binary files, and encrypted ddPCR archives.

Data flow:
  Parse  → ddpcr_parsers.py   (BioRadQLPParser / BioRadDdpcrParser)
  Quant  → ddpcr_quant.py     (compute_well_quant — BioRad-matched Poisson)
  Plot   → ddPCRvis.py        (plot_well_1d / plot_well_2d / 96-well overview)

Clustering priority:
  1. BioRad cluster assignments from QLP/ddPCR files (used directly)
  2. Thresholds from file metadata (QLP thresholds embedded in file)
  3. Shared thresholds pooled across conditions (for CSV files)
  4. Auto-detected thresholds from amplitude distributions (CSV fallback)
"""

import sys
import os
import math
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
    from scipy.signal import find_peaks       # used by ddpcr_quant.detect_thresholds_1d
    from scipy.stats import gaussian_kde      # used by ddpcr_quant.detect_thresholds_1d
except ImportError:
    print("Installing scipy...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scipy", "--break-system-packages"])
    from scipy.signal import find_peaks
    from scipy.stats import gaussian_kde

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


# Full-featured QLP and ddPCR parsers live in ddpcr_parsers.py (drop in same folder)
from ddpcr_parsers import BioRadQLPParser, BioRadDdpcrParser
from ddpcr_quant import (
    WellQuant, ChannelQuant,
    compute_well_quant, well_quant_to_stats_row,
    CLUSTER_COLORS, cluster_labels,
    detect_thresholds_1d, assign_clusters_from_thresholds,
    pos_neg_from_labels,
)


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
                    self.output_format
                )
            elif self.file_type == 'ddpcr':
                output_dir, stats_file = process_ddpcr_file(
                    self.files[0],
                    self.well_assignments,
                    self.include_filtered,
                    self.progress.emit,
                    self.output_format
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
        ddpcr_files = []
        
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            if file_path.endswith('.csv'):
                csv_files.append(file_path)
            elif file_path.endswith('.qlp'):
                qlp_files.append(file_path)
            elif file_path.endswith('.ddpcr'):
                ddpcr_files.append(file_path)
        
        if qlp_files:
            if len(qlp_files) > 1:
                QMessageBox.warning(None, "Warning", "Only one QLP file can be processed at a time. Using the first file.")
            self.files_dropped.emit([qlp_files[0]], 'qlp')
        elif ddpcr_files:
            if len(ddpcr_files) > 1:
                QMessageBox.warning(None, "Warning", "Only one ddPCR file can be processed at a time. Using the first file.")
            self.files_dropped.emit([ddpcr_files[0]], 'ddpcr')
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
        subtitle = QLabel("Drag and drop Bio-Rad CSV files, QLP binary files, or ddPCR encrypted files")
        subtitle.setFont(QFont("Arial", 12))
        subtitle.setAlignment(Qt.AlignCenter)
        subtitle.setStyleSheet("color: #7f8c8d; margin-bottom: 20px;")
        layout.addWidget(subtitle)
        
        # Drop zone
        self.drop_zone = DropZone("📁 Drop CSV, QLP, or ddPCR files here\n\nor click a 'Browse' button below")
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
        
        browse_ddpcr_btn = QPushButton("📂 Browse ddPCR File")
        browse_ddpcr_btn.setFont(QFont("Arial", 11))
        browse_ddpcr_btn.setStyleSheet("""
            QPushButton {
                background-color: #e67e22;
                color: white;
                border: none;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #d35400;
            }
        """)
        browse_ddpcr_btn.clicked.connect(self.browse_ddpcr_file)
        
        browse_layout.addStretch()
        browse_layout.addWidget(browse_csv_btn)
        browse_layout.addWidget(browse_qlp_btn)
        browse_layout.addWidget(browse_ddpcr_btn)
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
        self.log("Drop CSV files, QLP files, ddPCR files, or click a 'Browse' button to get started.")
        self.log("- CSV files: Will be plotted individually and combined")
        self.log("- QLP files: Will prompt for well-to-condition assignments")
        self.log("- ddPCR files: Decrypted automatically, then same workflow as QLP")
        
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
    
    def browse_ddpcr_file(self):
        """Open file dialog to select a ddPCR file"""
        file, _ = QFileDialog.getOpenFileName(
            self,
            "Select ddPCR File",
            "",
            "ddPCR Files (*.ddpcr)"
        )
        if file:
            self.process_files([file], 'ddpcr')
    
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
        
        if file_type in ('qlp', 'ddpcr'):
            # Parse file to get well IDs, then show assignment dialog
            try:
                if file_type == 'qlp':
                    self.log("Parsing QLP file...")
                    parser = BioRadQLPParser(files[0])
                else:
                    self.log("Decrypting and parsing ddPCR file...")
                    parser = BioRadDdpcrParser(files[0])
                
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
                self.log(f"Error parsing file: {str(e)}")
                QMessageBox.critical(self, "Error", f"Failed to parse file:\n{str(e)}")
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




# ══════════════════════════════════════════════════════════════════════════════
#  Utilities
# ══════════════════════════════════════════════════════════════════════════════

def format_population_stats(count, total):
    """Format count and percentage; scientific notation for rare events."""
    if total == 0:
        return "0 (0%)"
    pct = count / total * 100
    if count == 0:
        pct_str = "0%"
    elif pct < 0.01:
        pct_str = f"{pct:.1e}%"
    else:
        pct_str = f"{pct:.1f}%"
    return f"{count:,} ({pct_str})"


def _amp_col(df: pd.DataFrame, name: str, fallback: str) -> str:
    """Resolve the amplitude column name for a given dye/channel."""
    for c in [f'{name}_Amplitude', fallback]:
        if c in df.columns:
            return c
    return fallback


def _shared_axis_limits(well_dataframes: dict, channel_names):
    """Compute amplitude axis limits across all wells for consistent scaling."""
    limits = {}
    for i, name in enumerate(channel_names[:2]):
        fallback = f'Ch{i+1}_Amplitude'
        vals = []
        for df in well_dataframes.values():
            col = _amp_col(df, name, fallback)
            if col in df.columns:
                v = df[col].values
                vals.extend(v[np.isfinite(v)])
        if vals:
            lo, hi = np.percentile(vals, 0.1), np.percentile(vals, 99.9)
            pad = (hi - lo) * 0.05
            limits[name] = (lo - pad, hi + pad)
    return limits


# ══════════════════════════════════════════════════════════════════════════════
#  1D amplitude plot  (per-channel strip, coloured by cluster)
# ══════════════════════════════════════════════════════════════════════════════

def plot_well_1d(df: pd.DataFrame, wq: WellQuant,
                  shared_limits: dict = None) -> 'plt.Figure':
    """
    1D amplitude strip plot, one panel per channel, coloured by cluster.

    Positive droplets for each channel are identified from the Cluster_Label
    column (contains '{dye}+'), so this works for both 2-channel duplex and
    N-channel amplitude-multiplex assays.  Gated droplets shown in light grey.
    Threshold line drawn where available from file metadata.
    """
    n_ch = len(wq.channels)
    if n_ch == 0:
        fig, ax = plt.subplots(figsize=(5, 7))
        ax.text(0.5, 0.5, 'No channels', ha='center', va='center',
                transform=ax.transAxes)
        return fig

    # Palette: index 0 = ch1 colour, index 1 = ch2 colour, etc.
    CH_COLORS = ['#1560bd', '#e06000', '#009944', '#c0392b',
                 '#8e44ad', '#16a085']

    fig, axes = plt.subplots(1, n_ch, figsize=(5 * n_ch, 7), squeeze=False)
    axes = axes[0]

    has_labels = ('Cluster_Label' in df.columns and
                  df['Cluster_Label'].notna().any())

    np.random.seed(42)
    max_pts = 20_000
    n = len(df)
    idx = np.random.choice(n, min(n, max_pts), replace=False)
    df_s = df.iloc[idx]

    for ax_i, (ax, ch) in enumerate(zip(axes, wq.channels)):
        col = _amp_col(df, ch.name, f'Ch{ax_i+1}_Amplitude')
        if col not in df.columns:
            ax.axis('off')
            continue

        amps = df_s[col].values
        ch_color = CH_COLORS[ax_i % len(CH_COLORS)]
        pos_pat = f'{ch.name}+'

        if has_labels:
            labels_s = df_s['Cluster_Label'].fillna('')
            gated_mask = (df_s['Cluster'] == 0 if 'Cluster' in df_s.columns
                          else labels_s.isin(['Gated', 'Filtered']))
            pos_mask = labels_s.str.contains(pos_pat, regex=False, na=False)
            neg_mask = ~pos_mask & ~gated_mask

            if gated_mask.any():
                jx = np.random.uniform(-0.3, 0.3, gated_mask.sum())
                ax.scatter(jx, amps[gated_mask.values], c='#e0e0e0', s=2,
                           alpha=0.35, edgecolors='none', rasterized=True, zorder=1)
            if neg_mask.any():
                jx = np.random.uniform(-0.3, 0.3, neg_mask.sum())
                ax.scatter(jx, amps[neg_mask.values], c='#909090', s=3,
                           alpha=0.5, edgecolors='none', rasterized=True, zorder=2)
            if pos_mask.any():
                jx = np.random.uniform(-0.3, 0.3, pos_mask.sum())
                ax.scatter(jx, amps[pos_mask.values], c=ch_color, s=3,
                           alpha=0.65, edgecolors='none', rasterized=True, zorder=3)
        else:
            jx = np.random.uniform(-0.3, 0.3, len(amps))
            ax.scatter(jx, amps, c='#909090', s=3, alpha=0.5,
                       edgecolors='none', rasterized=True)

        # Threshold line
        if ch.threshold is not None:
            ax.axhline(y=ch.threshold, color='red', linestyle='--',
                       linewidth=1.5, alpha=0.85, zorder=5)

        if shared_limits and ch.name in shared_limits:
            ax.set_ylim(shared_limits[ch.name])

        ax.set_xlim(-0.5, 0.5)
        ax.set_xticks([])
        ax.set_ylabel(f'{ch.name} Amplitude', fontsize=11)
        ax.grid(True, alpha=0.2, axis='y')

        # ── Title with concentration ──────────────────────────────────────
        title_lines = [ch.name]
        if ch.concentration is not None and np.isfinite(ch.concentration):
            title_lines.append(f'{ch.concentration:.2f} copies/µL')
            if ch.ci_lower_95 is not None and ch.ci_upper_95 is not None:
                title_lines.append(
                    f'95% CI  [{ch.ci_lower_95:.2f},  {ch.ci_upper_95:.2f}]')
        ax.set_title('\n'.join(title_lines), fontsize=10, fontweight='bold')

        # ── Stats box ────────────────────────────────────────────────────
        lines = [f'Positive: {ch.positives:,}',
                 f'Negative: {ch.negatives:,}']
        if wq.gated_count:
            lines.append(f'Gated:    {wq.gated_count:,}')
        lines.append(f'Accepted: {wq.accepted_count:,}')
        lines.append(f'Total:    {wq.total_events:,}')
        ax.text(0.97, 0.03, '\n'.join(lines),
                transform=ax.transAxes, fontsize=8, ha='right', va='bottom',
                family='monospace',
                bbox=dict(boxstyle='round', facecolor='white',
                          alpha=0.85, linewidth=0.5))

    fig.suptitle(f'Well {wq.well_id}', fontsize=14, fontweight='bold')
    plt.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  2D scatter plot  (Ch1 vs Ch2, coloured by cluster)
# ══════════════════════════════════════════════════════════════════════════════

def plot_well_2d(df: pd.DataFrame, wq: WellQuant,
                  shared_limits: dict = None) -> 'plt.Figure':
    """
    2D scatter plot(s) of channel pairs, coloured by Cluster_Label.

    For duplex (2 channels): one scatter panel.
    For multiplex (N channels): a grid of all N*(N-1)/2 pairwise scatters
      (capped at 6 panels to keep figure manageable).

    Cluster colours are assigned from the Cluster_Label string so every
    distinct population gets its own colour regardless of the integer encoding.
    Threshold lines shown where available from file metadata.
    """
    n_ch = len(wq.channels)
    if n_ch < 2:
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.text(0.5, 0.5, 'Need at least 2 channels',
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        return fig

    # Build all channel pairs (up to 6 panels)
    pairs = []
    for i in range(n_ch):
        for j in range(i + 1, n_ch):
            pairs.append((i, j))
            if len(pairs) >= 6:
                break
        if len(pairs) >= 6:
            break

    n_panels = len(pairs)
    cols = min(n_panels, 3)
    rows = math.ceil(n_panels / cols)
    fig, axes_arr = plt.subplots(rows, cols,
                                  figsize=(7 * cols, 7 * rows),
                                  squeeze=False)
    axes_flat = axes_arr.flatten()

    # Build a unique colour map from all distinct Cluster_Label values
    has_labels = ('Cluster_Label' in df.columns and
                  df['Cluster_Label'].notna().any())
    unique_labels = (df['Cluster_Label'].dropna().unique().tolist()
                     if has_labels else [])
    PALETTE = ['#808080',   # NN / generic negative
               '#1560bd',   # Ch1+ / FAM+
               '#e06000',   # Ch2+ / HEX+
               '#009944',   # Double+
               '#c0392b', '#8e44ad', '#16a085',
               '#f39c12', '#2980b9', '#d35400']
    label_color: Dict[str, str] = {}
    color_idx = 0
    for lbl in sorted(unique_labels):
        if lbl in ('Gated', 'Filtered'):
            label_color[lbl] = '#d0d0d0'
        elif '−' in lbl and '+' not in lbl:   # pure negative
            label_color[lbl] = '#808080'
        else:
            label_color[lbl] = PALETTE[color_idx % len(PALETTE)]
            color_idx += 1

    np.random.seed(42)
    max_pts = 30_000
    n = len(df)
    idx = np.random.choice(n, min(n, max_pts), replace=False)
    df_s = df.iloc[idx]

    for panel_i, (i, j) in enumerate(pairs):
        ax = axes_flat[panel_i]
        ch1 = wq.channels[i]
        ch2 = wq.channels[j]
        col1 = _amp_col(df, ch1.name, f'Ch{i+1}_Amplitude')
        col2 = _amp_col(df, ch2.name, f'Ch{j+1}_Amplitude')

        if col1 not in df.columns or col2 not in df.columns:
            ax.axis('off')
            continue

        # Plot by label (gated small and first so positives render on top)
        if has_labels:
            labels_s = df_s['Cluster_Label'].fillna('Unknown')
            # Sort: gated first, pure negatives next, positives last
            def _sort_key(lbl):
                if lbl in ('Gated', 'Filtered'):
                    return 0
                elif '+' not in lbl:
                    return 1
                return 2
            ordered = sorted(unique_labels, key=_sort_key)
            for lbl in ordered:
                mask = labels_s == lbl
                if not mask.any():
                    continue
                sub = df_s[mask]
                color = label_color.get(lbl, '#999999')
                is_gated = lbl in ('Gated', 'Filtered')
                ax.scatter(sub[col1], sub[col2],
                           c=color,
                           s=1.5 if is_gated else 3,
                           alpha=0.25 if is_gated else 0.55,
                           edgecolors='none', rasterized=True,
                           label=f'{lbl}: {wq.cluster_counts.get(_label_to_int(lbl, wq), 0):,}'
                                 if n_panels == 1 else None)
        else:
            ax.scatter(df_s[col1], df_s[col2],
                       c='#909090', s=3, alpha=0.5,
                       edgecolors='none', rasterized=True)

        # Threshold lines
        if ch1.threshold is not None:
            ax.axvline(x=ch1.threshold, color='red', linestyle='--',
                       linewidth=0.9, alpha=0.7)
        if ch2.threshold is not None:
            ax.axhline(y=ch2.threshold, color='red', linestyle='--',
                       linewidth=0.9, alpha=0.7)

        if shared_limits:
            if ch1.name in shared_limits:
                ax.set_xlim(shared_limits[ch1.name])
            if ch2.name in shared_limits:
                ax.set_ylim(shared_limits[ch2.name])

        ax.set_xlabel(f'{ch1.name} Amplitude', fontsize=11)
        ax.set_ylabel(f'{ch2.name} Amplitude', fontsize=11)
        ax.grid(True, alpha=0.2)

        # Panel title
        conc_parts = []
        for ch in [ch1, ch2]:
            if ch.concentration is not None and np.isfinite(ch.concentration):
                conc_parts.append(f'{ch.name}: {ch.concentration:.2f} cp/µL')
        panel_title = '    '.join(conc_parts) if conc_parts else f'{ch1.name} vs {ch2.name}'
        ax.set_title(panel_title, fontsize=10, fontweight='bold')

        # Legend only on the first panel (single panel or first of many)
        if panel_i == 0 and has_labels and n_panels == 1:
            ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
                      fontsize=9, borderaxespad=0, framealpha=0.9)

    # Hide unused axes
    for k in range(n_panels, len(axes_flat)):
        axes_flat[k].axis('off')

    # If >1 panel, add a shared legend at figure level
    if has_labels and n_panels > 1:
        handles = [plt.scatter([], [], c=label_color.get(lbl, '#999'),
                               s=20, label=lbl)
                   for lbl in sorted(unique_labels)]
        fig.legend(handles=handles, loc='lower center',
                   ncol=min(len(unique_labels), 4),
                   fontsize=9, framealpha=0.9,
                   bbox_to_anchor=(0.5, 0))
        fig.subplots_adjust(bottom=0.12)

    fig.suptitle(f'Well {wq.well_id}', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.08 if n_panels > 1 else 0, 1, 0.97])
    return fig


def _label_to_int(label: str, wq: WellQuant) -> int:
    """Map a Cluster_Label string back to its integer for legend count lookup."""
    for cid, cnt in wq.cluster_counts.items():
        pass  # can't reverse-map without the original df
    return -1   # fallback; counts shown inline for single-panel


# ══════════════════════════════════════════════════════════════════════════════
#  96-well plate overview
# ══════════════════════════════════════════════════════════════════════════════

def plot_96well_layout(data_dict, output_base_path, log_func=print,
                       condition_assignments=None, use_shared_thresholds=True,
                       output_format='png'):
    """Create 96-well plate visualisation (1D strip page + 2D scatter page)."""

    log_func("\nCreating 96-well plate visualisations...")

    # ── Compute per-well WellQuant objects ────────────────────────────────
    # If conditions were assigned, build pooled thresholds per condition
    # (only for CSV files; QLP/ddPCR already have BioRad clusters)
    shared_thresh_by_condition: dict = {}
    if use_shared_thresholds and condition_assignments:
        from collections import defaultdict
        cond_wells = defaultdict(list)
        for wid, cond in condition_assignments.items():
            if wid in data_dict:
                cond_wells[cond].append(wid)
        for cond, wids in cond_wells.items():
            first_df = data_dict[wids[0]]
            ch_names = first_df.attrs.get('channel_names', ['Ch1', 'Ch2'])
            # Only compute shared thresholds if none of the wells have BioRad clusters
            needs_thresh = all(
                not (data_dict[w]['Cluster'].isin([1, 2, 3, 4]).any()
                     if 'Cluster' in data_dict[w].columns else False)
                for w in wids
            )
            if needs_thresh:
                pooled_thresh = {}
                for i, name in enumerate(ch_names[:2]):
                    fallback = f'Ch{i+1}_Amplitude'
                    vals = []
                    for w in wids:
                        col = _amp_col(data_dict[w], name, fallback)
                        if col in data_dict[w].columns:
                            vals.extend(data_dict[w][col].values)
                    if vals:
                        detected = detect_thresholds_1d(np.array(vals))
                        if detected:
                            pooled_thresh[name] = detected[0]
                shared_thresh_by_condition[cond] = pooled_thresh

    well_quants = {}
    for key, df in data_dict.items():
        # Resolve well ID
        m = re.search(r'([A-Ha-h])(\d{1,2})', key)
        if m:
            wid = f'{m.group(1).upper()}{int(m.group(2)):02d}'
        else:
            wid = key
        cond = (condition_assignments or {}).get(wid, 'All')
        sh_thresh = shared_thresh_by_condition.get(cond)
        try:
            wq = compute_well_quant(df, shared_thresholds=sh_thresh)
            wq = WellQuant(
                well_id=wid,
                channels=wq.channels,
                cluster_counts=wq.cluster_counts,
                droplet_volume_nL=wq.droplet_volume_nL,
                total_events=wq.total_events,
                has_biorad_clusters=wq.has_biorad_clusters,
            )
            well_quants[wid] = (df, wq)
        except Exception as e:
            log_func(f"  ⚠️  Skipping {wid}: {e}")
            continue

    if not well_quants:
        log_func("  ⚠️  No valid wells to plot")
        return

    # Shared axis limits across all wells
    all_ch_names = next(iter(well_quants.values()))[0].attrs.get(
        'channel_names', ['Ch1', 'Ch2'])
    shared_lim = _shared_axis_limits(
        {wid: df for wid, (df, _) in well_quants.items()}, all_ch_names)

    output_base_path = Path(output_base_path)

    if output_format == 'pdf':
        pdf_path = output_base_path.with_suffix('.pdf')
        with PdfPages(pdf_path) as pdf:
            _create_96well_1d_page(well_quants, shared_lim, all_ch_names, pdf)
            _create_96well_2d_page(well_quants, shared_lim, all_ch_names, pdf)
        log_func(f"  ✓ 96-well plate saved: {pdf_path.name}")
    else:
        p1d = output_base_path.parent / f"{output_base_path.stem}_1D_view.png"
        p2d = output_base_path.parent / f"{output_base_path.stem}_2D_view.png"
        _create_96well_1d_page(well_quants, shared_lim, all_ch_names, str(p1d))
        _create_96well_2d_page(well_quants, shared_lim, all_ch_names, str(p2d))
        log_func(f"  ✓ 96-well plate saved: {p1d.name}  {p2d.name}")


def _create_96well_1d_page(well_quants, shared_lim, ch_names, destination):
    """96-well grid of 1D amplitude strips coloured by cluster."""
    ch1_name = ch_names[0] if ch_names else 'Ch1'
    ch2_name = ch_names[1] if len(ch_names) > 1 else 'Ch2'
    clabels = cluster_labels(ch1_name, ch2_name)

    fig = plt.figure(figsize=(22, 14))
    np.random.seed(42)

    # Global axis limits
    lim1 = shared_lim.get(ch1_name)
    lim2 = shared_lim.get(ch2_name)

    for row in range(8):
        for col in range(12):
            wid = f'{chr(ord("A") + row)}{col + 1:02d}'
            if wid not in well_quants:
                continue
            df, wq = well_quants[wid]

            pos_ch1 = row * 24 + col * 2 + 1
            pos_ch2 = row * 24 + col * 2 + 2
            ax1 = plt.subplot(8, 24, pos_ch1)
            ax2 = plt.subplot(8, 24, pos_ch2)

            has_cl = ('Cluster' in df.columns and
                      df['Cluster'].isin([1, 2, 3, 4]).any())

            ch1_col = _amp_col(df, ch1_name, 'Ch1_Amplitude')
            ch2_col = _amp_col(df, ch2_name, 'Ch2_Amplitude')

            n = len(df)
            max_pts = 5_000
            idx = np.random.choice(n, min(n, max_pts), replace=False)

            for ax, col_name, pos_set, color in [
                (ax1, ch1_col, {2, 3}, CLUSTER_COLORS[2]),
                (ax2, ch2_col, {3, 4}, CLUSTER_COLORS[4]),
            ]:
                if col_name not in df.columns:
                    ax.axis('off')
                    continue
                amps = df[col_name].values[idx]

                if has_cl:
                    cl = df['Cluster'].values[idx].astype(int)
                    g = cl == 0
                    if g.any():
                        jx = np.random.uniform(-0.4, 0.4, g.sum())
                        ax.scatter(jx, amps[g], c='#e8e8e8', s=0.8,
                                   alpha=0.5, edgecolors='none', rasterized=True)
                    neg = ~np.isin(cl, list(pos_set)) & (cl != 0)
                    if neg.any():
                        jx = np.random.uniform(-0.4, 0.4, neg.sum())
                        ax.scatter(jx, amps[neg], c='#b0b0b0', s=1,
                                   alpha=0.6, edgecolors='none', rasterized=True)
                    pos_m = np.isin(cl, list(pos_set))
                    if pos_m.any():
                        jx = np.random.uniform(-0.4, 0.4, pos_m.sum())
                        ax.scatter(jx, amps[pos_m], c=color, s=1.5,
                                   alpha=0.7, edgecolors='none', rasterized=True)
                else:
                    jx = np.random.uniform(-0.4, 0.4, len(idx))
                    ax.scatter(jx, amps, c='#b0b0b0', s=1,
                               alpha=0.6, edgecolors='none', rasterized=True)

            # Threshold lines
            ch1q = wq.ch(0)
            ch2q = wq.ch(1)
            if ch1q and ch1q.threshold is not None:
                ax1.axhline(y=ch1q.threshold, color='red',
                            linestyle='--', linewidth=0.5, alpha=0.8)
            if ch2q and ch2q.threshold is not None:
                ax2.axhline(y=ch2q.threshold, color='red',
                            linestyle='--', linewidth=0.5, alpha=0.8)

            # Axis formatting
            for ax, lim in [(ax1, lim1), (ax2, lim2)]:
                ax.set_xlim(-0.5, 0.5)
                ax.set_xticks([])
                ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                if lim:
                    ax.set_ylim(lim)

            # Well label
            ax1.text(0, 1.08, wid, transform=ax1.transAxes,
                     fontsize=5.5, fontweight='bold', ha='center', va='bottom')

            # Concentration text
            def _fmt_conc(ch_q):
                if ch_q and ch_q.concentration is not None:
                    return f'{ch_q.concentration:.1f}'
                pos = ch_q.positives if ch_q else 0
                neg = ch_q.negatives if ch_q else 0
                return f'P:{pos}'

            ch1q = wq.ch(0)
            ch2q = wq.ch(1)
            for ax, chq in [(ax1, ch1q), (ax2, ch2q)]:
                if chq:
                    txt = _fmt_conc(chq)
                    ax.text(0.5, 0.02, txt, transform=ax.transAxes,
                            fontsize=3.5, ha='center', va='bottom',
                            family='monospace', fontweight='bold',
                            bbox=dict(facecolor='white', alpha=0.7,
                                      linewidth=0, pad=0.5))

    # Column / row headers
    for c in range(12):
        ax = plt.subplot(8, 24, c * 2 + 1)
        ax.text(1, 1.18, f'{c+1:02d}', transform=ax.transAxes,
                fontsize=6.5, ha='center', va='bottom', fontweight='bold')
    for r in range(8):
        ax = plt.subplot(8, 24, r * 24 + 1)
        ax.text(-0.35, 0.5, chr(ord('A') + r), transform=ax.transAxes,
                fontsize=6.5, ha='right', va='center', fontweight='bold')
    for r in range(8):
        for c in range(12):
            wid = f'{chr(ord("A") + r)}{c+1:02d}'
            if wid not in well_quants:
                for offset in [1, 2]:
                    try:
                        plt.subplot(8, 24, r * 24 + c * 2 + offset).axis('off')
                    except Exception:
                        pass

    plt.suptitle(
        f'96-Well Plate — 1D Amplitude  '
        f'(Blue = {ch1_name}+,  Orange = {ch2_name}+,  '
        f'Text = copies/µL)',
        fontsize=12, fontweight='bold', y=0.998)
    plt.tight_layout(rect=[0.02, 0, 1, 0.998])

    if isinstance(destination, (str, Path)):
        fig.savefig(str(destination), dpi=250, bbox_inches='tight')
    else:
        destination.savefig(fig, dpi=200, bbox_inches='tight')
    plt.close(fig)


def _create_96well_2d_page(well_quants, shared_lim, ch_names, destination):
    """96-well grid of 2D scatter plots coloured by cluster."""
    ch1_name = ch_names[0] if ch_names else 'Ch1'
    ch2_name = ch_names[1] if len(ch_names) > 1 else 'Ch2'
    clabels = cluster_labels(ch1_name, ch2_name)

    fig = plt.figure(figsize=(22, 14))
    np.random.seed(42)

    lim1 = shared_lim.get(ch1_name)
    lim2 = shared_lim.get(ch2_name)

    for row in range(8):
        for col in range(12):
            wid = f'{chr(ord("A") + row)}{col + 1:02d}'
            if wid not in well_quants:
                continue
            df, wq = well_quants[wid]

            pos = row * 12 + col + 1
            ax = plt.subplot(8, 12, pos)

            ch1_col = _amp_col(df, ch1_name, 'Ch1_Amplitude')
            ch2_col = _amp_col(df, ch2_name, 'Ch2_Amplitude')

            if ch1_col not in df.columns or ch2_col not in df.columns:
                ax.axis('off')
                continue

            has_cl = ('Cluster' in df.columns and
                      df['Cluster'].isin([1, 2, 3, 4]).any())

            n = len(df)
            max_pts = 4_000
            idx = np.random.choice(n, min(n, max_pts), replace=False)
            df_p = df.iloc[idx]

            stats_parts = []
            for cid in [0, 1, 4, 2, 3]:
                if has_cl:
                    mask = df_p['Cluster'] == cid
                elif cid == 1:
                    mask = pd.Series(True, index=df_p.index)
                else:
                    continue
                if not mask.any():
                    continue
                sub = df_p[mask]
                ax.scatter(sub[ch1_col], sub[ch2_col],
                           c=CLUSTER_COLORS[cid],
                           s=0.8 if cid == 0 else 1.5,
                           alpha=0.35 if cid == 0 else 0.6,
                           edgecolors='none', rasterized=True)
                cnt = wq.cluster_counts.get(cid, 0)
                if cnt > 0 and cid != 0:
                    stats_parts.append(f'{clabels[cid].split()[0]}: {cnt:,}')

            # Threshold lines
            ch1q, ch2q = wq.ch(0), wq.ch(1)
            if ch1q and ch1q.threshold is not None:
                ax.axvline(x=ch1q.threshold, color='red',
                           linestyle='--', linewidth=0.3, alpha=0.6)
            if ch2q and ch2q.threshold is not None:
                ax.axhline(y=ch2q.threshold, color='red',
                           linestyle='--', linewidth=0.3, alpha=0.6)

            if lim1:
                ax.set_xlim(lim1)
            if lim2:
                ax.set_ylim(lim2)
            ax.set_xticks([])
            ax.set_yticks([])

            ax.text(0.5, 1.04, wid, transform=ax.transAxes,
                    fontsize=5.5, fontweight='bold', ha='center', va='bottom')

            # Concentration overlay
            conc_parts = []
            for chq in [ch1q, ch2q]:
                if chq and chq.concentration is not None:
                    conc_parts.append(f'{chq.name}: {chq.concentration:.1f}')
            if conc_parts:
                ax.text(0.97, 0.97, '\n'.join(conc_parts),
                        transform=ax.transAxes, fontsize=3.5,
                        ha='right', va='top', family='monospace',
                        bbox=dict(boxstyle='round', facecolor='white',
                                  alpha=0.8, linewidth=0, pad=0.3))

    for c in range(12):
        ax = plt.subplot(8, 12, c + 1)
        ax.text(0.5, 1.14, f'{c+1:02d}', transform=ax.transAxes,
                fontsize=6.5, ha='center', va='bottom', fontweight='bold')
    for r in range(8):
        ax = plt.subplot(8, 12, r * 12 + 1)
        ax.text(-0.18, 0.5, chr(ord('A') + r), transform=ax.transAxes,
                fontsize=6.5, ha='right', va='center', fontweight='bold')
    for r in range(8):
        for c in range(12):
            wid = f'{chr(ord("A") + r)}{c+1:02d}'
            if wid not in well_quants:
                try:
                    plt.subplot(8, 12, r * 12 + c + 1).axis('off')
                except Exception:
                    pass

    plt.suptitle(
        f'96-Well Plate — 2D Scatter  (copies/µL in corner)',
        fontsize=12, fontweight='bold', y=0.998)
    plt.tight_layout(rect=[0.02, 0, 1, 0.998])

    if isinstance(destination, (str, Path)):
        fig.savefig(str(destination), dpi=250, bbox_inches='tight')
    else:
        destination.savefig(fig, dpi=200, bbox_inches='tight')
    plt.close(fig)

def process_ddpcr_files(csv_files, log_func=print, output_format='png'):
    """Process Bio-Rad CSV exports and generate plots and statistics."""

    if not csv_files:
        raise ValueError("No CSV files provided")

    first_file = Path(csv_files[0])
    output_dir = first_file.parent / "plots"
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

            # Normalise column names
            if 'Ch1 Amplitude' in df.columns:
                df = df.rename(columns={'Ch1 Amplitude': 'Ch1_Amplitude',
                                         'Ch2 Amplitude': 'Ch2_Amplitude'})

            amp_cols = [c for c in df.columns if c.endswith('_Amplitude')]
            if len(amp_cols) >= 2:
                ch_names = [c.replace('_Amplitude', '') for c in amp_cols[:2]]
            else:
                ch_names = ['Ch1', 'Ch2']
            df.attrs['channel_names'] = ch_names
            df.attrs['channel_map'] = {f'Ch{i+1}_Amplitude': n
                                        for i, n in enumerate(ch_names)}

            m = re.search(r'([A-H])(\d{1,2})', csv_path.stem)
            well_id = (f'{m.group(1)}{int(m.group(2)):02d}' if m
                       else csv_path.stem)
            df.attrs.setdefault('well_metadata', {})['well_id'] = well_id

            all_data[well_id] = df

            wq = compute_well_quant(df)
            log_func(f"  Well {well_id}: "
                     + '  '.join(
                         f"{ch.name} {ch.concentration:.2f} copies/µL"
                         if ch.concentration is not None else f"{ch.name} pos={ch.positives}"
                         for ch in wq.channels))

            shared_lim = {}
            fig = plot_well_1d(df, wq, shared_lim)
            fig.savefig(output_dir / f"{well_id}_1D_plot.{output_format}",
                        format=output_format, dpi=300, bbox_inches='tight')
            plt.close(fig)

            fig2 = plot_well_2d(df, wq, shared_lim)
            fig2.savefig(output_dir / f"{well_id}_2D_plot.{output_format}",
                         format=output_format, dpi=300, bbox_inches='tight')
            plt.close(fig2)

            all_stats_rows.append(well_quant_to_stats_row(wq))

        except Exception as e:
            log_func(f"  ⚠️  Error processing {csv_path.name}: {e}")
            import traceback
            log_func(traceback.format_exc())
            continue

    if all_data:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        plate_output = output_dir / f"96well_plate_{timestamp}"
        plot_96well_layout(all_data, plate_output, log_func,
                           output_format=output_format)

    stats_file = None
    if all_stats_rows:
        stats_df = pd.DataFrame(all_stats_rows)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        stats_file = output_dir / f"statistics_{timestamp}.csv"
        stats_df.to_csv(stats_file, index=False)
        log_func(f"\nStatistics saved to: {stats_file.name}")

    return output_dir, stats_file


def _process_well_dataframes(well_dataframes, source_path, well_assignments,
                               include_filtered, log_func, output_format):
    """
    Shared processing pipeline for QLP and ddPCR files.
    Computes quantification, exports per-well CSVs, and generates plots.
    """
    output_dir = source_path.parent / "plots"
    output_dir.mkdir(exist_ok=True)
    csv_dir = output_dir / "well_csvs"
    csv_dir.mkdir(exist_ok=True)

    log_func(f"\nOutput directory: {output_dir}")
    log_func(f"Well CSV directory: {csv_dir}")

    # Build shared thresholds for wells that lack BioRad clusters (CSV fallback)
    # For QLP/ddPCR files with BioRad clusters, this is a no-op.
    shared_thresh_by_cond: dict = {}
    if well_assignments:
        from collections import defaultdict
        cond_wells = defaultdict(list)
        for wid, cond in well_assignments.items():
            if wid in well_dataframes:
                cond_wells[cond].append(wid)
        for cond, wids in cond_wells.items():
            first = well_dataframes[wids[0]]
            ch_names = first.attrs.get('channel_names', ['Ch1', 'Ch2'])
            needs = all(
                not (well_dataframes[w]['Cluster'].isin([1, 2, 3, 4]).any()
                     if 'Cluster' in well_dataframes[w].columns else False)
                for w in wids
            )
            if needs:
                pooled = {}
                for i, name in enumerate(ch_names[:2]):
                    fallback = f'Ch{i+1}_Amplitude'
                    vals = []
                    for w in wids:
                        col = _amp_col(well_dataframes[w], name, fallback)
                        if col in well_dataframes[w].columns:
                            vals.extend(well_dataframes[w][col].values)
                    if vals:
                        detected = detect_thresholds_1d(np.array(vals))
                        if detected:
                            pooled[name] = detected[0]
                shared_thresh_by_cond[cond] = pooled

    # Shared axis limits for consistent scaling across all wells
    first_df = next(iter(well_dataframes.values()))
    ch_names_global = first_df.attrs.get('channel_names', ['Ch1', 'Ch2'])
    shared_lim = _shared_axis_limits(well_dataframes, ch_names_global)

    all_stats_rows = []

    for well_id, df in sorted(well_dataframes.items()):
        ch_names = df.attrs.get('channel_names', ['Ch1', 'Ch2'])
        cond = (well_assignments or {}).get(well_id, 'Unassigned')
        sh_thresh = shared_thresh_by_cond.get(cond)

        # Export CSV
        df.to_csv(csv_dir / f"{well_id}.csv", index=False)

        # Compute quantification
        try:
            wq = compute_well_quant(df, shared_thresholds=sh_thresh)
        except Exception as e:
            log_func(f"  ⚠️  Quantification failed for {well_id}: {e}")
            continue

        # Log summary
        conc_str = '  '.join(
            f"{ch.name}: {ch.concentration:.2f} cp/µL"
            if ch.concentration is not None else f"{ch.name}: pos={ch.positives}"
            for ch in wq.channels)
        log_func(f"  {well_id} ({cond}): {wq.total_events:,} droplets  {conc_str}")

        # 1D plot
        fig = plot_well_1d(df, wq, shared_lim)
        if cond != 'Unassigned':
            fig.suptitle(f'Well {well_id}  ({cond})',
                         fontsize=14, fontweight='bold')
        fig.savefig(output_dir / f"{well_id}_1D_plot.{output_format}",
                    format=output_format, dpi=300, bbox_inches='tight')
        plt.close(fig)

        # 2D plot
        fig2 = plot_well_2d(df, wq, shared_lim)
        if cond != 'Unassigned':
            fig2.axes[0].set_title(
                f'Well {well_id}  ({cond})', fontsize=11, fontweight='bold')
        fig2.savefig(output_dir / f"{well_id}_2D_plot.{output_format}",
                     format=output_format, dpi=300, bbox_inches='tight')
        plt.close(fig2)

        row = well_quant_to_stats_row(wq, condition=cond)
        all_stats_rows.append(row)

    log_func(f"\nGenerated {len(well_dataframes)} individual well plots (1D and 2D)")

    # 96-well overview
    log_func("\nCreating 96-well plate overview...")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    plate_output = output_dir / f"{source_path.stem}_96well_plate_{timestamp}"
    plot_96well_layout(well_dataframes, plate_output, log_func,
                       condition_assignments=well_assignments,
                       use_shared_thresholds=bool(well_assignments),
                       output_format=output_format)

    # Statistics CSV
    stats_file = None
    if all_stats_rows:
        stats_df = pd.DataFrame(all_stats_rows)
        filtered_suffix = "_all" if include_filtered else "_accepted"
        stats_file = (output_dir /
                      f"{source_path.stem}_statistics{filtered_suffix}_{timestamp}.csv")
        stats_df.to_csv(stats_file, index=False)
        log_func(f"\nStatistics saved to: {stats_file.name}")
        log_func(f"Individual well CSVs in: {csv_dir.name}/")

    return output_dir, stats_file


def process_qlp_file(qlp_file, well_assignments, include_filtered,
                     log_func=print, output_format='png'):
    """Process a QLP binary file."""
    qlp_path = Path(qlp_file)
    log_func(f"\nParsing QLP file: {qlp_path.name}")
    log_func(f"Include filtered droplets (Cluster 0): {include_filtered}")

    parser = BioRadQLPParser(qlp_file)
    well_dataframes = parser.parse_to_dataframe(include_filtered=include_filtered)
    log_func(f"Extracted {len(well_dataframes)} wells")

    return _process_well_dataframes(well_dataframes, qlp_path, well_assignments,
                                     include_filtered, log_func, output_format)


def process_ddpcr_file(ddpcr_file, well_assignments, include_filtered,
                       log_func=print, output_format='png'):
    """Process a .ddpcr encrypted archive."""
    ddpcr_path = Path(ddpcr_file)
    log_func(f"\nDecrypting and parsing ddPCR file: {ddpcr_path.name}")
    log_func(f"Include filtered droplets (Cluster 0): {include_filtered}")

    parser = BioRadDdpcrParser(ddpcr_file)
    well_dataframes = parser.parse_to_dataframe(include_filtered=include_filtered)
    log_func(f"Extracted {len(well_dataframes)} wells")

    return _process_well_dataframes(well_dataframes, ddpcr_path, well_assignments,
                                     include_filtered, log_func, output_format)

def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()