#!/usr/bin/env python3
"""
Divine Pixel Tool - G-S Compatible Version with Cascading Divination and Multithreading
Advanced Image Pixel Divination System with G-S Alignment Data Support and Cascading Processing

Author: Angledcrystals
Date: 2025-06-10
Time: 10:33:03 UTC

This tool can "divine" pixels using the same alignment data format
as the G-S Divine Stereo Viewer, with full compatibility for
S_x, S_y, G_theta, G_phi, and Hadit vector data.
Now includes cascading divination where subsequent pixels are calculated
using previously calculated pixels.

MINIMAL THREADING: Added only essential threading support while preserving ALL original functionality.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import cv2
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy import ndimage
from scipy.interpolate import griddata, RectBivariateSpline
from scipy.spatial import cKDTree
import time
from datetime import datetime
import json
import sys
import traceback
from collections import deque
import threading
import os

class GSCompatibleDivinePixelTool:
    def __init__(self, root):
        self.root = root
        self.root.title("🔮 Divine Pixel Tool v1.3 - G-S Compatible with Cascading + Threading")
        self.root.geometry("1600x1000")
        
        # Image data
        self.original_image = None
        self.divined_image = None
        self.mask_image = None
        self.gs_alignment_data = None  # G-S compatible alignment data
        self.gs_coordinate_map = None  # Generated G-S coordinate maps
        
        # Threading support (MINIMAL ADDITION)
        self.processing_thread = None
        self.cancel_event = threading.Event()
        
        # Divination parameters (matching G-S viewer)
        self.divination_method = tk.StringVar(value="s_coordinate_flow")
        self.expansion_factor = tk.DoubleVar(value=1.5)
        self.divine_intensity = tk.DoubleVar(value=0.8)
        self.divine_precision = tk.DoubleVar(value=0.1)
        self.divine_search_radius = tk.IntVar(value=8)
        self.divine_similarity_threshold = tk.DoubleVar(value=10.0)
        self.edge_extension = tk.IntVar(value=50)
        
        # G-S specific parameters
        self.use_gs_alignment_data = tk.BooleanVar(value=False)
        self.gs_coordinate_precision = tk.DoubleVar(value=0.1)
        self.hadit_enhancement = tk.BooleanVar(value=False)
        
        # Advanced parameters
        self.use_flow_analysis = tk.BooleanVar(value=True)
        self.use_similarity_matching = tk.BooleanVar(value=True)
        self.use_geometric_extrapolation = tk.BooleanVar(value=True)
        self.use_gradient_synthesis = tk.BooleanVar(value=True)
        self.divine_outside_bounds = tk.BooleanVar(value=True)
        
        # Cascading divination parameters
        self.enable_cascading = tk.BooleanVar(value=True)
        self.cascading_order = tk.StringVar(value="nearest_first")
        self.cascading_batch_size = tk.IntVar(value=1)
        self.cascading_max_distance = tk.IntVar(value=50)
        self.cascading_update_frequency = tk.StringVar(value="immediate")
        self.show_cascading_progress = tk.BooleanVar(value=True)
        
        # Processing options
        self.noise_reduction = tk.BooleanVar(value=True)
        self.preserve_structures = tk.BooleanVar(value=True)
        self.auto_coordinate_generation = tk.BooleanVar(value=True)
        
        # Debug options
        self.debug_mode = tk.BooleanVar(value=False)
        self.show_divination_map = tk.BooleanVar(value=False)
        self.show_confidence_map = tk.BooleanVar(value=False)
        
        # Cascading tracking
        self.cascading_stats = {
            'total_pixels': 0,
            'processed_pixels': 0,
            'successful_pixels': 0,
            'cascade_generations': 0
        }
        
        self.setup_gui()
        
    def setup_gui(self):
        """Setup the main GUI layout."""
        # Control panel
        control_frame = ttk.Frame(self.root)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        
        # Visualization panel
        viz_frame = ttk.Frame(self.root)
        viz_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.setup_controls(control_frame)
        self.setup_visualization(viz_frame)
        
    def setup_controls(self, parent):
        """Setup control panel."""
        
        # Title
        title_label = ttk.Label(parent, text="🔮 Divine Pixel Tool", font=("Arial", 16, "bold"))
        title_label.pack(pady=(0, 10))
        
        subtitle_label = ttk.Label(parent, text="G-S Compatible with Cascading Divination + Threading", font=("Arial", 10))
        subtitle_label.pack(pady=(0, 15))
        
        # Create scrollable frame
        canvas = tk.Canvas(parent, highlightthickness=0, width=480)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Mouse wheel scrolling
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind_all("<MouseWheel>", _on_mousewheel)
        
        self.setup_control_sections(scrollable_frame)
        
    def setup_control_sections(self, parent):
        """Setup all control sections."""
        
        # Data Loading Section
        load_frame = ttk.LabelFrame(parent, text="Step 1: Load G-S Compatible Data", padding=10)
        load_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Button(load_frame, text="🖼️ Load Image to Divine", 
                  command=self.load_image, width=40).pack(fill=tk.X, pady=2)
        
        ttk.Button(load_frame, text="📊 Load G-S Alignment Data (JSON)", 
                  command=self.load_gs_alignment_data, width=40).pack(fill=tk.X, pady=2)
        
        ttk.Button(load_frame, text="🎭 Load Mask (Optional)", 
                  command=self.load_mask, width=40).pack(fill=tk.X, pady=2)
        
        ttk.Button(load_frame, text="🗺️ Load Coordinate Map (Optional)", 
                  command=self.load_coordinate_map, width=40).pack(fill=tk.X, pady=2)
        
        self.image_status = ttk.Label(load_frame, text="No image loaded", foreground="red")
        self.image_status.pack(pady=5)
        
        self.alignment_status = ttk.Label(load_frame, text="No G-S alignment data loaded", foreground="orange")
        self.alignment_status.pack(pady=2)
        
        # G-S Alignment Section
        gs_frame = ttk.LabelFrame(parent, text="Step 2: G-S Alignment Configuration", padding=10)
        gs_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Checkbutton(gs_frame, text="🔗 Use G-S Alignment Data for Enhanced Divination", 
                       variable=self.use_gs_alignment_data).pack(anchor=tk.W, pady=2)
        
        ttk.Label(gs_frame, text="G-S Coordinate Precision:").pack(anchor=tk.W)
        gs_precision_frame = ttk.Frame(gs_frame)
        gs_precision_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(gs_precision_frame, from_=0.01, to=1.0, variable=self.gs_coordinate_precision, 
                 orient=tk.HORIZONTAL).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(gs_precision_frame, textvariable=self.gs_coordinate_precision, width=8).pack(side=tk.RIGHT)
        
        ttk.Checkbutton(gs_frame, text="🌟 Enable Hadit Vector Enhancement", 
                       variable=self.hadit_enhancement).pack(anchor=tk.W, pady=2)
        
        self.generate_gs_coords_button = ttk.Button(gs_frame, text="🧮 Generate G-S Coordinate Maps", 
                                                   command=self.generate_gs_coordinate_maps, 
                                                   state="disabled")
        self.generate_gs_coords_button.pack(fill=tk.X, pady=5)
        
        # Cascading Divination Section
        cascade_frame = ttk.LabelFrame(parent, text="Step 2.5: Cascading Divination", padding=10)
        cascade_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Checkbutton(cascade_frame, text="🌊 Enable Cascading Divination", 
                       variable=self.enable_cascading).pack(anchor=tk.W, pady=2)
        
        ttk.Label(cascade_frame, text="Cascading Order:").pack(anchor=tk.W)
        cascade_order_combo = ttk.Combobox(cascade_frame, textvariable=self.cascading_order,
                                          values=["nearest_first", "flow_based", "coordinate_based", 
                                                 "spiral_out", "random", "distance_weighted"])
        cascade_order_combo.pack(fill=tk.X, pady=2)
        
        ttk.Label(cascade_frame, text="Batch Size (pixels per update):").pack(anchor=tk.W)
        batch_frame = ttk.Frame(cascade_frame)
        batch_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(batch_frame, from_=1, to=20, variable=self.cascading_batch_size, 
                 orient=tk.HORIZONTAL).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(batch_frame, textvariable=self.cascading_batch_size, width=8).pack(side=tk.RIGHT)
        
        ttk.Label(cascade_frame, text="Max Cascade Distance:").pack(anchor=tk.W)
        distance_frame = ttk.Frame(cascade_frame)
        distance_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(distance_frame, from_=10, to=100, variable=self.cascading_max_distance, 
                 orient=tk.HORIZONTAL).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(distance_frame, textvariable=self.cascading_max_distance, width=8).pack(side=tk.RIGHT)
        
        ttk.Label(cascade_frame, text="Update Frequency:").pack(anchor=tk.W)
        update_combo = ttk.Combobox(cascade_frame, textvariable=self.cascading_update_frequency,
                                   values=["immediate", "batch", "generation"])
        update_combo.pack(fill=tk.X, pady=2)
        
        ttk.Checkbutton(cascade_frame, text="📊 Show Cascading Progress", 
                       variable=self.show_cascading_progress).pack(anchor=tk.W, pady=2)
        
        # Divination Method Section
        method_frame = ttk.LabelFrame(parent, text="Step 3: Divination Method", padding=10)
        method_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(method_frame, text="Primary Divination Algorithm:").pack(anchor=tk.W)
        method_combo = ttk.Combobox(method_frame, textvariable=self.divination_method,
                                   values=["s_coordinate_flow", "s_coordinate_similarity", 
                                          "geometric_extrapolation", "combined_divine", 
                                          "gs_enhanced_divine", "omniscient_divine",
                                          "cascading_divine"])
        method_combo.pack(fill=tk.X, pady=2)
        
        # Divination Parameters (matching G-S viewer)
        params_frame = ttk.LabelFrame(parent, text="Step 4: Divination Parameters", padding=10)
        params_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Divine Precision (matching G-S naming)
        ttk.Label(params_frame, text="Divine Precision (S-coordinate sensitivity):").pack(anchor=tk.W)
        precision_frame = ttk.Frame(params_frame)
        precision_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(precision_frame, from_=0.01, to=1.0, variable=self.divine_precision, 
                 orient=tk.HORIZONTAL).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(precision_frame, textvariable=self.divine_precision, width=8).pack(side=tk.RIGHT)
        
        # Search Radius (matching G-S naming)
        ttk.Label(params_frame, text="Search Radius (pixels):").pack(anchor=tk.W)
        search_frame = ttk.Frame(params_frame)
        search_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(search_frame, from_=3, to=20, variable=self.divine_search_radius, 
                 orient=tk.HORIZONTAL).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(search_frame, textvariable=self.divine_search_radius, width=8).pack(side=tk.RIGHT)
        
        # Similarity Threshold (matching G-S naming)
        ttk.Label(params_frame, text="Similarity Threshold (degrees):").pack(anchor=tk.W)
        similarity_frame = ttk.Frame(params_frame)
        similarity_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(similarity_frame, from_=1.0, to=30.0, variable=self.divine_similarity_threshold, 
                 orient=tk.HORIZONTAL).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(similarity_frame, textvariable=self.divine_similarity_threshold, width=8).pack(side=tk.RIGHT)
        
        # Expansion Factor
        ttk.Label(params_frame, text="Expansion Factor (how far to divine):").pack(anchor=tk.W)
        expansion_frame = ttk.Frame(params_frame)
        expansion_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(expansion_frame, from_=1.0, to=3.0, variable=self.expansion_factor, 
                 orient=tk.HORIZONTAL).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(expansion_frame, textvariable=self.expansion_factor, width=8).pack(side=tk.RIGHT)
        
        # Edge Extension
        ttk.Label(params_frame, text="Edge Extension (pixels):").pack(anchor=tk.W)
        edge_frame = ttk.Frame(params_frame)
        edge_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(edge_frame, from_=10, to=200, variable=self.edge_extension, 
                 orient=tk.HORIZONTAL).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(edge_frame, textvariable=self.edge_extension, width=8).pack(side=tk.RIGHT)
        
        # Advanced Divination Options
        advanced_frame = ttk.LabelFrame(parent, text="Advanced Divination Methods", padding=10)
        advanced_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Checkbutton(advanced_frame, text="🌊 S-Coordinate Flow Analysis", 
                       variable=self.use_flow_analysis).pack(anchor=tk.W, pady=1)
        
        ttk.Checkbutton(advanced_frame, text="🔍 G-S Similarity Matching", 
                       variable=self.use_similarity_matching).pack(anchor=tk.W, pady=1)
        
        ttk.Checkbutton(advanced_frame, text="📐 Geometric Extrapolation", 
                       variable=self.use_geometric_extrapolation).pack(anchor=tk.W, pady=1)
        
        ttk.Checkbutton(advanced_frame, text="🧮 Gradient Synthesis", 
                       variable=self.use_gradient_synthesis).pack(anchor=tk.W, pady=1)
        
        ttk.Checkbutton(advanced_frame, text="🚀 Divine Outside Image Bounds", 
                       variable=self.divine_outside_bounds).pack(anchor=tk.W, pady=1)
        
        # Processing Options
        process_frame = ttk.LabelFrame(parent, text="Processing Options", padding=10)
        process_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Checkbutton(process_frame, text="🧹 Noise Reduction", 
                       variable=self.noise_reduction).pack(anchor=tk.W, pady=1)
        
        ttk.Checkbutton(process_frame, text="🏗️ Preserve Structures", 
                       variable=self.preserve_structures).pack(anchor=tk.W, pady=1)
        
        ttk.Checkbutton(process_frame, text="🗺️ Auto-Generate Coordinates", 
                       variable=self.auto_coordinate_generation).pack(anchor=tk.W, pady=1)
        
        # Debug Options
        debug_frame = ttk.LabelFrame(parent, text="Debug & Visualization", padding=10)
        debug_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Checkbutton(debug_frame, text="Debug Divine Process", 
                       variable=self.debug_mode).pack(anchor=tk.W, pady=1)
        
        ttk.Checkbutton(debug_frame, text="Show Divination Map", 
                       variable=self.show_divination_map).pack(anchor=tk.W, pady=1)
        
        ttk.Checkbutton(debug_frame, text="Show Confidence Map", 
                       variable=self.show_confidence_map).pack(anchor=tk.W, pady=1)
        
        # Divine Button (MODIFIED FOR THREADING)
        divine_frame = ttk.LabelFrame(parent, text="Step 5: Perform Divination", padding=10)
        divine_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.divine_button = ttk.Button(divine_frame, text="🔮 DIVINE PIXELS (THREADED)", 
                                       command=self.divine_pixels_threaded, 
                                       state="disabled")
        self.divine_button.pack(fill=tk.X, pady=5)
        
        # NEW: Cancel button for threading
        self.cancel_button = ttk.Button(divine_frame, text="⏹️ Cancel Operation", 
                                       command=self.cancel_processing, 
                                       state="disabled")
        self.cancel_button.pack(fill=tk.X, pady=2)
        
        self.test_button = ttk.Button(divine_frame, text="🧪 Test Divine Process", 
                                     command=self.test_divine_process, 
                                     state="disabled")
        self.test_button.pack(fill=tk.X, pady=2)
        
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(divine_frame, variable=self.progress_var, 
                                          maximum=100)
        self.progress_bar.pack(fill=tk.X, pady=5)
        
        self.status_label = ttk.Label(divine_frame, text="Load an image to begin", 
                                     foreground="orange")
        self.status_label.pack(pady=5)
        
        # Cascading Progress
        self.cascade_progress_var = tk.DoubleVar()
        self.cascade_progress_bar = ttk.Progressbar(divine_frame, variable=self.cascade_progress_var, 
                                                   maximum=100)
        self.cascade_progress_bar.pack(fill=tk.X, pady=2)
        
        self.cascade_status_label = ttk.Label(divine_frame, text="", 
                                             foreground="blue", font=("Arial", 8))
        self.cascade_status_label.pack(pady=2)
        
        # Export Section
        export_frame = ttk.LabelFrame(parent, text="Step 6: Export Results", padding=10)
        export_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.export_button = ttk.Button(export_frame, text="💾 Export Divined Image", 
                                       command=self.export_result, 
                                       state="disabled")
        self.export_button.pack(fill=tk.X, pady=2)
        
        self.export_maps_button = ttk.Button(export_frame, text="🗺️ Export G-S Divine Maps", 
                                            command=self.export_maps, 
                                            state="disabled")
        self.export_maps_button.pack(fill=tk.X, pady=2)
        
        # Statistics Section
        stats_frame = ttk.LabelFrame(parent, text="Divine Enhancement Statistics", padding=5)
        stats_frame.pack(fill=tk.X, pady=(10, 0))
        
        self.stats_label = ttk.Label(stats_frame, text="No divination performed yet", 
                                    font=("Courier", 8))
        self.stats_label.pack(fill=tk.X)

    def setup_visualization(self, parent):
        """Setup visualization panel."""
        self.fig = Figure(figsize=(14, 10), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, parent)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Toolbar
        toolbar_frame = ttk.Frame(parent)
        toolbar_frame.pack(fill=tk.X)
        
        from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()
        
        self.show_welcome_display()
        
    def show_welcome_display(self):
        """Show welcome message."""
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        
        welcome_text = """🔮 DIVINE PIXEL TOOL v1.3 + THREADING
G-S Compatible Advanced Image Pixel Divination
with Cascading Processing + Background Threading
by Angledcrystals

🌟 G-S Compatible Features:
• Full compatibility with G-S Stereo Viewer alignment data
• S_x, S_y coordinate support 
• G_theta, G_phi coordinate integration
• Hadit_theta, Hadit_phi vector enhancement
• S-coordinate flow-based pixel divination
• Geometric extrapolation beyond image bounds

🌊 Cascading Divination Features:
• Sequential pixel processing using previously calculated pixels
• Multiple cascading orders: nearest-first, flow-based, coordinate-based
• Configurable batch processing and update frequencies
• Real-time cascading progress visualization
• Enhanced continuity and pattern propagation

🧵 NEW: Threading Support:
• Background processing keeps GUI responsive
• Cancellable operations
• Real-time progress updates
• Non-blocking interface during intensive processing

📊 Cascading Order Options:
• Nearest-First: Process pixels closest to known data first
• Flow-Based: Follow S-coordinate flow directions
• Coordinate-Based: Process by coordinate similarity
• Spiral-Out: Expand from center in spiral pattern
• Distance-Weighted: Prioritize by multiple distance metrics

Load an image and configure cascading parameters to unlock
the enhanced G-S compatible threaded cascading divine pixel experience!"""
        
        ax.text(0.5, 0.5, welcome_text, ha='center', va='center', 
                fontsize=9, transform=ax.transAxes, 
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightcyan", alpha=0.9))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        ax.set_title("G-S Compatible Divine Pixel Tool with Cascading + Threading", fontsize=16, weight='bold')
        
        self.canvas.draw()

    # === THREADING WRAPPER FUNCTIONS (MINIMAL ADDITION) ===
    
    def divine_pixels_threaded(self):
        """Threaded wrapper for the original divine_pixels function."""
        if self.original_image is None:
            messagebox.showwarning("No Data", "Load an image first")
            return
        
        # Disable divine button, enable cancel
        self.divine_button.config(state="disabled")
        self.cancel_button.config(state="normal")
        self.cancel_event.clear()
        
        # Start processing in background thread
        self.processing_thread = threading.Thread(target=self.divine_pixels_worker, daemon=True)
        self.processing_thread.start()
    
    def divine_pixels_worker(self):
        """Background worker thread that runs the original divine_pixels logic."""
        try:
            # Call the ORIGINAL divine_pixels method with threading-safe updates
            self.divine_pixels()
        except Exception as e:
            # Thread-safe error handling
            self.root.after(0, lambda: self.handle_processing_error(str(e)))
        finally:
            # Thread-safe cleanup
            self.root.after(0, self.finish_processing)
    
    def cancel_processing(self):
        """Cancel the current processing operation."""
        self.cancel_event.set()
        self.root.after(0, self.finish_processing)
        self.root.after(0, lambda: self.status_label.config(text="❌ Operation cancelled", foreground="red"))
    
    def finish_processing(self):
        """Re-enable controls after processing finishes."""
        self.divine_button.config(state="normal")
        self.cancel_button.config(state="disabled")
    
    def handle_processing_error(self, error_msg):
        """Handle processing errors in thread-safe manner."""
        self.status_label.config(text=f"❌ Error: {error_msg}", foreground="red")
        self.finish_processing()
        messagebox.showerror("Processing Error", f"Divination failed: {error_msg}")
    
    def update_progress_threaded(self, value):
        """Thread-safe progress update."""
        self.root.after(0, lambda: self.progress_var.set(value))
    
    def update_cascade_progress_threaded(self, value, text):
        """Thread-safe cascade progress update."""
        self.root.after(0, lambda: self.cascade_progress_var.set(value))
        self.root.after(0, lambda: self.cascade_status_label.config(text=text))
    
    def update_status_threaded(self, text, color="blue"):
        """Thread-safe status update."""
        self.root.after(0, lambda: self.status_label.config(text=text, foreground=color))

    # === ALL ORIGINAL METHODS PRESERVED EXACTLY ===
    
    def cascading_divine_process(self, image, mask):
        """Main cascading divination process where subsequent pixels use previously calculated pixels."""
        if self.debug_mode.get():
            print("🌊 Starting cascading divine process...")
        
        result = image.copy()
        remaining_mask = mask.copy()
        
        # Initialize cascading statistics
        self.cascading_stats = {
            'total_pixels': np.sum(mask),
            'processed_pixels': 0,
            'successful_pixels': 0,
            'cascade_generations': 0
        }
        
        if self.cascading_stats['total_pixels'] == 0:
            return result
        
        # Order pixels for cascading processing
        ordered_coords = self.order_pixels_for_cascading(remaining_mask, self.cascading_order.get())
        
        if self.debug_mode.get():
            print(f"   Ordered {len(ordered_coords)} pixels for cascading")
        
        batch_size = self.cascading_batch_size.get()
        update_frequency = self.cascading_update_frequency.get()
        
        # Process pixels in cascading order
        generation = 0
        total_processed = 0
        
        while len(ordered_coords) > 0 and generation < 1000 and not self.cancel_event.is_set():  # Safety limit + cancel check
            generation += 1
            batch_success = 0
            batch_coords = []
            
            # Process current batch
            for i in range(min(batch_size, len(ordered_coords))):
                if self.cancel_event.is_set():  # Check for cancellation
                    break
                    
                coord_idx, coord_y, coord_x = ordered_coords.pop(0)
                
                # Skip if pixel was already processed in this generation
                if not remaining_mask[coord_y, coord_x]:
                    continue
                
                # Divine this pixel using current state of the image
                divined_pixel = self.divine_single_pixel_cascading(
                    coord_y, coord_x, result, remaining_mask)
                
                if divined_pixel is not None:
                    result[coord_y, coord_x] = divined_pixel
                    remaining_mask[coord_y, coord_x] = False
                    batch_success += 1
                    batch_coords.append((coord_y, coord_x))
                    self.cascading_stats['successful_pixels'] += 1
                
                self.cascading_stats['processed_pixels'] += 1
                total_processed += 1
                
                # Update progress for immediate mode (THREAD-SAFE)
                if update_frequency == "immediate" and self.show_cascading_progress.get():
                    progress = (self.cascading_stats['processed_pixels'] / 
                              self.cascading_stats['total_pixels']) * 100
                    status_text = f"Gen {generation}: {self.cascading_stats['successful_pixels']}/{self.cascading_stats['total_pixels']} pixels"
                    self.update_cascade_progress_threaded(progress, status_text)
            
            # Early exit if cancelled
            if self.cancel_event.is_set():
                break
            
            # Update coordinates for next generation if any pixels were successfully processed
            if batch_success > 0 and len(ordered_coords) > 0:
                # Re-order remaining pixels based on new state
                if update_frequency == "batch" or (generation % 5 == 0):  # Re-order every few generations
                    remaining_coords_mask = np.zeros_like(remaining_mask)
                    for _, coord_y, coord_x in ordered_coords:
                        if remaining_mask[coord_y, coord_x]:
                            remaining_coords_mask[coord_y, coord_x] = True
                    
                    if np.any(remaining_coords_mask):
                        ordered_coords = self.order_pixels_for_cascading(
                            remaining_coords_mask, self.cascading_order.get())
            
            self.cascading_stats['cascade_generations'] = generation
            
            # Update progress for batch/generation modes (THREAD-SAFE)
            if update_frequency in ["batch", "generation"] and self.show_cascading_progress.get():
                progress = (self.cascading_stats['processed_pixels'] / 
                          self.cascading_stats['total_pixels']) * 100
                status_text = f"Gen {generation}: {self.cascading_stats['successful_pixels']}/{self.cascading_stats['total_pixels']} pixels"
                self.update_cascade_progress_threaded(progress, status_text)
            
            if self.debug_mode.get():
                print(f"   Generation {generation}: {batch_success} pixels processed, {len(ordered_coords)} remaining")
            
            # Stop if no progress in this generation
            if batch_success == 0:
                break
        
        if self.debug_mode.get():
            print(f"   Cascading complete: {self.cascading_stats['successful_pixels']}/{self.cascading_stats['total_pixels']} pixels successful in {generation} generations")
        
        return result

    def order_pixels_for_cascading(self, mask, order_method):
        """Order pixels based on cascading strategy."""
        mask_coords = np.column_stack(np.where(mask))
        if len(mask_coords) == 0:
            return []
        
        height, width = mask.shape
        
        if order_method == "nearest_first":
            return self.order_nearest_first(mask_coords, mask)
        elif order_method == "flow_based":
            return self.order_flow_based(mask_coords, mask)
        elif order_method == "coordinate_based":
            return self.order_coordinate_based(mask_coords, mask)
        elif order_method == "spiral_out":
            return self.order_spiral_out(mask_coords, mask)
        elif order_method == "distance_weighted":
            return self.order_distance_weighted(mask_coords, mask)
        elif order_method == "random":
            np.random.shuffle(mask_coords)
            return [(i, y, x) for i, (y, x) in enumerate(mask_coords)]
        else:
            # Default to nearest first
            return self.order_nearest_first(mask_coords, mask)

    def order_nearest_first(self, mask_coords, mask):
        """Order pixels by distance to nearest known (non-masked) pixel."""
        height, width = mask.shape
        
        # Find all known pixels (not masked)
        known_mask = ~mask
        if not np.any(known_mask):
            # If no known pixels, start from center
            center_y, center_x = height // 2, width // 2
            distances = np.sqrt((mask_coords[:, 0] - center_y)**2 + (mask_coords[:, 1] - center_x)**2)
        else:
            known_coords = np.column_stack(np.where(known_mask))
            
            # Calculate minimum distance to any known pixel
            distances = []
            for y, x in mask_coords:
                min_dist = np.min(np.sqrt((known_coords[:, 0] - y)**2 + (known_coords[:, 1] - x)**2))
                distances.append(min_dist)
            distances = np.array(distances)
        
        # Sort by distance (nearest first)
        sorted_indices = np.argsort(distances)
        ordered_coords = [(i, mask_coords[sorted_indices[i]][0], mask_coords[sorted_indices[i]][1]) 
                         for i in range(len(mask_coords))]
        
        return ordered_coords

    def order_flow_based(self, mask_coords, mask):
        """Order pixels based on S-coordinate flow direction."""
        if self.gs_coordinate_map is None:
            return self.order_nearest_first(mask_coords, mask)
        
        height, width = mask.shape
        
        # Get coordinate maps
        if 's_x' in self.gs_coordinate_map and 's_y' in self.gs_coordinate_map:
            s_x = self.gs_coordinate_map['s_x']
            s_y = self.gs_coordinate_map['s_y']
        else:
            return self.order_nearest_first(mask_coords, mask)
        
        # Resize if needed
        if s_x.shape != (height, width):
            s_x = cv2.resize(s_x, (width, height), interpolation=cv2.INTER_LINEAR)
            s_y = cv2.resize(s_y, (width, height), interpolation=cv2.INTER_LINEAR)
        
        # Calculate flow vectors
        s_x_grad_x = np.gradient(s_x, axis=1)
        s_x_grad_y = np.gradient(s_x, axis=0)
        s_y_grad_x = np.gradient(s_y, axis=1)
        s_y_grad_y = np.gradient(s_y, axis=0)
        
        flow_x = s_x_grad_x + s_y_grad_x
        flow_y = s_x_grad_y + s_y_grad_y
        flow_magnitude = np.sqrt(flow_x**2 + flow_y**2)
        
        # Order by flow magnitude (highest flow first)
        flow_values = []
        for y, x in mask_coords:
            if y < height and x < width:
                flow_values.append(flow_magnitude[y, x])
            else:
                flow_values.append(0)
        
        flow_values = np.array(flow_values)
        sorted_indices = np.argsort(-flow_values)  # Descending order
        
        ordered_coords = [(i, mask_coords[sorted_indices[i]][0], mask_coords[sorted_indices[i]][1]) 
                         for i in range(len(mask_coords))]
        
        return ordered_coords

    def order_coordinate_based(self, mask_coords, mask):
        """Order pixels based on coordinate similarity to existing pixels."""
        if self.gs_coordinate_map is None:
            return self.order_nearest_first(mask_coords, mask)
        
        height, width = mask.shape
        
        # Get coordinate maps
        if 's_x' in self.gs_coordinate_map and 's_y' in self.gs_coordinate_map:
            s_x = self.gs_coordinate_map['s_x']
            s_y = self.gs_coordinate_map['s_y']
        else:
            return self.order_nearest_first(mask_coords, mask)
        
        # Resize if needed
        if s_x.shape != (height, width):
            s_x = cv2.resize(s_x, (width, height), interpolation=cv2.INTER_LINEAR)
            s_y = cv2.resize(s_y, (width, height), interpolation=cv2.INTER_LINEAR)
        
        # Find known pixels
        known_mask = ~mask
        if not np.any(known_mask):
            return self.order_nearest_first(mask_coords, mask)
        
        known_coords = np.column_stack(np.where(known_mask))
        known_s_coords = np.column_stack([s_x[known_mask], s_y[known_mask]])
        
        # Build KD-tree for coordinate lookups
        from scipy.spatial import cKDTree
        coord_tree = cKDTree(known_s_coords)
        
        # Calculate coordinate distances for mask pixels
        coord_distances = []
        for y, x in mask_coords:
            if y < height and x < width:
                pixel_s_coord = np.array([s_x[y, x], s_y[y, x]])
                min_dist, _ = coord_tree.query(pixel_s_coord)
                coord_distances.append(min_dist)
            else:
                coord_distances.append(float('inf'))
        
        coord_distances = np.array(coord_distances)
        sorted_indices = np.argsort(coord_distances)  # Ascending order (closest coordinates first)
        
        ordered_coords = [(i, mask_coords[sorted_indices[i]][0], mask_coords[sorted_indices[i]][1]) 
                         for i in range(len(mask_coords))]
        
        return ordered_coords

    def order_spiral_out(self, mask_coords, mask):
        """Order pixels in expanding spiral from center."""
        height, width = mask.shape
        center_y, center_x = height // 2, width // 2
        
        # Calculate angle and distance from center
        angles_distances = []
        for y, x in mask_coords:
            dy, dx = y - center_y, x - center_x
            angle = np.arctan2(dy, dx)
            distance = np.sqrt(dy**2 + dx**2)
            angles_distances.append((angle, distance))
        
        # Sort by distance first, then by angle for spiral effect
        angles_distances = np.array(angles_distances)
        distances = angles_distances[:, 1]
        angles = angles_distances[:, 0]
        
        # Create spiral ordering: sort by distance, then by angle within distance bands
        distance_order = np.argsort(distances)
        
        ordered_coords = [(i, mask_coords[distance_order[i]][0], mask_coords[distance_order[i]][1]) 
                         for i in range(len(mask_coords))]
        
        return ordered_coords

    def order_distance_weighted(self, mask_coords, mask):
        """Order pixels using weighted combination of multiple distance metrics."""
        height, width = mask.shape
        
        # Combine multiple ordering methods with weights
        nearest_order = self.order_nearest_first(mask_coords, mask)
        flow_order = self.order_flow_based(mask_coords, mask)
        coord_order = self.order_coordinate_based(mask_coords, mask)
        
        # Create weighted scores
        scores = np.zeros(len(mask_coords))
        
        for i, (_, y, x) in enumerate(nearest_order):
            scores[i] += 0.5 * (len(mask_coords) - i)  # Distance weight
        
        for i, (_, y, x) in enumerate(flow_order):
            scores[i] += 0.3 * (len(mask_coords) - i)  # Flow weight
        
        for i, (_, y, x) in enumerate(coord_order):
            scores[i] += 0.2 * (len(mask_coords) - i)  # Coordinate weight
        
        # Sort by combined score
        sorted_indices = np.argsort(-scores)  # Descending order
        
        ordered_coords = [(i, mask_coords[sorted_indices[i]][0], mask_coords[sorted_indices[i]][1]) 
                         for i in range(len(mask_coords))]
        
        return ordered_coords

    def divine_single_pixel_cascading(self, pixel_y, pixel_x, current_image, current_mask):
        """Divine a single pixel using the current state of the image (including previously divined pixels)."""
        
        # Use the appropriate divination method based on current settings
        method = self.divination_method.get()
        
        if method == "gs_enhanced_divine" and self.gs_alignment_data is not None:
            return self.divine_pixel_gs_enhanced_cascading(pixel_y, pixel_x, current_image, current_mask)
        elif method == "s_coordinate_flow":
            return self.divine_pixel_flow_cascading(pixel_y, pixel_x, current_image, current_mask)
        elif method == "s_coordinate_similarity":
            return self.divine_pixel_similarity_cascading(pixel_y, pixel_x, current_image, current_mask)
        elif method == "geometric_extrapolation":
            return self.divine_pixel_geometric_cascading(pixel_y, pixel_x, current_image, current_mask)
        else:
            # Default to combined approach
            return self.divine_pixel_combined_cascading(pixel_y, pixel_x, current_image, current_mask)

    def divine_pixel_gs_enhanced_cascading(self, pixel_y, pixel_x, current_image, current_mask):
        """Divine single pixel using G-S enhanced method with cascading."""
        if self.gs_coordinate_map is None:
            return self.divine_pixel_combined_cascading(pixel_y, pixel_x, current_image, current_mask)
        
        height, width = current_image.shape[:2]
        
        if 'g_theta' in self.gs_coordinate_map and 'g_phi' in self.gs_coordinate_map:
            s_x = self.gs_coordinate_map['s_x']
            s_y = self.gs_coordinate_map['s_y']
            g_theta = self.gs_coordinate_map['g_theta']
            g_phi = self.gs_coordinate_map['g_phi']
            
            # Resize if needed
            if s_x.shape != (height, width):
                s_x = cv2.resize(s_x, (width, height), interpolation=cv2.INTER_LINEAR)
                s_y = cv2.resize(s_y, (width, height), interpolation=cv2.INTER_LINEAR)
                g_theta = cv2.resize(g_theta, (width, height), interpolation=cv2.INTER_LINEAR)
                g_phi = cv2.resize(g_phi, (width, height), interpolation=cv2.INTER_LINEAR)
            
            return self.gs_coordinate_pixel_synthesis(pixel_y, pixel_x, s_x, s_y, g_theta, g_phi, 
                                                     current_image, current_mask)
        
        return self.divine_pixel_combined_cascading(pixel_y, pixel_x, current_image, current_mask)

    def divine_pixel_flow_cascading(self, pixel_y, pixel_x, current_image, current_mask):
        """Divine single pixel using flow analysis with cascading."""
        if self.gs_coordinate_map is None:
            return self.divine_pixel_geometric_cascading(pixel_y, pixel_x, current_image, current_mask)
        
        height, width = current_image.shape[:2]
        
        # Get coordinate maps
        if 's_x' in self.gs_coordinate_map and 's_y' in self.gs_coordinate_map:
            s_x = self.gs_coordinate_map['s_x']
            s_y = self.gs_coordinate_map['s_y']
        else:
            y_grid, x_grid = np.mgrid[0:height, 0:width]
            s_x = (x_grid / width) * 2.0 - 1.0
            s_y = (y_grid / height) * 2.0 - 1.0
        
        # Resize if needed
        if s_x.shape != (height, width):
            s_x = cv2.resize(s_x, (width, height), interpolation=cv2.INTER_LINEAR)
            s_y = cv2.resize(s_y, (width, height), interpolation=cv2.INTER_LINEAR)
        
        # Calculate flow field
        s_x_grad_x = np.gradient(s_x, axis=1)
        s_x_grad_y = np.gradient(s_x, axis=0)
        s_y_grad_x = np.gradient(s_y, axis=1)
        s_y_grad_y = np.gradient(s_y, axis=0)
        
        flow_x = s_x_grad_x + s_y_grad_x
        flow_y = s_x_grad_y + s_y_grad_y
        
        # Normalize flow vectors
        flow_magnitude = np.sqrt(flow_x**2 + flow_y**2)
        flow_magnitude[flow_magnitude == 0] = 1e-6
        
        flow_x_norm = flow_x / flow_magnitude
        flow_y_norm = flow_y / flow_magnitude
        
        # Trace flow to find source pixel
        return self.trace_s_coordinate_flow_source(pixel_x, pixel_y, flow_x_norm, flow_y_norm, 
                                                  s_x, s_y, current_image, current_mask)

    def divine_pixel_similarity_cascading(self, pixel_y, pixel_x, current_image, current_mask):
        """Divine single pixel using similarity matching with cascading."""
        height, width = current_image.shape[:2]
        
        # Get coordinate maps
        if self.gs_coordinate_map is not None and 's_x' in self.gs_coordinate_map:
            s_x = self.gs_coordinate_map['s_x']
            s_y = self.gs_coordinate_map['s_y']
        else:
            y_grid, x_grid = np.mgrid[0:height, 0:width]
            s_x = (x_grid / width) * 2.0 - 1.0
            s_y = (y_grid / height) * 2.0 - 1.0
        
        # Resize if needed
        if s_x.shape != (height, width):
            s_x = cv2.resize(s_x, (width, height), interpolation=cv2.INTER_LINEAR)
            s_y = cv2.resize(s_y, (width, height), interpolation=cv2.INTER_LINEAR)
        
        # Find valid pixels (not currently masked)
        valid_mask = ~current_mask
        valid_coords = np.column_stack(np.where(valid_mask))
        
        if len(valid_coords) == 0:
            return None
        
        # S-coordinates of valid pixels
        valid_s_coords = np.column_stack([s_x[valid_mask], s_y[valid_mask]])
        
        # Build KD-tree
        try:
            from scipy.spatial import cKDTree
            s_coord_tree = cKDTree(valid_s_coords)
            
            # Get S-coordinates at pixel location
            pixel_s_x = s_x[pixel_y, pixel_x]
            pixel_s_y = s_y[pixel_y, pixel_x]
            pixel_s_coord = np.array([pixel_s_x, pixel_s_y])
            
            # Find closest matches
            num_neighbors = min(5, len(valid_s_coords))
            distances, indices = s_coord_tree.query(pixel_s_coord, k=num_neighbors)
            
            # Check if closest match is within threshold
            if distances[0] <= self.divine_similarity_threshold.get():
                # Weight by similarity
                weights = 1.0 / (distances + 1e-6)
                weights = weights / np.sum(weights)
                
                # Weighted average
                divined_pixel = np.zeros(3)
                for i, idx in enumerate(indices):
                    pixel_y_valid, pixel_x_valid = valid_coords[idx]
                    pixel_value = current_image[pixel_y_valid, pixel_x_valid]
                    divined_pixel += pixel_value * weights[i]
                
                return divined_pixel.astype(np.uint8)
                
        except Exception as e:
            if self.debug_mode.get():
                print(f"   Similarity matching failed: {e}")
        
        return None

    def divine_pixel_geometric_cascading(self, pixel_y, pixel_x, current_image, current_mask):
        """Divine single pixel using geometric extrapolation with cascading."""
        height, width = current_image.shape[:2]
        
        # Find valid pixels (not currently masked)
        valid_mask = ~current_mask
        valid_coords = np.column_stack(np.where(valid_mask))
        
        if len(valid_coords) < 3:
            return None
        
        valid_y, valid_x = valid_coords[:, 0], valid_coords[:, 1]
        
        # Find closest valid pixels
        distances = np.sqrt((valid_y - pixel_y)**2 + (valid_x - pixel_x)**2)
        
        # Use closest pixels for extrapolation
        search_radius = self.divine_search_radius.get()
        close_mask = distances <= search_radius
        
        if not np.any(close_mask):
            # If no pixels within search radius, use closest few
            num_closest = min(8, len(distances))
            closest_indices = np.argpartition(distances, num_closest)[:num_closest]
            close_mask = np.zeros_like(distances, dtype=bool)
            close_mask[closest_indices] = True
        
        close_coords = valid_coords[close_mask]
        close_distances = distances[close_mask]
        
        # Weight by distance
        weights = 1.0 / (close_distances + 1e-6)
        weights = weights / np.sum(weights)
        
        # Extrapolate pixel value
        close_pixels = current_image[close_coords[:, 0], close_coords[:, 1]]
        divined_pixel = np.average(close_pixels, weights=weights, axis=0)
        
        return divined_pixel.astype(np.uint8)

    def divine_pixel_combined_cascading(self, pixel_y, pixel_x, current_image, current_mask):
        """Divine single pixel using combined methods with cascading."""
        
        # Try multiple methods in order of preference
        methods = [
            self.divine_pixel_similarity_cascading,
            self.divine_pixel_flow_cascading,
            self.divine_pixel_geometric_cascading
        ]
        
        for method in methods:
            try:
                result = method(pixel_y, pixel_x, current_image, current_mask)
                if result is not None:
                    return result
            except Exception as e:
                if self.debug_mode.get():
                    print(f"   Method {method.__name__} failed: {e}")
                continue
        
        return None

    # === ALL OTHER ORIGINAL METHODS PRESERVED EXACTLY ===
    
    def load_image(self):
        """Load image to divine."""
        file_path = filedialog.askopenfilename(
            title="Select Image to Divine",
            filetypes=[("Image files", "*.png *.jpg *.jpeg *.bmp *.tiff"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                pil_image = Image.open(file_path)
                self.original_image = np.array(pil_image.convert('RGB'))
                
                print(f"🔮 Loaded image for divination: {self.original_image.shape}")
                self.update_status()
                self.show_loaded_image()
                
                # Auto-generate coordinate map if enabled and no G-S data
                if self.auto_coordinate_generation.get() and self.gs_alignment_data is None:
                    self.generate_coordinate_map()
                    
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load image: {str(e)}")

    def load_gs_alignment_data(self):
        """Load G-S alignment data file (same format as G-S Stereo Viewer)."""
        file_path = filedialog.askopenfilename(
            title="Select G-S Alignment Data",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                with open(file_path, 'r') as f:
                    self.gs_alignment_data = json.load(f)
                
                print(f"📊 Loaded G-S alignment data: {len(self.gs_alignment_data)} alignments")
                
                # Validate G-S alignment data format
                if self.gs_alignment_data:
                    sample = self.gs_alignment_data[0]
                    print(f"   Sample alignment keys: {list(sample.keys())}")
                    
                    # Check for required G-S fields
                    required_fields = ['S_x', 'S_y', 'G_theta', 'G_phi']
                    missing_fields = [field for field in required_fields if field not in sample]
                    
                    if missing_fields:
                        messagebox.showwarning("Incomplete G-S Data", 
                                             f"Missing required G-S fields: {missing_fields}\n" +
                                             "Divination may be limited.")
                    else:
                        print("   ✅ All required G-S alignment fields present")
                        
                        # Check for optional Hadit fields
                        if 'Hadit_theta' in sample and 'Hadit_phi' in sample:
                            print("   🌟 Hadit vector data available")
                        else:
                            print("   ⚠️ No Hadit vector data (optional)")
                
                self.update_status()
                
                # Enable coordinate generation button
                if self.original_image is not None:
                    self.generate_gs_coords_button.config(state="normal")
                    
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load G-S alignment data: {str(e)}")

    def load_mask(self):
        """Load mask image (optional)."""
        file_path = filedialog.askopenfilename(
            title="Select Mask Image",
            filetypes=[("Image files", "*.png *.jpg *.jpeg *.bmp *.tiff"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                pil_image = Image.open(file_path)
                if pil_image.mode == 'RGB':
                    gray_image = pil_image.convert('L')
                    self.mask_image = np.array(gray_image) > 128  # Convert to boolean
                else:
                    self.mask_image = np.array(pil_image) > 128
                
                print(f"🎭 Loaded mask: {self.mask_image.shape}")
                self.update_status()
                    
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load mask: {str(e)}")
                
    def load_coordinate_map(self):
        """Load coordinate map (optional)."""
        file_path = filedialog.askopenfilename(
            title="Select Coordinate Map",
            filetypes=[("NumPy files", "*.npz"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                self.gs_coordinate_map = np.load(file_path)
                print(f"🗺️ Loaded coordinate map with keys: {list(self.gs_coordinate_map.keys())}")
                self.update_status()
                    
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load coordinate map: {str(e)}")

    def generate_gs_coordinate_maps(self):
        """Generate G-S coordinate maps from alignment data (G-S Stereo Viewer compatible)."""
        if self.gs_alignment_data is None:
            messagebox.showwarning("No Data", "Load G-S alignment data first")
            return
            
        if self.original_image is None:
            messagebox.showwarning("No Image", "Load an image first")
            return
            
        try:
            print("🧮 Generating G-S coordinate maps from alignment data...")
            
            height, width = self.original_image.shape[:2]
            
            # Extract coordinate data from G-S alignments (exact same method as G-S viewer)
            s_x_points = [a['S_x'] for a in self.gs_alignment_data]
            s_y_points = [a['S_y'] for a in self.gs_alignment_data]
            g_theta_points = [a['G_theta'] for a in self.gs_alignment_data]
            g_phi_points = [a['G_phi'] for a in self.gs_alignment_data]
            
            # Create grid coordinates
            y_grid, x_grid = np.mgrid[0:height, 0:width]
            
            # Normalize grid coordinates to match S-coordinate space (G-S method)
            s_x_min, s_x_max = min(s_x_points), max(s_x_points)
            s_y_min, s_y_max = min(s_y_points), max(s_y_points)
            
            if s_x_max > s_x_min and s_y_max > s_y_min:
                x_normalized = (x_grid / width) * (s_x_max - s_x_min) + s_x_min
                y_normalized = (y_grid / height) * (s_y_max - s_y_min) + s_y_min
            else:
                x_normalized = (x_grid / width) * 2.0 - 1.0
                y_normalized = (y_grid / height) * 2.0 - 1.0
            
            # Create point arrays for interpolation (G-S method)
            alignment_points = np.column_stack([s_x_points, s_y_points])
            grid_points = np.column_stack([x_normalized.ravel(), y_normalized.ravel()])
            
            # Interpolate G-S coordinates (same as G-S viewer)
            try:
                g_theta_interp = griddata(alignment_points, g_theta_points, grid_points, 
                                        method='linear', fill_value=0).reshape((height, width))
                g_phi_interp = griddata(alignment_points, g_phi_points, grid_points, 
                                      method='linear', fill_value=0).reshape((height, width))
            except:
                # Fallback to nearest neighbor if linear fails
                g_theta_interp = griddata(alignment_points, g_theta_points, grid_points, 
                                        method='nearest', fill_value=0).reshape((height, width))
                g_phi_interp = griddata(alignment_points, g_phi_points, grid_points, 
                                      method='nearest', fill_value=0).reshape((height, width))
            
            # Store coordinate maps (G-S compatible format)
            self.gs_coordinate_map = {
                'g_theta': g_theta_interp,
                'g_phi': g_phi_interp,
                's_x': x_normalized,
                's_y': y_normalized,
                'x_grid': x_grid,
                'y_grid': y_grid,
                'alignment_source': True
            }
            
            # Add Hadit data if available (G-S enhancement)
            if self.hadit_enhancement.get() and len(self.gs_alignment_data) > 0:
                sample = self.gs_alignment_data[0]
                if 'Hadit_theta' in sample and 'Hadit_phi' in sample:
                    hadit_theta_points = [a['Hadit_theta'] for a in self.gs_alignment_data]
                    hadit_phi_points = [a['Hadit_phi'] for a in self.gs_alignment_data]
                    
                    try:
                        hadit_theta_interp = griddata(alignment_points, hadit_theta_points, grid_points, 
                                                    method='linear', fill_value=0).reshape((height, width))
                        hadit_phi_interp = griddata(alignment_points, hadit_phi_points, grid_points, 
                                                  method='linear', fill_value=0).reshape((height, width))
                        
                        self.gs_coordinate_map['hadit_theta'] = hadit_theta_interp
                        self.gs_coordinate_map['hadit_phi'] = hadit_phi_interp
                        
                        print("   🌟 Hadit vector enhancement included")
                    except Exception as e:
                        print(f"   ⚠️ Failed to interpolate Hadit data: {e}")
            
            print(f"✅ Generated G-S coordinate maps from {len(self.gs_alignment_data)} alignments")
            print(f"   G_theta range: [{g_theta_interp.min():.1f}°, {g_theta_interp.max():.1f}°]")
            print(f"   G_phi range: [{g_phi_interp.min():.1f}°, {g_phi_interp.max():.1f}°]")
            print(f"   S_x range: [{x_normalized.min():.3f}, {x_normalized.max():.3f}]")
            print(f"   S_y range: [{y_normalized.min():.3f}, {y_normalized.max():.3f}]")
            
            # Visualize coordinate maps
            self.visualize_gs_coordinate_maps()
            
            messagebox.showinfo("Success", f"🧮 Generated G-S coordinate maps from {len(self.gs_alignment_data)} alignments")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate G-S coordinates: {str(e)}")
            print(f"Error details: {str(e)}")

    def generate_coordinate_map(self):
        """Generate basic S-coordinate map for divination (fallback when no G-S data)."""
        if self.original_image is None:
            return
            
        height, width = self.original_image.shape[:2]
        
        # Create coordinate grids
        y_grid, x_grid = np.mgrid[0:height, 0:width]
        
        # Normalize coordinates to [-1, 1] range (S-coordinate space)
        s_x = (x_grid / width) * 2.0 - 1.0
        s_y = (y_grid / height) * 2.0 - 1.0
        
        # Create coordinate map
        self.gs_coordinate_map = {
            's_x': s_x,
            's_y': s_y,
            'x_grid': x_grid,
            'y_grid': y_grid,
            'alignment_source': False
        }
        
        if self.debug_mode.get():
            print(f"🗺️ Generated basic coordinate map: {s_x.shape}")
            print(f"   S_x range: [{s_x.min():.3f}, {s_x.max():.3f}]")
            print(f"   S_y range: [{s_y.min():.3f}, {s_y.max():.3f}]")

    # === MAIN DIVINATION FUNCTION (PRESERVED WITH THREADING UPDATES) ===
    
    def divine_pixels(self):
        """Main divine pixels function with G-S enhancement and cascading support (PRESERVED)."""
        if self.original_image is None:
            messagebox.showwarning("No Data", "Load an image first")
            return
            
        try:
            start_time = time.time()
            
            print("🔮 Starting G-S Compatible Divine Pixel Process with Cascading...")
            
            # Step 1: Generate coordinate map if not available
            if self.gs_coordinate_map is None:
                if self.gs_alignment_data is not None and self.use_gs_alignment_data.get():
                    self.generate_gs_coordinate_maps()
                elif self.auto_coordinate_generation.get():
                    self.generate_coordinate_map()
            
            self.update_progress_threaded(20)
            
            # Step 2: Create expanded canvas if divine outside bounds is enabled
            if self.divine_outside_bounds.get():
                self.divined_image = self.create_expanded_canvas()
            else:
                self.divined_image = self.original_image.copy()
            
            self.update_progress_threaded(40)
            
            # Step 3: Create divination mask
            divination_mask = self.create_divination_mask()
            
            self.update_progress_threaded(60)
            
            # Step 4: Apply chosen divination method
            method = self.divination_method.get()
            
            # Use cascading if enabled or if method is specifically cascading
            if self.enable_cascading.get() or method == "cascading_divine":
                print("🌊 Using cascading divination process...")
                self.divined_image = self.cascading_divine_process(self.divined_image, divination_mask)
            else:
                # Use traditional non-cascading methods
                if method == "s_coordinate_flow":
                    self.divined_image = self.s_coordinate_flow_divination(
                        self.divined_image, divination_mask)
                elif method == "s_coordinate_similarity":
                    self.divined_image = self.s_coordinate_similarity_divination(
                        self.divined_image, divination_mask)
                elif method == "geometric_extrapolation":
                    self.divined_image = self.geometric_extrapolation_divination(
                        self.divined_image, divination_mask)
                elif method == "combined_divine":
                    self.divined_image = self.combined_divine_process(
                        self.divined_image, divination_mask)
                elif method == "gs_enhanced_divine":
                    self.divined_image = self.gs_enhanced_divine_process(
                        self.divined_image, divination_mask)
                elif method == "omniscient_divine":
                    self.divined_image = self.omniscient_divine_process(
                        self.divined_image, divination_mask)
            
            self.update_progress_threaded(80)
            
            # Step 5: Post-processing
            if self.noise_reduction.get():
                self.divined_image = self.apply_noise_reduction(self.divined_image)
            
            processing_time = time.time() - start_time
            
            # Update status with cascading info (THREAD-SAFE)
            cascade_info = ""
            if self.enable_cascading.get() or method == "cascading_divine":
                cascade_info = f" (Cascading: {self.cascading_stats['successful_pixels']}/{self.cascading_stats['total_pixels']} in {self.cascading_stats['cascade_generations']} gen)"
            
            status_text = f"✅ G-S Divination completed ({processing_time:.2f}s){cascade_info}"
            self.update_status_threaded(status_text, "green")
            
            self.export_button.config(state="normal")
            self.export_maps_button.config(state="normal")
            
            # Update statistics
            self.update_statistics()
            
            print(f"🔮 G-S compatible divine pixel process completed in {processing_time:.2f} seconds")
            if cascade_info:
                print(f"🌊 Cascading statistics: {self.cascading_stats}")
            self.visualize_results()
            
        except Exception as e:
            error_msg = f"Failed to divine pixels: {str(e)}"
            print(f"Error details: {str(e)}")
            traceback.print_exc()
            raise Exception(error_msg)
        finally:
            self.update_progress_threaded(100)

    # === ALL EXISTING DIVINATION METHODS (PRESERVED EXACTLY) ===
    
    def gs_enhanced_divine_process(self, image, mask):
        """G-S enhanced divine process using alignment data."""
        if self.gs_alignment_data is None or not self.use_gs_alignment_data.get():
            if self.debug_mode.get():
                print("⚠️ No G-S alignment data available, falling back to combined divine")
            return self.combined_divine_process(image, mask)
        
        if self.debug_mode.get():
            print("📊 Performing G-S enhanced divine process...")
        
        result = image.copy()
        height, width = result.shape[:2]
        
        # Use coordinate maps generated from G-S alignment data
        if self.gs_coordinate_map is not None and self.gs_coordinate_map.get('alignment_source', False):
            s_x = self.gs_coordinate_map['s_x']
            s_y = self.gs_coordinate_map['s_y']
            g_theta = self.gs_coordinate_map['g_theta']
            g_phi = self.gs_coordinate_map['g_phi']
            
            # Process mask pixels using G-S enhanced methods
            mask_coords = np.column_stack(np.where(mask))
            divine_success = 0
            
            for mask_y, mask_x in mask_coords:
                # Check for cancellation
                if self.cancel_event.is_set():
                    break
                    
                # Use G-S alignment data for precise pixel synthesis
                divined_pixel = self.gs_coordinate_pixel_synthesis(
                    mask_y, mask_x, s_x, s_y, g_theta, g_phi, image, mask)
                
                if divined_pixel is not None:
                    result[mask_y, mask_x] = divined_pixel
                    divine_success += 1
            
            if self.debug_mode.get():
                print(f"   G-S enhanced divine success: {divine_success}/{len(mask_coords)} pixels")
        else:
            if self.debug_mode.get():
                print("   No G-S-generated coordinates available, using combined method")
            result = self.combined_divine_process(result, mask)
        
        return result

    def gs_coordinate_pixel_synthesis(self, mask_y, mask_x, s_x, s_y, g_theta, g_phi, source_image, mask):
        """Synthesize pixel using G-S coordinate data."""
        height, width = s_x.shape
        
        # Get coordinates at mask location
        mask_s_x = s_x[mask_y, mask_x]
        mask_s_y = s_y[mask_y, mask_x]
        mask_g_theta = g_theta[mask_y, mask_x]
        mask_g_phi = g_phi[mask_y, mask_x]
        
        # Find similar coordinates in G-S alignment data
        best_match = self.find_best_gs_alignment_match(mask_s_x, mask_s_y, mask_g_theta, mask_g_phi)
        
        if best_match is not None:
            # Use G-S alignment data to guide pixel synthesis
            return self.synthesize_from_gs_alignment_match(
                mask_y, mask_x, best_match, s_x, s_y, source_image, mask)
        
        # Fallback to coordinate-based synthesis
        return self.coordinate_based_synthesis(mask_y, mask_x, s_x, s_y, source_image, mask)

    def find_best_gs_alignment_match(self, target_s_x, target_s_y, target_g_theta, target_g_phi):
        """Find best matching G-S alignment from data."""
        if self.gs_alignment_data is None:
            return None
        
        best_distance = float('inf')
        best_match = None
        
        precision = self.gs_coordinate_precision.get()
        
        for alignment in self.gs_alignment_data:
            # Calculate coordinate distances
            s_distance = np.sqrt((alignment['S_x'] - target_s_x)**2 + 
                               (alignment['S_y'] - target_s_y)**2)
            
            # Handle theta wraparound (G-S method)
            theta_diff = abs(alignment['G_theta'] - target_g_theta)
            theta_diff = min(theta_diff, 360 - theta_diff)
            
            phi_diff = abs(alignment['G_phi'] - target_g_phi)
            
            g_distance = np.sqrt(theta_diff**2 + phi_diff**2)
            
            # Combined distance (weight S-coordinates more heavily)
            combined_distance = s_distance * 0.7 + g_distance * 0.3 / 180.0  # Normalize G distance
            
            if combined_distance < best_distance and combined_distance <= precision:
                best_distance = combined_distance
                best_match = alignment
        
        return best_match

    def synthesize_from_gs_alignment_match(self, mask_y, mask_x, alignment, s_x, s_y, source_image, mask):
        """Synthesize pixel based on G-S alignment match."""
        height, width = s_x.shape
        
        # Use G-S alignment data to find similar coordinate locations
        target_s_x = alignment['S_x']
        target_s_y = alignment['S_y']
        
        # Find pixels with similar S-coordinates
        s_distances = np.sqrt((s_x - target_s_x)**2 + (s_y - target_s_y)**2)
        
        # Find valid pixels (not masked) with similar coordinates
        valid_mask = ~mask
        valid_similar = valid_mask & (s_distances < self.divine_precision.get())
        
        if np.any(valid_similar):
            # Weight by coordinate similarity
            similar_distances = s_distances[valid_similar]
            weights = 1.0 / (similar_distances + 1e-6)
            weights = weights / np.sum(weights)
            
            # Weighted average of similar pixels
            similar_pixels = source_image[valid_similar]
            synthesized_pixel = np.average(similar_pixels, weights=weights, axis=0)
            
            # Apply Hadit enhancement if available
            if (self.hadit_enhancement.get() and 
                'hadit_theta' in self.gs_coordinate_map and 
                'Hadit_theta' in alignment):
                synthesized_pixel = self.apply_hadit_enhancement(
                    synthesized_pixel, alignment, mask_y, mask_x)
            
            return synthesized_pixel.astype(np.uint8)
        
        return None

    def apply_hadit_enhancement(self, pixel, alignment, mask_y, mask_x):
        """Apply Hadit vector enhancement to synthesized pixel."""
        try:
            # Get Hadit vectors from alignment and coordinate map
            alignment_hadit_theta = alignment.get('Hadit_theta', 0) * np.pi / 180.0
            alignment_hadit_phi = alignment.get('Hadit_phi', 0) * np.pi / 180.0
            
            if 'hadit_theta' in self.gs_coordinate_map:
                local_hadit_theta = self.gs_coordinate_map['hadit_theta'][mask_y, mask_x] * np.pi / 180.0
                local_hadit_phi = self.gs_coordinate_map['hadit_phi'][mask_y, mask_x] * np.pi / 180.0
                
                # Calculate Hadit vectors
                alignment_hadit = np.array([
                    np.sin(alignment_hadit_phi) * np.cos(alignment_hadit_theta),
                    np.sin(alignment_hadit_phi) * np.sin(alignment_hadit_theta),
                    np.cos(alignment_hadit_phi)
                ])
                
                local_hadit = np.array([
                    np.sin(local_hadit_phi) * np.cos(local_hadit_theta),
                    np.sin(local_hadit_phi) * np.sin(local_hadit_theta),
                    np.cos(local_hadit_phi)
                ])
                
                # Calculate enhancement factor based on Hadit similarity
                hadit_similarity = np.dot(alignment_hadit, local_hadit)
                enhancement_factor = 1.0 + (hadit_similarity * 0.1)  # 10% max enhancement
                
                # Apply enhancement to pixel
                enhanced_pixel = pixel * enhancement_factor
                enhanced_pixel = np.clip(enhanced_pixel, 0, 255)
                
                return enhanced_pixel
        except Exception as e:
            if self.debug_mode.get():
                print(f"   Hadit enhancement failed: {e}")
        
        return pixel

    def create_expanded_canvas(self):
        """Create expanded canvas for divining outside image bounds."""
        if self.original_image is None:
            return None
            
        height, width = self.original_image.shape[:2]
        expansion = self.edge_extension.get()
        
        # Create larger canvas
        new_height = height + 2 * expansion
        new_width = width + 2 * expansion
        
        expanded_image = np.zeros((new_height, new_width, 3), dtype=np.uint8)
        
        # Place original image in center
        expanded_image[expansion:expansion+height, expansion:expansion+width] = self.original_image
        
        if self.debug_mode.get():
            print(f"🚀 Created expanded canvas: {expanded_image.shape}")
        
        return expanded_image

    def create_divination_mask(self):
        """Create mask indicating which pixels need divination."""
        if self.divined_image is None:
            return None
            
        height, width = self.divined_image.shape[:2]
        
        if self.mask_image is not None:
            # Use provided mask, resize if needed
            if self.mask_image.shape != (height, width):
                mask = cv2.resize(self.mask_image.astype(np.uint8), (width, height), 
                                interpolation=cv2.INTER_NEAREST).astype(bool)
            else:
                mask = self.mask_image.copy()
        else:
            # Create mask for pixels that need divination
            if self.divine_outside_bounds.get():
                # Mask everything outside original image area
                orig_height, orig_width = self.original_image.shape[:2]
                expansion = self.edge_extension.get()
                
                mask = np.ones((height, width), dtype=bool)
                mask[expansion:expansion+orig_height, expansion:expansion+orig_width] = False
            else:
                # Mask black pixels or use entire image
                mask = np.all(self.divined_image == 0, axis=2)
        
        if self.debug_mode.get():
            print(f"🎭 Created divination mask: {np.sum(mask)} pixels to divine")
        
        return mask

    def s_coordinate_flow_divination(self, image, mask):
        """Divine pixels using S-coordinate flow analysis (G-S method)."""
        if self.gs_coordinate_map is None:
            if self.debug_mode.get():
                print("⚠️ No G-S coordinate map available for flow analysis")
            return image
            
        if self.debug_mode.get():
            print("🌊 Performing S-coordinate flow divination...")
        
        result = image.copy()
        height, width = result.shape[:2]
        
        # Get coordinate maps
        if 's_x' in self.gs_coordinate_map and 's_y' in self.gs_coordinate_map:
            s_x = self.gs_coordinate_map['s_x']
            s_y = self.gs_coordinate_map['s_y']
        else:
            # Generate basic coordinates
            y_grid, x_grid = np.mgrid[0:height, 0:width]
            s_x = (x_grid / width) * 2.0 - 1.0
            s_y = (y_grid / height) * 2.0 - 1.0
        
        # Resize if needed
        if s_x.shape != (height, width):
            s_x = cv2.resize(s_x, (width, height), interpolation=cv2.INTER_LINEAR)
            s_y = cv2.resize(s_y, (width, height), interpolation=cv2.INTER_LINEAR)
        
        # Calculate S-coordinate flow field (gradients) - G-S method
        s_x_grad_x = np.gradient(s_x, axis=1)
        s_x_grad_y = np.gradient(s_x, axis=0)
        s_y_grad_x = np.gradient(s_y, axis=1)
        s_y_grad_y = np.gradient(s_y, axis=0)
        
        # Flow field vectors
        flow_x = s_x_grad_x + s_y_grad_x
        flow_y = s_x_grad_y + s_y_grad_y
        
        # Normalize flow vectors
        flow_magnitude = np.sqrt(flow_x**2 + flow_y**2)
        flow_magnitude[flow_magnitude == 0] = 1e-6
        
        flow_x_norm = flow_x / flow_magnitude
        flow_y_norm = flow_y / flow_magnitude
        
        # Process mask pixels
        mask_coords = np.column_stack(np.where(mask))
        divine_success = 0
        
        for mask_y, mask_x in mask_coords:
            # Check for cancellation
            if self.cancel_event.is_set():
                break
                
            # Trace flow to find source pixel
            source_pixel = self.trace_s_coordinate_flow_source(
                mask_x, mask_y, flow_x_norm, flow_y_norm, 
                s_x, s_y, image, mask)
            
            if source_pixel is not None:
                result[mask_y, mask_x] = source_pixel
                divine_success += 1
        
        if self.debug_mode.get():
            print(f"   Flow divination success: {divine_success}/{len(mask_coords)} pixels")
        
        return result

    def trace_s_coordinate_flow_source(self, hole_x, hole_y, flow_dir_x, flow_dir_y, 
                                      s_x, s_y, original_texture, holes_mask):
        """Trace along S-coordinate flow to find source pixel (G-S method)."""
        
        height, width = s_x.shape
        
        # Current S-coordinate values at hole
        if hole_y >= height or hole_x >= width:
            return None
            
        target_s_x = s_x[hole_y, hole_x]
        target_s_y = s_y[hole_y, hole_x]
        
        # Trace backwards along flow line
        max_trace_distance = self.divine_search_radius.get()
        step_size = 0.5
        
        current_x = float(hole_x)
        current_y = float(hole_y)
        
        for step in range(int(max_trace_distance / step_size)):
            # Move backwards along flow
            if hole_y < height and hole_x < width:
                current_x -= flow_dir_x[hole_y, hole_x] * step_size
                current_y -= flow_dir_y[hole_y, hole_x] * step_size
            
            # Check bounds
            if current_x < 0 or current_x >= width-1 or current_y < 0 or current_y >= height-1:
                break
            
            # Get sample coordinates
            sample_x = int(np.round(current_x))
            sample_y = int(np.round(current_y))
            
            # Check if this location has valid data
            if not holes_mask[sample_y, sample_x]:
                # Calculate S-coordinate similarity
                source_s_x = s_x[sample_y, sample_x]
                source_s_y = s_y[sample_y, sample_x]
                
                s_distance = np.sqrt((source_s_x - target_s_x)**2 + (source_s_y - target_s_y)**2)
                
                # If S-coordinates are similar enough, use this pixel
                if s_distance <= self.divine_precision.get():
                    return original_texture[sample_y, sample_x]
        
        return None

    def s_coordinate_similarity_divination(self, image, mask):
        """Divine pixels using S-coordinate similarity matching (G-S method)."""
        if self.debug_mode.get():
            print("🔍 Performing G-S coordinate similarity divination...")
        
        result = image.copy()
        height, width = result.shape[:2]
        
        # Get coordinate maps
        if self.gs_coordinate_map is not None and 's_x' in self.gs_coordinate_map:
            s_x = self.gs_coordinate_map['s_x']
            s_y = self.gs_coordinate_map['s_y']
        else:
            y_grid, x_grid = np.mgrid[0:height, 0:width]
            s_x = (x_grid / width) * 2.0 - 1.0
            s_y = (y_grid / height) * 2.0 - 1.0
        
        # Resize if needed
        if s_x.shape != (height, width):
            s_x = cv2.resize(s_x, (width, height), interpolation=cv2.INTER_LINEAR)
            s_y = cv2.resize(s_y, (width, height), interpolation=cv2.INTER_LINEAR)
        
        # Build KD-tree for fast lookups
        valid_mask = ~mask
        valid_coords = np.column_stack(np.where(valid_mask))
        
        if len(valid_coords) == 0:
            if self.debug_mode.get():
                print("   Warning: No valid pixels for similarity matching")
            return result
        
        # S-coordinates of valid pixels
        valid_s_coords = np.column_stack([
            s_x[valid_mask].ravel(),
            s_y[valid_mask].ravel()
        ])
        
        # Build KD-tree
        s_coord_tree = cKDTree(valid_s_coords)
        
        # Process mask pixels
        mask_coords = np.column_stack(np.where(mask))
        divine_success = 0
        
        for mask_y, mask_x in mask_coords:
            # Check for cancellation
            if self.cancel_event.is_set():
                break
                
            # Get S-coordinates at mask location
            mask_s_x = s_x[mask_y, mask_x]
            mask_s_y = s_y[mask_y, mask_x]
            mask_s_coord = np.array([mask_s_x, mask_s_y])
            
            # Find closest matches
            num_neighbors = min(5, len(valid_s_coords))
            try:
                distances, indices = s_coord_tree.query(mask_s_coord, k=num_neighbors)
                
                # Check if closest match is within threshold
                if distances[0] <= self.divine_similarity_threshold.get():
                    # Weight by similarity
                    weights = 1.0 / (distances + 1e-6)
                    weights = weights / np.sum(weights)
                    
                    # Weighted average
                    divined_pixel = np.zeros(3)
                    for i, idx in enumerate(indices):
                        pixel_y, pixel_x = valid_coords[idx]
                        pixel_value = image[pixel_y, pixel_x]
                        divined_pixel += pixel_value * weights[i]
                    
                    result[mask_y, mask_x] = divined_pixel.astype(np.uint8)
                    divine_success += 1
                    
            except Exception as e:
                if self.debug_mode.get():
                    print(f"   KDTree query failed: {e}")
                continue
        
        if self.debug_mode.get():
            print(f"   G-S similarity divination success: {divine_success}/{len(mask_coords)} pixels")
        
        return result

    def geometric_extrapolation_divination(self, image, mask):
        """Divine pixels using geometric extrapolation."""
        if self.debug_mode.get():
            print("📐 Performing geometric extrapolation divination...")
        
        result = image.copy()
        height, width = result.shape[:2]
        
        # Get valid pixels
        valid_mask = ~mask
        
        if np.sum(valid_mask) < 4:
            if self.debug_mode.get():
                print("   Warning: Insufficient valid pixels for extrapolation")
            return result
        
        valid_y, valid_x = np.where(valid_mask)
        valid_points = np.column_stack([valid_y, valid_x])
        valid_texture = image[valid_mask]
        
        # Process mask pixels
        mask_coords = np.column_stack(np.where(mask))
        divine_success = 0
        
        for mask_y, mask_x in mask_coords:
            # Check for cancellation
            if self.cancel_event.is_set():
                break
                
            # Find closest valid pixels
            distances = np.sqrt((valid_y - mask_y)**2 + (valid_x - mask_x)**2)
            
            # Use closest pixels for extrapolation
            num_closest = min(8, len(distances))
            closest_indices = np.argpartition(distances, num_closest)[:num_closest]
            
            # Weight by distance
            weights = 1.0 / (distances[closest_indices] + 1e-6)
            weights = weights / np.sum(weights)
            
            # Extrapolate pixel value
            divined_pixel = np.average(valid_texture[closest_indices], 
                                     weights=weights, axis=0)
            
            result[mask_y, mask_x] = divined_pixel.astype(np.uint8)
            divine_success += 1
        
        if self.debug_mode.get():
            print(f"   Geometric extrapolation success: {divine_success}/{len(mask_coords)} pixels")
        
        return result

    def combined_divine_process(self, image, mask):
        """Multi-method divine process."""
        if self.debug_mode.get():
            print("🔮 Performing combined divine process...")
        
        result = image.copy()
        remaining_mask = mask.copy()
        
        methods = []
        if self.use_flow_analysis.get():
            methods.append(("flow", self.s_coordinate_flow_divination))
        if self.use_similarity_matching.get():
            methods.append(("similarity", self.s_coordinate_similarity_divination))
        if self.use_geometric_extrapolation.get():
            methods.append(("geometric", self.geometric_extrapolation_divination))
        
        for method_name, method_func in methods:
            if np.any(remaining_mask) and not self.cancel_event.is_set():
                result = method_func(result, remaining_mask)
                # Update remaining mask
                remaining_mask = np.all(result == 0, axis=2) & mask
                if self.debug_mode.get():
                    print(f"   After {method_name}: {np.sum(remaining_mask)} pixels remaining")
        
        return result

    def omniscient_divine_process(self, image, mask):
        """Ultimate omniscient divine process using all methods including G-S."""
        if self.debug_mode.get():
            print("🌟 Performing omniscient divine process...")
        
        # First apply G-S enhanced if available
        if self.gs_alignment_data is not None and self.use_gs_alignment_data.get():
            result = self.gs_enhanced_divine_process(image, mask)
        else:
            result = self.combined_divine_process(image, mask)
        
        # Then apply gradient synthesis for remaining pixels
        remaining_mask = np.all(result == 0, axis=2) & mask
        if np.any(remaining_mask) and self.use_gradient_synthesis.get() and not self.cancel_event.is_set():
            result = self.gradient_synthesis_divination(result, remaining_mask)
        
        return result

    def gradient_synthesis_divination(self, image, mask):
        """Divine pixels using gradient synthesis."""
        if self.debug_mode.get():
            print("🧮 Performing gradient synthesis divination...")
        
        result = image.copy()
        height, width = result.shape[:2]
        
        # Calculate image gradients
        gray_image = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
        grad_x = cv2.Sobel(gray_image, cv2.CV_64F, 1, 0, ksize=3)
        grad_y = cv2.Sobel(gray_image, cv2.CV_64F, 0, 1, ksize=3)
        
        # Process mask pixels
        mask_coords = np.column_stack(np.where(mask))
        divine_success = 0
        
        for mask_y, mask_x in mask_coords:
            # Check for cancellation
            if self.cancel_event.is_set():
                break
                
            # Find gradient pattern in neighborhood
            search_radius = self.divine_search_radius.get() // 2
            y_min = max(0, mask_y - search_radius)
            y_max = min(height, mask_y + search_radius + 1)
            x_min = max(0, mask_x - search_radius)
            x_max = min(width, mask_x + search_radius + 1)
            
            # Extract neighborhood
            region_image = image[y_min:y_max, x_min:x_max]
            region_mask = mask[y_min:y_max, x_min:x_max]
            
            # Find valid pixels in region
            valid_region = ~region_mask
            if np.any(valid_region):
                # Use median of valid pixels as synthesis
                valid_pixels = region_image[valid_region]
                synthesized_pixel = np.median(valid_pixels, axis=0)
                
                result[mask_y, mask_x] = synthesized_pixel.astype(np.uint8)
                divine_success += 1
        
        if self.debug_mode.get():
            print(f"   Gradient synthesis success: {divine_success}/{len(mask_coords)} pixels")
        
        return result

    def coordinate_based_synthesis(self, mask_y, mask_x, s_x, s_y, source_image, mask):
        """Fallback coordinate-based synthesis."""
        # Find nearby valid pixels
        search_radius = self.divine_search_radius.get()
        y_min = max(0, mask_y - search_radius)
        y_max = min(s_x.shape[0], mask_y + search_radius + 1)
        x_min = max(0, mask_x - search_radius)
        x_max = min(s_x.shape[1], mask_x + search_radius + 1)
        
        # Extract neighborhood
        region_mask = mask[y_min:y_max, x_min:x_max]
        region_image = source_image[y_min:y_max, x_min:x_max]
        
        # Find valid pixels in region
        valid_region = ~region_mask
        if np.any(valid_region):
            valid_pixels = region_image[valid_region]
            return np.median(valid_pixels, axis=0).astype(np.uint8)
        
        return None

    def apply_noise_reduction(self, image):
        """Apply noise reduction to divined image."""
        if self.debug_mode.get():
            print("🧹 Applying noise reduction...")
        
        # Apply bilateral filter to reduce noise while preserving edges
        return cv2.bilateralFilter(image, 5, 50, 50)

    # === TEST FUNCTIONS ===
    
    def test_divine_process(self):
        """Test divine process with artificial holes."""
        if self.original_image is None:
            messagebox.showwarning("No Data", "Load an image first")
            return
        
        try:
            print("🧪 Testing G-S divine process with cascading...")
            
            # Create test image with artificial holes
            test_image = self.original_image.copy()
            height, width = test_image.shape[:2]
            
            # Create test mask (random holes)
            test_mask = np.zeros((height, width), dtype=bool)
            num_holes = 10
            hole_size = 15
            
            for i in range(num_holes):
                y_center = np.random.randint(hole_size, height - hole_size)
                x_center = np.random.randint(hole_size, width - hole_size)
                
                test_mask[y_center-hole_size//2:y_center+hole_size//2,
                         x_center-hole_size//2:x_center+hole_size//2] = True
                test_image[y_center-hole_size//2:y_center+hole_size//2,
                          x_center-hole_size//2:x_center+hole_size//2] = 0
            
            # Generate coordinate map for test
            if self.gs_coordinate_map is None:
                if self.gs_alignment_data is not None and self.use_gs_alignment_data.get():
                    self.generate_gs_coordinate_maps()
                else:
                    self.generate_coordinate_map()
            
            # Apply divine process
            method = self.divination_method.get()
            
            # Test cascading if enabled
            if self.enable_cascading.get() or method == "cascading_divine":
                print("   Testing cascading divine process...")
                divined_test = self.cascading_divine_process(test_image, test_mask)
            elif method == "gs_enhanced_divine":
                divined_test = self.gs_enhanced_divine_process(test_image, test_mask)
            elif method == "s_coordinate_flow":
                divined_test = self.s_coordinate_flow_divination(test_image, test_mask)
            elif method == "s_coordinate_similarity":
                divined_test = self.s_coordinate_similarity_divination(test_image, test_mask)
            elif method == "combined_divine":
                divined_test = self.combined_divine_process(test_image, test_mask)
            else:
                divined_test = self.omniscient_divine_process(test_image, test_mask)
            
            # Visualize test results
            self.visualize_test_results(self.original_image, test_image, divined_test, test_mask)
            
            print("🧪 G-S divine process test completed")
            
        except Exception as e:
            messagebox.showerror("Test Error", f"Test failed: {str(e)}")
            traceback.print_exc()

    # === VISUALIZATION METHODS (ALL PRESERVED) ===
    
    def show_loaded_image(self):
        """Show loaded image."""
        if self.original_image is None:
            return
            
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.imshow(self.original_image)
        ax.set_title(f"Loaded Image ({self.original_image.shape[1]}x{self.original_image.shape[0]})")
        ax.axis('off')
        
        info_text = ""
        if self.gs_alignment_data is not None:
            info_text += f"G-S Alignments: {len(self.gs_alignment_data)}\n"
        if self.mask_image is not None:
            info_text += f"Mask: {np.sum(self.mask_image)} pixels\n"
        if self.gs_coordinate_map is not None:
            source = "G-S Enhanced" if self.gs_coordinate_map.get('alignment_source', False) else "Generated"
            info_text += f"Coordinates: {source}\n"
        
        # Add threading info
        info_text += "Threading: ENABLED\n"
        
        # Add cascading info
        if self.enable_cascading.get():
            info_text += f"Cascading: {self.cascading_order.get()}\n"
            info_text += f"Batch Size: {self.cascading_batch_size.get()}"
        
        if info_text:
            ax.text(0.02, 0.98, info_text, transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.8))
        
        self.canvas.draw()

    def visualize_gs_coordinate_maps(self):
        """Visualize generated G-S coordinate maps."""
        if self.gs_coordinate_map is None:
            return
            
        self.fig.clear()
        
        # Create subplot layout
        if 'g_theta' in self.gs_coordinate_map:
            # Full G-S coordinate visualization
            gs = self.fig.add_gridspec(2, 3, hspace=0.4, wspace=0.3)
            
            ax1 = self.fig.add_subplot(gs[0, 0])
            im1 = ax1.imshow(self.gs_coordinate_map['g_theta'], cmap='hsv')
            ax1.set_title("G_theta (degrees)")
            ax1.axis('off')
            self.fig.colorbar(im1, ax=ax1, fraction=0.046)
            
            ax2 = self.fig.add_subplot(gs[0, 1])
            im2 = ax2.imshow(self.gs_coordinate_map['g_phi'], cmap='plasma')
            ax2.set_title("G_phi (degrees)")
            ax2.axis('off')
            self.fig.colorbar(im2, ax=ax2, fraction=0.046)
            
            ax3 = self.fig.add_subplot(gs[0, 2])
            im3 = ax3.imshow(self.gs_coordinate_map['s_x'], cmap='coolwarm')
            ax3.set_title("S_x coordinates")
            ax3.axis('off')
            self.fig.colorbar(im3, ax=ax3, fraction=0.046)
            
            ax4 = self.fig.add_subplot(gs[1, 0])
            im4 = ax4.imshow(self.gs_coordinate_map['s_y'], cmap='coolwarm')
            ax4.set_title("S_y coordinates")
            ax4.axis('off')
            self.fig.colorbar(im4, ax=ax4, fraction=0.046)
            
            # Show alignment points if available
            if self.gs_alignment_data is not None:
                ax5 = self.fig.add_subplot(gs[1, 1:])
                s_x_points = [a['S_x'] for a in self.gs_alignment_data]
                s_y_points = [a['S_y'] for a in self.gs_alignment_data]
                ax5.scatter(s_x_points, s_y_points, c='red', s=20, alpha=0.7)
                ax5.set_xlabel("S_x")
                ax5.set_ylabel("S_y")
                ax5.set_title(f"G-S Alignment Points ({len(self.gs_alignment_data)})")
                ax5.grid(True, alpha=0.3)
            
            source_text = "Generated from G-S alignment data" if self.gs_coordinate_map.get('alignment_source', False) else "Basic coordinate generation"
            self.fig.suptitle(f"G-S Compatible Coordinate Maps - {source_text} (THREADED)", fontsize=14, weight='bold')
        else:
            # Basic S-coordinate visualization
            ax1 = self.fig.add_subplot(121)
            im1 = ax1.imshow(self.gs_coordinate_map['s_x'], cmap='coolwarm')
            ax1.set_title("S_x coordinates")
            ax1.axis('off')
            self.fig.colorbar(im1, ax=ax1, fraction=0.046)
            
            ax2 = self.fig.add_subplot(122)
            im2 = ax2.imshow(self.gs_coordinate_map['s_y'], cmap='coolwarm')
            ax2.set_title("S_y coordinates")
            ax2.axis('off')
            self.fig.colorbar(im2, ax=ax2, fraction=0.046)
            
            self.fig.suptitle("Basic S-Coordinate Maps (THREADED)", fontsize=14, weight='bold')
        
        self.canvas.draw()

    def visualize_results(self):
        """Visualize divination results."""
        if self.divined_image is None:
            return
            
        self.fig.clear()
        
        if self.divine_outside_bounds.get():
            # Show original and expanded result
            gs = self.fig.add_gridspec(2, 2, hspace=0.3, wspace=0.2)
            
            ax1 = self.fig.add_subplot(gs[0, 0])
            ax1.imshow(self.original_image)
            ax1.set_title("Original Image")
            ax1.axis('off')
            
            ax2 = self.fig.add_subplot(gs[0, 1])
            ax2.imshow(self.divined_image)
            ax2.set_title("🔮 G-S Divined Image (Expanded)")
            ax2.axis('off')
            
            ax3 = self.fig.add_subplot(gs[1, :])
            ax3.imshow(self.divined_image)
            ax3.set_title("🔮 Full G-S Divine Result")
            ax3.axis('off')
        else:
            # Show side by side comparison
            ax1 = self.fig.add_subplot(121)
            ax1.imshow(self.original_image)
            ax1.set_title("Original Image")
            ax1.axis('off')
            
            ax2 = self.fig.add_subplot(122)
            ax2.imshow(self.divined_image)
            ax2.set_title("🔮 G-S Divined Image (THREADED)")
            ax2.axis('off')
        
        method_text = f"Method: {self.divination_method.get()}"
        if self.gs_alignment_data is not None and self.use_gs_alignment_data.get():
            method_text += f" • G-S Enhanced ({len(self.gs_alignment_data)} alignments)"
        if self.divine_outside_bounds.get():
            method_text += f" • Expanded: {self.edge_extension.get()}px"
        if self.enable_cascading.get():
            method_text += f" • Cascading: {self.cascading_order.get()}"
        method_text += " • THREADED"
        
        self.fig.suptitle(f"🔮 G-S Compatible Divine Pixel Results - {method_text}", 
                         fontsize=14, weight='bold')
        self.canvas.draw()

    def visualize_test_results(self, original, with_holes, divined, holes_mask):
        """Visualize test results."""
        self.fig.clear()
        
        gs = self.fig.add_gridspec(2, 2, hspace=0.3, wspace=0.2)
        
        # Original
        ax1 = self.fig.add_subplot(gs[0, 0])
        ax1.imshow(original)
        ax1.set_title("Original Image")
        ax1.axis('off')
        
        # With holes
        ax2 = self.fig.add_subplot(gs[0, 1])
        ax2.imshow(with_holes)
        ax2.set_title("With Test Holes")
        ax2.axis('off')
        
        # Divined
        ax3 = self.fig.add_subplot(gs[1, 0])
        ax3.imshow(divined)
        ax3.set_title("🔮 After G-S Divination")
        ax3.axis('off')
        
        # Difference
        ax4 = self.fig.add_subplot(gs[1, 1])
        difference = np.abs(original.astype(float) - divined.astype(float))
        difference_gray = np.mean(difference, axis=2)
        im4 = ax4.imshow(difference_gray, cmap='hot', vmin=0, vmax=50)
        ax4.set_title("Difference Map")
        ax4.axis('off')
        self.fig.colorbar(im4, ax=ax4, fraction=0.046)
        
        # Add method info
        method_text = f"Method: {self.divination_method.get()}"
        if self.gs_alignment_data is not None and self.use_gs_alignment_data.get():
            method_text += f"\nG-S Enhanced: {len(self.gs_alignment_data)} alignments"
        if self.enable_cascading.get():
            method_text += f"\nCascading: {self.cascading_order.get()}"
            method_text += f"\nGenerations: {self.cascading_stats.get('cascade_generations', 0)}"
        method_text += "\nTHREADED"
        
        ax3.text(0.02, 0.98, method_text, transform=ax3.transAxes, 
                verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.8))
        
        self.fig.suptitle("🧪 G-S Compatible Divine Process Test Results (THREADED)", fontsize=14, weight='bold')
        self.canvas.draw()

    # === STATISTICS AND STATUS ===
    
    def update_statistics(self):
        """Update divine enhancement statistics."""
        stats_text = "🔮 G-S DIVINE PIXEL STATISTICS + THREADING\n"
        stats_text += "=" * 40 + "\n"
        
        if self.original_image is not None:
            stats_text += f"Original: {self.original_image.shape}\n"
        
        if self.divined_image is not None:
            stats_text += f"Divined: {self.divined_image.shape}\n"
        
        if self.gs_alignment_data is not None and len(self.gs_alignment_data) > 0:
            stats_text += f"G-S Alignments: {len(self.gs_alignment_data)}\n"
            if self.use_gs_alignment_data.get():
                stats_text += "G-S Enhancement: ENABLED\n"
            else:
                stats_text += "G-S Enhancement: DISABLED\n"
            
            # Check for Hadit data
            sample = self.gs_alignment_data[0]
            if 'Hadit_theta' in sample and 'Hadit_phi' in sample:
                if self.hadit_enhancement.get():
                    stats_text += "Hadit Enhancement: ENABLED\n"
                else:
                    stats_text += "Hadit Enhancement: DISABLED\n"
            else:
                stats_text += "Hadit Data: NONE\n"
        else:
            stats_text += "G-S Alignments: NONE\n"
        
        # Threading info
        stats_text += "Threading: ENABLED\n"
        
        # Cascading statistics
        if self.enable_cascading.get():
            stats_text += f"Cascading: ENABLED ({self.cascading_order.get()})\n"
            if self.cascading_stats['total_pixels'] > 0:
                success_rate = (self.cascading_stats['successful_pixels'] / 
                              self.cascading_stats['total_pixels']) * 100
                stats_text += f"Cascade Success: {success_rate:.1f}%\n"
                stats_text += f"Generations: {self.cascading_stats['cascade_generations']}\n"
                stats_text += f"Batch Size: {self.cascading_batch_size.get()}\n"
        else:
            stats_text += "Cascading: DISABLED\n"
        
        stats_text += f"Method: {self.divination_method.get()}\n"
        stats_text += f"Divine Precision: {self.divine_precision.get():.3f}\n"
        stats_text += f"Search Radius: {self.divine_search_radius.get()}px\n"
        stats_text += f"Similarity Threshold: {self.divine_similarity_threshold.get():.1f}°\n"
        
        if self.gs_coordinate_map is not None:
            coord_source = "G-S Data" if self.gs_coordinate_map.get('alignment_source', False) else "Generated"
            stats_text += f"Coordinates: {coord_source}\n"
        
        stats_text += f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        
        self.stats_label.config(text=stats_text)

    def update_status(self):
        """Update status display."""
        if self.original_image is not None:
            self.image_status.config(text=f"✅ Image loaded: {self.original_image.shape}", 
                                   foreground="green")
            self.divine_button.config(state="normal")
            self.test_button.config(state="normal")
            self.status_label.config(text="Ready for G-S divine pixel process with threading & cascading", foreground="blue")
        else:
            self.image_status.config(text="❌ No image loaded", foreground="red")
        
        if self.gs_alignment_data is not None:
            self.alignment_status.config(text=f"✅ G-S data loaded: {len(self.gs_alignment_data)} alignments", 
                                       foreground="green")
            if self.original_image is not None:
                self.generate_gs_coords_button.config(state="normal")
        else:
            self.alignment_status.config(text="❌ No G-S alignment data loaded", foreground="orange")

    # === EXPORT METHODS ===
    
    def export_result(self):
        """Export divined image with G-S metadata and threading info."""
        if self.divined_image is None:
            messagebox.showwarning("No Data", "Perform divination first")
            return
        
        file_path = filedialog.asksaveasfilename(
            title="Save G-S Divined Image",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("JPEG files", "*.jpg"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                pil_image = Image.fromarray(self.divined_image.astype(np.uint8))
                pil_image.save(file_path)
                
                # Save enhanced metadata including G-S info and threading
                metadata_path = file_path.rsplit('.', 1)[0] + "_gs_threaded_metadata.txt"
                with open(metadata_path, 'w') as f:
                    f.write(f"🔮 G-S Compatible Divine Pixel Tool v1.3 + Threading Export\n")
                    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}\n")
                    f.write(f"Author: Angledcrystals\n")
                    f.write(f"Version: Divine Pixel Tool v1.3 (G-S Compatible with Threading & Cascading)\n\n")
                    
                    f.write(f"THREADING: ENABLED\n")
                    f.write(f"Background Processing: YES\n")
                    f.write(f"GUI Responsiveness: MAINTAINED\n")
                    f.write(f"Cancellation Support: YES\n\n")
                    
                    if self.original_image is not None:
                        f.write(f"Original Image: {self.original_image.shape}\n")
                    f.write(f"Divined Image: {self.divined_image.shape}\n")
                    f.write(f"Divination Method: {self.divination_method.get()}\n")
                    
                    # Cascading information
                    f.write(f"\nCASCADING DIVINATION:\n")
                    f.write(f"Cascading Enabled: {'YES' if self.enable_cascading.get() else 'NO'}\n")
                    if self.enable_cascading.get():
                        f.write(f"Cascading Order: {self.cascading_order.get()}\n")
                        f.write(f"Batch Size: {self.cascading_batch_size.get()}\n")
                        f.write(f"Max Distance: {self.cascading_max_distance.get()}\n")
                        f.write(f"Update Frequency: {self.cascading_update_frequency.get()}\n")
                        
                        if self.cascading_stats['total_pixels'] > 0:
                            success_rate = (self.cascading_stats['successful_pixels'] / 
                                          self.cascading_stats['total_pixels']) * 100
                            f.write(f"Total Pixels Processed: {self.cascading_stats['total_pixels']}\n")
                            f.write(f"Successfully Divined: {self.cascading_stats['successful_pixels']}\n")
                            f.write(f"Success Rate: {success_rate:.1f}%\n")
                            f.write(f"Cascade Generations: {self.cascading_stats['cascade_generations']}\n")
                    
                    if self.gs_alignment_data is not None and len(self.gs_alignment_data) > 0:
                        f.write(f"\nG-S ALIGNMENT DATA:\n")
                        f.write(f"Alignments Loaded: {len(self.gs_alignment_data)}\n")
                        f.write(f"G-S Enhancement: {'ENABLED' if self.use_gs_alignment_data.get() else 'DISABLED'}\n")
                        f.write(f"G-S Coordinate Precision: {self.gs_coordinate_precision.get():.3f}\n")
                        f.write(f"Hadit Enhancement: {'ENABLED' if self.hadit_enhancement.get() else 'DISABLED'}\n")
                    else:
                        f.write(f"\nG-S ALIGNMENT DATA: NONE\n")
                    
                    f.write(f"\nDIVINATION PARAMETERS:\n")
                    f.write(f"Divine Precision: {self.divine_precision.get():.3f}\n")
                    f.write(f"Search Radius: {self.divine_search_radius.get()}px\n")
                    f.write(f"Similarity Threshold: {self.divine_similarity_threshold.get():.1f}°\n")
                    f.write(f"Expansion Factor: {self.expansion_factor.get():.1f}\n")
                    
                    if self.divine_outside_bounds.get():
                        f.write(f"Edge Extension: {self.edge_extension.get()}px\n")
                        f.write("Divine Outside Bounds: ENABLED\n")
                    else:
                        f.write("Divine Outside Bounds: DISABLED\n")
                    
                    f.write(f"\nENABLED METHODS:\n")
                    f.write(f"Flow Analysis: {'✅' if self.use_flow_analysis.get() else '❌'}\n")
                    f.write(f"Similarity Matching: {'✅' if self.use_similarity_matching.get() else '❌'}\n")
                    f.write(f"Geometric Extrapolation: {'✅' if self.use_geometric_extrapolation.get() else '❌'}\n")
                    f.write(f"Gradient Synthesis: {'✅' if self.use_gradient_synthesis.get() else '❌'}\n")
                    f.write(f"Noise Reduction: {'✅' if self.noise_reduction.get() else '❌'}\n")
                    f.write(f"Cascading Divination: {'✅' if self.enable_cascading.get() else '❌'}\n")
                    f.write(f"Threading: ✅ ENABLED\n")
                    
                    if self.gs_coordinate_map is not None:
                        f.write(f"\nCOORDINATE MAP INFO:\n")
                        coord_source = "G-S Alignment Data" if self.gs_coordinate_map.get('alignment_source', False) else "Auto-Generated"
                        f.write(f"Source: {coord_source}\n")
                        f.write(f"Available Maps: {list(self.gs_coordinate_map.keys())}\n")
                
                messagebox.showinfo("Export Success", 
                                  f"🔮 G-S divined image with threading saved:\n{file_path}\n\nG-S metadata saved:\n{metadata_path}")
                
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to save image: {str(e)}")

    def export_maps(self):
        """Export G-S divine maps and coordinate data with threading info."""
        if self.gs_coordinate_map is None:
            messagebox.showwarning("No Data", "No G-S coordinate maps to export")
            return
        
        file_path = filedialog.asksaveasfilename(
            title="Save G-S Divine Maps",
            defaultextension=".npz",
            filetypes=[("NumPy files", "*.npz"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                # Include divination parameters and G-S data in export
                export_data = dict(self.gs_coordinate_map)
                export_data['divination_method'] = self.divination_method.get()
                export_data['divine_precision'] = self.divine_precision.get()
                export_data['divine_search_radius'] = self.divine_search_radius.get()
                export_data['divine_similarity_threshold'] = self.divine_similarity_threshold.get()
                export_data['expansion_factor'] = self.expansion_factor.get()
                export_data['edge_extension'] = self.edge_extension.get()
                
                # Threading info
                export_data['threading_enabled'] = True
                export_data['background_processing'] = True
                export_data['gui_responsive'] = True
                export_data['cancellation_support'] = True
                
                # Cascading parameters
                export_data['enable_cascading'] = self.enable_cascading.get()
                export_data['cascading_order'] = self.cascading_order.get()
                export_data['cascading_batch_size'] = self.cascading_batch_size.get()
                export_data['cascading_max_distance'] = self.cascading_max_distance.get()
                export_data['cascading_update_frequency'] = self.cascading_update_frequency.get()
                
                # Cascading statistics
                if self.cascading_stats['total_pixels'] > 0:
                    export_data['cascading_total_pixels'] = self.cascading_stats['total_pixels']
                    export_data['cascading_successful_pixels'] = self.cascading_stats['successful_pixels']
                    export_data['cascading_generations'] = self.cascading_stats['cascade_generations']
                    export_data['cascading_success_rate'] = (self.cascading_stats['successful_pixels'] / 
                                                            self.cascading_stats['total_pixels']) * 100
                
                # G-S alignment data info
                if self.gs_alignment_data is not None and len(self.gs_alignment_data) > 0:
                    export_data['gs_alignment_count'] = len(self.gs_alignment_data)
                    export_data['use_gs_alignment_data'] = self.use_gs_alignment_data.get()
                    export_data['gs_coordinate_precision'] = self.gs_coordinate_precision.get()
                    export_data['hadit_enhancement'] = self.hadit_enhancement.get()
                    
                    # Store sample alignment for reference
                    sample_alignment = self.gs_alignment_data[0]
                    for key, value in sample_alignment.items():
                        export_data[f'sample_gs_alignment_{key}'] = value
                
                # Processing options
                export_data['use_flow_analysis'] = self.use_flow_analysis.get()
                export_data['use_similarity_matching'] = self.use_similarity_matching.get()
                export_data['use_geometric_extrapolation'] = self.use_geometric_extrapolation.get()
                export_data['use_gradient_synthesis'] = self.use_gradient_synthesis.get()
                export_data['divine_outside_bounds'] = self.divine_outside_bounds.get()
                export_data['noise_reduction'] = self.noise_reduction.get()
                
                # Metadata
                export_data['export_timestamp'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')
                export_data['tool_version'] = 'Divine Pixel Tool v1.3 + Threading (G-S Compatible with Cascading)'
                export_data['author'] = 'Angledcrystals'
                export_data['cpu_cores'] = os.cpu_count()
                
                np.savez_compressed(file_path, **export_data)
                
                messagebox.showinfo("Export Success", f"🗺️ G-S divine maps with threading data saved:\n{file_path}")
                
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to save G-S maps: {str(e)}")

def main():
    """Main application entry point."""
    print("🔮 Starting G-S Compatible Divine Pixel Tool v1.3 + Threading")
    print("Advanced Image Pixel Divination System with G-S Stereo Viewer Compatibility")
    print("Now featuring Cascading Divination with Background Threading Support")
    print("Author: Angledcrystals")
    print("Date: 2025-06-10 10:46:18 UTC")
    print("=" * 80)
    
    try:
        root = tk.Tk()
        app = GSCompatibleDivinePixelTool(root)
        
        # Center window
        root.update_idletasks()
        width = root.winfo_width()
        height = root.winfo_height()
        x = (root.winfo_screenwidth() // 2) - (width // 2)
        y = (root.winfo_screenheight() // 2) - (height // 2)
        root.geometry(f"{width}x{height}+{x}+{y}")
        
        print("🔮 G-S Compatible Divine Pixel Tool with Threading launched successfully!")
        print("Ready for G-S enhanced threaded cascading pixel divination...")
        print("\n📊 G-S Compatibility Features:")
        print("   • Full compatibility with G-S Stereo Viewer alignment format")
        print("   • S_x, S_y, G_theta, G_phi coordinate support")  
        print("   • Hadit_theta, Hadit_phi vector enhancement")
        print("   • G-S enhanced divination methods")
        print("   • Complete metadata export with G-S alignment info")
        print("   • Identical coordinate generation methods as G-S Stereo Viewer")
        print("\n🧵 NEW Threading Features:")
        print("   • Background processing keeps GUI responsive during divination")
        print("   • Cancellable operations - stop processing at any time")
        print("   • Thread-safe progress updates for real-time feedback")
        print("   • Non-blocking interface during intensive processing")
        print("   • Graceful error handling across threads")
        print("\n🌊 Enhanced Cascading Features:")
        print("   • Sequential pixel processing using previously calculated pixels")
        print("   • Multiple cascading orders: nearest-first, flow-based, coordinate-based")
        print("   • Configurable batch processing and update frequencies")
        print("   • Real-time cascading progress visualization")
        print("   • Enhanced continuity and pattern propagation")
        print("   • Spiral-out and distance-weighted ordering algorithms")
        print("\n🚀 Threading Benefits:")
        print("   • Responsive GUI: Interface remains interactive during processing")
        print("   • Cancellable operations: Users can stop long-running processes")
        print("   • Real-time updates: Live progress feedback during divination")
        print("   • Better user experience: No more frozen interface")
        print("   • Maintained functionality: All original features preserved")
        print("\n🔧 Key Improvements:")
        print("   • MINIMAL CODE CHANGES: Only threading layer added")
        print("   • ALL ORIGINAL FUNCTIONALITY PRESERVED")
        print("   • Thread-safe cascading progress updates")
        print("   • Cancellation checks in processing loops")
        print("   • Background worker thread for intensive operations")
        print("   • Enhanced error handling and user feedback")
        
        # Add exception handling for better error reporting
        def show_error(exc_type, exc_value, exc_traceback):
            error_msg = ''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))
            print(f"❌ ERROR: {error_msg}")
            messagebox.showerror("Divine Pixel Tool Error", 
                              f"An error occurred:\n\n{exc_value}\n\nCheck console for details.")
            
        sys.excepthook = show_error
        
        root.mainloop()
        
    except Exception as e:
        print(f"❌ Failed to start G-S Compatible Divine Pixel Tool: {e}")
        traceback.print_exc()
        messagebox.showerror("Startup Error", f"Failed to start application: {str(e)}")
    
    print("🔮 G-S Compatible Divine Pixel Tool with Threading closed")

if __name__ == "__main__":
    main()
