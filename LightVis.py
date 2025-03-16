import sys
import time
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                            QLabel, QSlider, QPushButton, QComboBox, QGroupBox, QFrame,
                            QCheckBox, QTabWidget, QRadioButton, QButtonGroup, QSizePolicy)
from PyQt5.QtCore import Qt, QTimer, QSize
from PyQt5.QtGui import QPainter,QFont, QColor, QPen, QBrush
import pyqtgraph as pg

class PlayPauseButton(QPushButton):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setCheckable(True)
        self.setFixedSize(80, 40)
        self.setCursor(Qt.PointingHandCursor)
        self.setStyleSheet("""
            QPushButton {
                border-radius: 10px;
                border: none;
            }
            QPushButton:hover {
                background-color: #1C97EA;
            }
            QPushButton:checked {
                background-color: #007ACC;
            }
            QPushButton:checked:hover {
                background-color: #1C97EA;
            }
            QPushButton:!checked {
                background-color: #E74C3C;
            }
            QPushButton:!checked:hover {
                background-color: #FF6B6B;
            }
        """)
        self.update_icon(False)
        
    def update_icon(self, paused):
        if paused:
            # Play icon (triangle)
            self.setToolTip("Play Animation")
            self.setText("▶")
        else:
            # Pause icon (two vertical bars)
            self.setToolTip("Pause Animation")
            self.setText("||")
        
        # Set text color and font
        self.setStyleSheet(self.styleSheet() + """
            QPushButton {
                color: white;
                font-size: 18px;
                font-weight: bold;
            }
        """)

# Custom colored slider for wavelength selection
class ColoredSlider(QSlider):
    def __init__(self, parent=None):
        super().__init__(Qt.Horizontal, parent)
        self.setStyleSheet("""
            QSlider::groove:horizontal {
                border: 1px solid #5c5c5c;
                height: 8px;
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0.000 #ff0000,    /* 400 THz - red */
                    stop:0.200 #ff8800,    /* 460 THz - orange */
                    stop:0.350 #ffff00,    /* 510 THz - yellow */
                    stop:0.500 #00ff00,    /* 540 THz - green */
                    stop:0.650 #0088ff,    /* 580 THz - light blue */
                    stop:0.800 #4400ff,    /* 680 THz - deep blue */
                    stop:1.000 #8800ff);   /* 790 THz - violet */
                margin: 2px 0;
                border-radius: 5px;
            }
            QSlider::handle:horizontal {
                background: white;
                border: 1px solid #999999;
                width: 18px;
                height: 18px;
                margin: -6px 0;
                border-radius: 9px;
            }
        """)

# Define frequency to RGB color mapping
def wavelength_to_rgb(wavelength):
    """Convert wavelength (in nm) to RGB color values"""
    # Ensure wavelength is within visible spectrum (approximately 380-750 nm)
    if wavelength < 380:
        wavelength = 380  # Cap at violet
    elif wavelength > 750:
        wavelength = 750  # Cap at red
    
    
    # Convert to RGB using a simplified approximation of visible spectrum
    if 380 <= wavelength < 440:
        # Violet to blue
        r = ((440 - wavelength) / (440 - 380)) * 0.8
        g = 0.0
        b = 1.0
    elif 440 <= wavelength < 490:
        # Blue to cyan
        r = 0.0
        g = ((wavelength - 440) / (490 - 440))
        b = 1.0
    elif 490 <= wavelength < 510:
        # Cyan to green
        r = 0.0
        g = 1.0
        b = ((510 - wavelength) / (510 - 490))
    elif 510 <= wavelength < 580:
        # Green to yellow
        r = ((wavelength - 510) / (580 - 510))
        g = 1.0
        b = 0.0
    elif 580 <= wavelength < 645:
        # Yellow to red
        r = 1.0
        g = ((645 - wavelength) / (645 - 580))
        b = 0.0
    else:  # 645-750
        # Red
        r = 1.0
        g = 0.0
        b = 0.0
    
    # Scale RGB values to 0-255 range
    r = int(r * 255)
    g = int(g * 255)
    b = int(b * 255)
    
    return r, g, b

def frequency_to_rgb(frequency_THz):
    """Convert frequency in THz to RGB color"""
    # Convert frequency to wavelength in nm
    # c = f*λ, where c is speed of light (299,792,458 m/s)
    # λ(nm) = 299792.458 / f(THz)
    wavelength_nm = 299792.458 / frequency_THz
    
    # Now use the existing wavelength to RGB conversion
    return wavelength_to_rgb(wavelength_nm)

class WaveSimulationWidget(pg.PlotWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.angle_of_incidence = 0
        self.time = 0.0
        self.visualization_scale = 0.1
        # Initialize colors for mediums
        self.medium1_color = '#22222259'  
        self.medium2_color = '#1E90FF59'  
        self.medium3_color = '#88DDFF80'  
        
        self.medium1_name = 'Air'
        self.medium2_name = 'Water'
        self.medium3_name = 'Glass (Crown)'

        self.medium_presets = {}
        # Initialize parameters
        self.frequency = 545  # nm
        self.amplitude = 1.0
        self.speed = 1.0
        self.angle = 0.0  # degrees
        self.ray_target = 0.0
        self.show_interference = False
        self.show_ray_mode = False
        self.white_light = False
        self.paused = False
        self.prism_mode = False
        self.superposition_enabled = False
        self.superposition_wave = None

        self.prism_frequencies = np.linspace(400, 790, 10)
        self.prism_wavelengths = [299792.458 / freq for freq in self.prism_frequencies]
        # Initialize medium properties
        self.n1 = 1.0  
        self.n2 = 1.33  
        self.n3 = 1.5  
        
        # Initialize boundaries
        self.boundary1 = 1000  # First boundary position
        self.boundary2 = 2000  # Second boundary position
        
        # Create the plot with a transparent background
        self.setBackground(None)  
        self.getPlotItem().setTitle("")  
        self.getPlotItem().hideAxis('left') 
        self.getPlotItem().hideAxis('bottom')  

        # Configure a thin, visible grid
        grid_pen = pg.mkPen(color=(255, 255, 255, 80), width=0.5, style=Qt.DotLine)
        
        self.getPlotItem().getAxis('left').setPen(grid_pen)
        self.getPlotItem().getAxis('bottom').setPen(grid_pen)
        
        
        # Disable mouse interaction and context menu
        self.setMouseEnabled(x=False, y=False)
        self.setMenuEnabled(False)
        
        # Set fixed view range
        self.setXRange(0, 3000, padding=0)
        self.setYRange(-2, 2, padding=0)
        
        # Initialize lists for plot items
        self.medium_rects = []
        self.medium_labels = []
        
        # Create wave curve
        self.x = np.linspace(0, 3000, 6000)
 
        
        # Initialize all container lists first
        self.wave_curves = []   # Initialize wave_curves list
        

        # Set up the plot
        self.setAntialiasing(True)

        # Create custom grid lines that will be drawn on top
        self.grid_lines_x = []
        self.grid_lines_y = []
        
       # Horizontal gridlines
        for y in np.arange(-2, 3, 1):
            line = pg.InfiniteLine(pos=y, angle=0, pen=pg.mkPen(color=(200, 200, 200, 150), width=0.5, style=Qt.DotLine))
            self.addItem(line)
            self.grid_lines_y.append(line)
        
        # Vertical gridlines
        for x in np.arange(0, 3001, 500):
            line = pg.InfiniteLine(pos=x, angle=90, pen=pg.mkPen(color=(200, 200, 200, 150), width=0.5, style=Qt.DotLine))
            self.addItem(line)
            self.grid_lines_x.append(line)
        # Remove the plot border
        self.getPlotItem().setContentsMargins(0, 0, 0, 0)
        self.getPlotItem().layout.setContentsMargins(0, 0, 0, 0)
        self.getPlotItem().vb.setContentsMargins(0, 0, 0, 0)

        # Remove the border and axis rectangle
        self.getPlotItem().getViewBox().setBackgroundColor(None)
        self.getPlotItem().getViewBox().setBorder(None)
        
        # Hide the axis box and remove borders
        for axis in ['left', 'bottom', 'top', 'right']:
            if self.getPlotItem().getAxis(axis):
                self.getPlotItem().getAxis(axis).setStyle(showValues=True, tickLength=5)
                if axis in ['top', 'right']:
                    self.getPlotItem().showAxis(axis, False)
        


        # Set grid line spacing
        self.getPlotItem().getAxis('left').setTicks([[(i, str(i)) for i in range(-2, 3, 1)]])
        self.getPlotItem().getAxis('bottom').setTicks([[(i, str(i)) for i in range(0, 3001, 500)]])

        # Set up the x-axis
        self.x = np.linspace(0, 3000, 3000)

        # Create boundary lines early
        self.boundary1_line = pg.InfiniteLine(pos=self.boundary1, angle=90, pen=pg.mkPen('w', width=2, style=Qt.DashLine))
        self.boundary2_line = pg.InfiniteLine(pos=self.boundary2, angle=90, pen=pg.mkPen('w', width=2, style=Qt.DashLine))
        self.addItem(self.boundary1_line)
        self.addItem(self.boundary2_line)
        
        self.medium_rects = []  # Initialize as empty list
        self.medium_rects.append(self.boundary1_line)
        self.medium_rects.append(self.boundary2_line)
        

        # Initial parameters
        self.wavelength = 550
        self.amplitude = 5
        self.speed = 4
        self.n1 = 1.0003  # Air
        self.n2 = 1.33    # Water
        self.n3 = 1.52    # Glass (Crown)
        self.time = 0
        self.angle_of_incidence = 0  # Default angle is 0 (straight line)
        self.ray_target_y = 0  # Y-coordinate where ray hits first boundary
        
        # Initialize flags
        self.show_interference = False
        self.show_ray_mode = False
        self.white_light = False
        self.prism_mode = False

        
               # Create angle labels with improved styling
        label_html_style = """<div style="font-family: Arial; font-size: 16pt; font-weight: bold; color: white;">"""
        self.angle1_label = pg.TextItem(html=f"{label_html_style}θ₁: 0°</div>", anchor=(0.5, 0.5))
        self.angle2_label = pg.TextItem(html=f"{label_html_style}θ₂: 0°</div>", anchor=(0.5, 0.5))
        self.angle3_label = pg.TextItem(html=f"{label_html_style}θ₃: 0°</div>", anchor=(0.5, 0.5))
        self.addItem(self.angle1_label)
        self.addItem(self.angle2_label)
        self.addItem(self.angle3_label)
        
        # Position the angle labels in the center of each medium
        self.angle1_label.setPos(self.boundary1/2, -1.5)
        self.angle2_label.setPos(self.boundary1 + (self.boundary2-self.boundary1)/2, -1.5)
        self.angle3_label.setPos(self.boundary2 + (3000-self.boundary2)/2, -1.5)

        # Add interference wave curves with specific colors
        self.interference_wave1 = self.plot(self.x, np.zeros_like(self.x), 
                                          pen=pg.mkPen('r', width=2, style=Qt.DashLine))
        self.interference_wave2 = self.plot(self.x, np.zeros_like(self.x), 
                                          pen=pg.mkPen('g', width=2, style=Qt.DashLine))
        
        wave = self.calculate_wave(self.wavelength)
        r, g, b = wavelength_to_rgb(self.wavelength)
        wave_color = QColor(r, g, b)
        self.wave_curve = self.plot(self.x, wave, pen=pg.mkPen(wave_color, width=4))

        # Initialize interference waves with proper data
        wave1, wave2 = self.calculate_interference_waves(self.wavelength)
        self.interference_wave1.setData(self.x, wave1)
        self.interference_wave2.setData(self.x, wave2)
        
        # Hide interference waves initially
        self.interference_wave1.setVisible(False)
        self.interference_wave2.setVisible(False)
        self.update_plot()
        
        # Create ray lines to visualize refraction directions
        self.ray_incident = pg.PlotDataItem([], [], pen=pg.mkPen('#ffffff', width=4, style=Qt.DashLine))
        self.ray_refracted1 = pg.PlotDataItem([], [], pen=pg.mkPen('#00aaff', width=4, style=Qt.DashLine))
        self.ray_refracted2 = pg.PlotDataItem([], [], pen=pg.mkPen('#22ff22', width=4, style=Qt.DashLine))
        
        # Add reflection markers
        self.reflection_marker1 = pg.ScatterPlotItem(size=12, brush='y', symbol='o')
        self.reflection_marker2 = pg.ScatterPlotItem(size=12, brush='y', symbol='o')
        
        # Add ray lines to plot but keep them hidden initially
        self.addItem(self.ray_incident)
        self.addItem(self.ray_refracted1)
        self.addItem(self.ray_refracted2)
        self.addItem(self.reflection_marker1)
        self.addItem(self.reflection_marker2)
        self.ray_incident.setVisible(False)
        self.ray_refracted1.setVisible(False)
        self.ray_refracted2.setVisible(False)
        self.reflection_marker1.setVisible(False)
        self.reflection_marker2.setVisible(False)

 
        # Add frame rate control
        self.last_update_time = time.time()
        self.frame_time = 1/60  # Target 60 FPS
        self.skip_frames = 0


        # Set up the animation timer
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_animation)
        self.timer.start(16)  # 16ms interval (60 fps)

        self.update_wavelength_scale()
        
    def calculate_interference_waves(self, wavelength):
        """Calculate the component waves that create interference"""
        if not self.show_interference:
            return np.zeros_like(self.x), np.zeros_like(self.x)
        
        # Use reduced resolution for interference waves
        # Take every 3rd point for better performance
        reduced_x = self.x[::3]
        
        # Get the actual refracted wave with reduced points
        actual_wave = self.calculate_wave(self.frequency)
        reduced_actual = actual_wave[::3]
        
        # Wave 1 (red) - original wave as if it was coming from vacuum (n=1.0)
        # Calculate vacuum wavelength
        vacuum_wavelength = 299792.458 / self.frequency
        
        # Calculate wave number for vacuum (n=1.0)
        k_vacuum = 2 * np.pi / vacuum_wavelength
        
        # Calculate phase with time component to ensure animation
        # Use the same speed scaling as in calculate_wave
        omega_vacuum = 2 * np.pi * self.frequency * (self.speed / 50) / 100
        phase_vacuum = omega_vacuum * self.time
        
        # For all media, calculate the wave as if it was in vacuum (n=1.0)
        # Use reduced resolution
        reduced_wave1 = self.amplitude * self.visualization_scale * 1.5 * np.sin(k_vacuum * reduced_x - phase_vacuum)
        
        # Wave 2 (green) - the difference between actual wave and vacuum wave
        # This represents the effect of the media on the wave
        reduced_wave2 = reduced_actual - reduced_wave1
        
        return reduced_x, reduced_wave1, reduced_wave2
        
    def calculate_superposition_interference(self):
        """Calculate interference waves based on the superposition of all wavelengths"""

        # Create cache key based on current state
        cache_key = (self.time, self.n1, self.n2, self.n3)
        # Check if we have this calculation cached
        if hasattr(self, '_interference_cache') and cache_key in self._interference_cache:
            return self._interference_cache[cache_key]
        # Use reduced resolution for interference waves
        reduced_x = self.x[::3]
        
        # Start with zeros
        superposition = np.zeros_like(reduced_x)
        
        # Add all individual waves to create the superposition
        for curve, freq in self.wave_curves:
            wave = self.calculate_wave(freq)
            reduced_wave = wave[::3]  # Take every 3rd point
            superposition += reduced_wave
            
        # Scale the superposition to keep it within reasonable amplitude
        if len(self.wave_curves) > 0:
            superposition = superposition / len(self.wave_curves)
        
        # For interference visualization, we need two components:
        # 1. The original wave (red) - this should be the superposition of vacuum waves
        vacuum_superposition = np.zeros_like(reduced_x)
        
        # Calculate a reference wave with n=1.0 for all media
        for curve, freq in self.wave_curves:
            # Calculate vacuum wavelength
            vacuum_wavelength = 299792.458 / freq
            
            # Calculate wave number for vacuum
            k_vacuum = 2 * np.pi / vacuum_wavelength
            
            # Calculate phase with time component
            # Use the same speed scaling as in calculate_wave
            omega_vacuum = 2 * np.pi * freq * (self.speed / 50) / 100
            phase_vacuum = omega_vacuum * self.time
            
            # Calculate vacuum wave with reduced resolution
            vacuum_wave = self.amplitude * self.visualization_scale * 1.5 * np.sin(k_vacuum * reduced_x - phase_vacuum)
            
            # Add to vacuum superposition
            vacuum_superposition += vacuum_wave
        
        # Scale the vacuum superposition
        if len(self.wave_curves) > 0:
            vacuum_superposition = vacuum_superposition / len(self.wave_curves)
        
        # Wave 1 (red) - the vacuum superposition
        wave1 = vacuum_superposition
        
        # Wave 2 (green) - the difference between actual and vacuum superposition
        # This represents the effect of the media on the wave
        wave2 = superposition - vacuum_superposition
            
        # Cache the result
        if not hasattr(self, '_interference_cache'):
            self._interference_cache = {}
        self._interference_cache[cache_key] = (reduced_x, wave1, wave2)
        
        # Limit cache size
        if len(self._interference_cache) > 50:
            oldest_keys = list(self._interference_cache.keys())[:10]
            for key in oldest_keys:
                self._interference_cache.pop(key)

        return reduced_x, wave1, wave2

    def calculate_wave(self, frequency, additional_phase=0):
        """Calculate wave values for the given frequency"""
        # Convert frequency to wavelength in nm
        wavelength_nm = 299792.458 / frequency
        
        # Use cache key that includes additional_phase
        cache_key = (frequency, self.time, self.n1, self.n2, self.n3, 
                     self.amplitude, self.angle_of_incidence, additional_phase)
        
        # Check if we have this calculation cached
        if hasattr(self, '_wave_calc_cache') and cache_key in self._wave_calc_cache:
            return self._wave_calc_cache[cache_key]
        
        # Initialize cache if it doesn't exist
        if not hasattr(self, '_wave_calc_cache'):
            self._wave_calc_cache = {}
        
        # Calculate wavelengths in each medium (nm)
        wavelength_m1 = wavelength_nm / self.n1
        wavelength_m2 = wavelength_nm / self.n2
        wavelength_m3 = wavelength_nm / self.n3
        
        # Calculate wave numbers (radians per nm)
        k1 = 2 * np.pi / wavelength_m1
        k2 = 2 * np.pi / wavelength_m2
        k3 = 2 * np.pi / wavelength_m3
        
        # Calculate angular frequency (radians per time unit)
        omega = 2 * np.pi * frequency * (self.speed / 50) / 100
        
        # Calculate phase at each point
        phase = omega * self.time + additional_phase
        
        # Pre-allocate array for better performance
        y = np.zeros_like(self.x)
        
        # Apply proper scaling to amplitude
        scaled_amplitude = self.amplitude * self.visualization_scale * 1.5 
        
        # Vectorized calculation for each region
        mask1 = self.x < self.boundary1
        mask2 = (self.x >= self.boundary1) & (self.x < self.boundary2)
        mask3 = self.x >= self.boundary2
        
        # Calculate phase shifts
        phase_shift2 = (k1 - k2) * self.boundary1
        phase_shift3 = (k1 - k2) * self.boundary1 + (k2 - k3) * self.boundary2
        
        # Vectorized calculations (much faster than loops)
        y[mask1] = scaled_amplitude * np.sin(k1 * self.x[mask1] - phase)
        y[mask2] = scaled_amplitude * np.sin(k2 * self.x[mask2] - phase + phase_shift2)
        y[mask3] = scaled_amplitude * np.sin(k3 * self.x[mask3] - phase + phase_shift3)
        
        # Cache the result
        self._wave_calc_cache[cache_key] = y
        
        # Limit cache size to prevent memory issues
        if len(self._wave_calc_cache) > 100:
            # Remove oldest items first (more efficient)
            keys = list(self._wave_calc_cache.keys())
            oldest_keys = keys[:10]
            for key in oldest_keys:
                self._wave_calc_cache.pop(key)
        
        return y

    
    def update_animation(self):
        """Update the animation for each timer tick"""
        # Implement frame skipping for performance
        current_time = time.time()
        elapsed = current_time - self.last_update_time

        if elapsed > 0.033 and self.skip_frames < 2:  # 33ms = ~30 FPS
            self.skip_frames += 1
            return

        self.skip_frames = 0
        self.last_update_time = current_time

        if not self.paused:  # Only update time if not paused
            self.time += 0.01
        
        # Only update UI if something has changed
        needs_update = False
        
        # Store previous state to check for changes
        if not hasattr(self, '_prev_state'):
            self._prev_state = {
                'time': self.time - 1,  # Force initial update
                'white_light': self.white_light,
                'show_interference': self.show_interference,
                'frequency': self.frequency,
                'n1': self.n1,
                'n2': self.n2,
                'n3': self.n3,
                'superposition_enabled': self.superposition_enabled
            }
            needs_update = True
        
        # Check if any relevant state has changed
        if (abs(self._prev_state['time'] - self.time) > 0.009 or
            self._prev_state['white_light'] != self.white_light or
            self._prev_state['show_interference'] != self.show_interference or
            self._prev_state['frequency'] != self.frequency or
            self._prev_state['superposition_enabled'] != self.superposition_enabled or
            abs(self._prev_state['n1'] - self.n1) > 0.0001 or
            abs(self._prev_state['n2'] - self.n2) > 0.0001 or
            abs(self._prev_state['n3'] - self.n3) > 0.0001):
            needs_update = True
        
        if needs_update:
            # For white light mode, optimize updates
            if self.white_light:
                # Only update wavelength scale once, using the main frequency
                if (self._prev_state['frequency'] != self.frequency or
                    abs(self._prev_state['n1'] - self.n1) > 0.0001 or
                    abs(self._prev_state['n2'] - self.n2) > 0.0001 or
                    abs(self._prev_state['n3'] - self.n3) > 0.0001):
                    self.update_wavelength_scale()
                
                if self.superposition_enabled:
                    # Update superposition wave
                    self.update_superposition_wave()
                else:
                    # Only update component waves if they're visible
                    visible_components = False
                    for curve_item in self.wave_curves:
                        if isinstance(curve_item, tuple) and len(curve_item) == 2:
                            curve, _ = curve_item
                            if hasattr(curve, 'setVisible') and curve.isVisible():
                                visible_components = True
                                break
                    
                    if visible_components:
                        # Batch update for white light mode with reduced processing
                        self.setUpdatesEnabled(False)  # Disable updates while making changes
                        
                        # Pre-calculate all waves with reduced data
                        reduced_x = self.x[::20]  # Take every 20th point
                        
                        for curve, freq in self.wave_curves:
                            if curve.isVisible():
                                wave = self.calculate_wave(freq)
                                reduced_wave = wave[::20]  # Take every 20th point
                                curve.setData(reduced_x, reduced_wave)
                        
                        self.setUpdatesEnabled(True)  # Re-enable updates after all changes
            else:
                # Update single wave
                wave = self.calculate_wave(self.frequency)
                self.wave_curve.setData(self.x, wave)
                
                # Update wavelength scale for single wave mode
                if (self._prev_state['frequency'] != self.frequency or
                    abs(self._prev_state['n1'] - self.n1) > 0.0001 or
                    abs(self._prev_state['n2'] - self.n2) > 0.0001 or
                    abs(self._prev_state['n3'] - self.n3) > 0.0001):
                    self.update_wavelength_scale()
            
            # Update interference waves if enabled - use reduced resolution
            if self.show_interference and needs_update:
                if self.white_light:
                    # Always use superposition-based interference for white light mode
                    reduced_x, wave1, wave2 = self.calculate_superposition_interference()
                else:
                    # Use normal interference with reduced resolution for single wave mode
                    reduced_x, wave1, wave2 = self.calculate_interference_waves(self.wavelength)
                
                # Update with reduced data
                self.interference_wave1.setData(reduced_x, wave1)
                self.interference_wave2.setData(reduced_x, wave2)
            
            # Update previous state
            self._prev_state = {
                'time': self.time,
                'white_light': self.white_light,
                'show_interference': self.show_interference,
                'frequency': self.frequency,
                'n1': self.n1,
                'n2': self.n2,
                'n3': self.n3,
                'superposition_enabled': self.superposition_enabled
            }

    def toggle_pause(self, paused):
        """Toggle pause state"""
        self.paused = paused
        if hasattr(self, 'timer'):
            if paused:
                self.timer.stop()
            else:
                self.timer.start(16)



    def toggle_interference(self, enabled):
        """Toggle visibility of interference waves"""
        self.show_interference = enabled

        # Set visibility based on enabled state
        self.interference_wave1.setVisible(enabled)
        self.interference_wave2.setVisible(enabled)
        
        if enabled:
            # Update pen colors to make them more visible
            self.interference_wave1.setPen(pg.mkPen('r', width=2, style=Qt.DashLine))
            self.interference_wave2.setPen(pg.mkPen('g', width=2, style=Qt.DashLine))
            
            # Force a significant time update to ensure waves are animated
            old_time = self.time
            self.time += 0.1  # Add a time offset
            
            # Force recalculation with the new time
            if self.white_light:
                # Always use superposition-based interference for white light
                reduced_x, wave1, wave2 = self.calculate_superposition_interference()
            else:
                # Calculate normal interference with reduced resolution for single wave
                reduced_x, wave1, wave2 = self.calculate_interference_waves(self.wavelength)
            
            # Update the interference waves with reduced data
            self.interference_wave1.setData(reduced_x, wave1)
            self.interference_wave2.setData(reduced_x, wave2)
            
            # Reset time to original plus a small increment to keep animation flowing
            self.time = old_time + 0.01
        else: 
            self.interference_wave1.setVisible(False)
            self.interference_wave2.setVisible(False)

    def toggle_superposition(self, enabled):
        """Toggle between showing individual waves or a single superposition wave"""
        if self.white_light:
            # In white light mode, invert the behavior
            # enabled=True means "Show Components"
            # enabled=False means "Show Superposition" (default)
            self.superposition_enabled = not enabled
            
            # Create superposition wave if it doesn't exist
            if self.superposition_wave is None:
                self.superposition_wave = self.plot(self.x, np.zeros_like(self.x), 
                                                   pen=pg.mkPen('w', width=4))
            
            # Show/hide appropriate waves
            if self.superposition_wave:
                self.superposition_wave.setVisible(not enabled)
            
            # If showing components, recreate them with reduced data
            if enabled:
                # Recreate with reduced data for better performance
                self.create_white_light_curves()
                
                # Now show the components
                for curve_item in self.wave_curves:
                    if isinstance(curve_item, tuple) and len(curve_item) == 2:
                        curve, _ = curve_item
                        if hasattr(curve, 'setVisible'):
                            curve.setVisible(True)
                self._update_component_visibility(enabled)            
            else:
                # Hide all component waves
                for curve_item in self.wave_curves:
                    if isinstance(curve_item, tuple) and len(curve_item) == 2:
                        curve, _ = curve_item
                        if hasattr(curve, 'setVisible'):
                            curve.setVisible(False)
                
                # Update the superposition wave
                self.update_superposition_wave()
        else:
            # In normal mode, keep original behavior
            self.superposition_enabled = enabled
            
            # Create superposition wave if it doesn't exist and enabled
            if enabled and self.superposition_wave is None:
                self.superposition_wave = self.plot(self.x, np.zeros_like(self.x), 
                                                   pen=pg.mkPen('w', width=4))
            
            # Show/hide superposition wave
            if self.superposition_wave:
                self.superposition_wave.setVisible(enabled)
            
            # Update the superposition wave if enabled
            if enabled:
                self.update_superposition_wave()

    def update_plot(self):
        """Update the plot with current parameters"""
        try:
            # Defer updates until all changes are made
            self.setUpdatesEnabled(False)
            
            # Store current y range to restore after updates
            y_range = self.getViewBox().viewRange()[1]
            
            # Initialize medium rectangles and labels if they don't exist
            if not hasattr(self, '_plot_initialized'):
                self._create_initial_plot_elements()
                self._plot_initialized = True
            
            # Update medium rectangles instead of recreating them
            self._update_medium_rectangles()
            
            # Update medium labels
            self._update_medium_labels()
            
            # Update the wave curve with current parameters
            if not self.white_light:
                # Get color from wavelength
                r, g, b = wavelength_to_rgb(299792.458 / self.frequency)
                wave_color = QColor(r, g, b)
                self.wave_curve.setPen(pg.mkPen(wave_color, width=4))
                
                # Update wave data
                wave = self.calculate_wave(self.frequency)
                self.wave_curve.setData(self.x, wave)
            else:
                # Batch update white light curves
                # Pre-calculate all waves
                all_waves = []
                for curve, freq in self.wave_curves:
                    if curve.isVisible():  # Only update visible curves
                        wave = self.calculate_wave(freq)
                        all_waves.append((curve, wave))
                
                # Update all curves at once
                for curve, wave in all_waves:
                    curve.setData(self.x, wave)
            
            # Update ray lines if ray mode is enabled
            if self.show_ray_mode:
                self.update_ray_lines()
                
            # Restore y range to prevent auto-scaling
            self.setYRange(y_range[0], y_range[1], padding=0)
            
            # Re-enable updates
            self.setUpdatesEnabled(True)
            
        except Exception as e:
            print(f"Error in update_plot: {str(e)}")
            self.setUpdatesEnabled(True)  # Make sure updates are re-enabled
            return
            
    def _create_initial_plot_elements(self):
        """Create initial plot elements that will be reused"""
        # Create medium rectangles
        self.medium1_rect = pg.FillBetweenItem(
            pg.PlotCurveItem([0, self.boundary1], [2, 2]),
            pg.PlotCurveItem([0, self.boundary1], [-2, -2]),
            brush=pg.mkBrush(self.medium1_color))
        self.addItem(self.medium1_rect)
        
        self.medium2_rect = pg.FillBetweenItem(
            pg.PlotCurveItem([self.boundary1, self.boundary2], [2, 2]),
            pg.PlotCurveItem([self.boundary1, self.boundary2], [-2, -2]),
            brush=pg.mkBrush(self.medium2_color))
        self.addItem(self.medium2_rect)
        
        self.medium3_rect = pg.FillBetweenItem(
            pg.PlotCurveItem([self.boundary2, 3000], [2, 2]),
            pg.PlotCurveItem([self.boundary2, 3000], [-2, -2]),
            brush=pg.mkBrush(self.medium3_color))
        self.addItem(self.medium3_rect)
        
        # Create medium labels
        self.medium1_color_label = pg.TextItem("", anchor=(0.5, 0), color='white')
        self.medium1_name_label = pg.TextItem("", anchor=(0.5, 1), color='white')
        self.addItem(self.medium1_color_label)
        self.addItem(self.medium1_name_label)
        
        self.medium2_color_label = pg.TextItem("", anchor=(0.5, 0), color='white')
        self.medium2_name_label = pg.TextItem("", anchor=(0.5, 1), color='white')
        self.addItem(self.medium2_color_label)
        self.addItem(self.medium2_name_label)
        
        self.medium3_color_label = pg.TextItem("", anchor=(0.5, 0), color='white')
        self.medium3_name_label = pg.TextItem("", anchor=(0.5, 1), color='white')
        self.addItem(self.medium3_color_label)
        self.addItem(self.medium3_name_label)
        
    def _update_medium_rectangles(self):
        """Update medium rectangles without recreating them"""
        # Update medium 1
        self.medium1_rect.setCurves(
            pg.PlotCurveItem([0, self.boundary1], [2, 2]),
            pg.PlotCurveItem([0, self.boundary1], [-2, -2]))
        self.medium1_rect.setBrush(pg.mkBrush(self.medium1_color))
        
        # Update medium 2
        self.medium2_rect.setCurves(
            pg.PlotCurveItem([self.boundary1, self.boundary2], [2, 2]),
            pg.PlotCurveItem([self.boundary1, self.boundary2], [-2, -2]))
        self.medium2_rect.setBrush(pg.mkBrush(self.medium2_color))
        
        # Update medium 3
        self.medium3_rect.setCurves(
            pg.PlotCurveItem([self.boundary2, 3000], [2, 2]),
            pg.PlotCurveItem([self.boundary2, 3000], [-2, -2]))
        self.medium3_rect.setBrush(pg.mkBrush(self.medium3_color))
        
        # Update boundary lines
        self.boundary1_line.setValue(self.boundary1)
        self.boundary2_line.setValue(self.boundary2)
        
    def _update_medium_labels(self):
        """Update medium labels without recreating them"""
        # Update medium 1 labels
        self.medium1_color_label.setText(f"n₁: {self.n1:.4f}")
        self.medium1_color_label.setPos(self.boundary1 / 2, -1.2)
        self.medium1_name_label.setText(f"{self.medium1_name}")
        self.medium1_name_label.setPos(self.boundary1 / 2, 1.8)
        
        # Update medium 2 labels
        self.medium2_color_label.setText(f"n₂: {self.n2:.4f}")
        self.medium2_color_label.setPos(self.boundary1 + (self.boundary2 - self.boundary1) / 2, -1.2)
        self.medium2_name_label.setText(f"{self.medium2_name}")
        self.medium2_name_label.setPos(self.boundary1 + (self.boundary2 - self.boundary1) / 2, 1.8)
        
        # Update medium 3 labels
        self.medium3_color_label.setText(f"n₃: {self.n3:.4f}")
        self.medium3_color_label.setPos(self.boundary2 + (3000 - self.boundary2) / 2, -1.2)
        self.medium3_name_label.setText(f"{self.medium3_name}")
        self.medium3_name_label.setPos(self.boundary2 + (3000 - self.boundary2) / 2, 1.8)

    def toggle_ray_mode(self, enabled):
        """Toggle ray visualization mode"""
        self.show_ray_mode = enabled
        self.ray_incident.setVisible(enabled)
        self.ray_refracted1.setVisible(enabled)
        self.ray_refracted2.setVisible(enabled)
        
        # Initially hide reflection markers (they'll be shown if TIR occurs)
        self.reflection_marker1.setVisible(False)
        self.reflection_marker2.setVisible(False)
        
        if enabled:
            self.update_ray_lines()


    def update_wavelength_scale(self):
        """Update the wavelength scale indicators on the grid"""
        scale_key = (self.frequency, self.n1, self.n2, self.n3)
        
        # Skip update if nothing has changed
        if hasattr(self, '_scale_key') and self._scale_key == scale_key:
            return
        
        # Store new key
        self._scale_key = scale_key

        # Calculate values
        vacuum_wavelength_nm = 299792.458 / self.frequency
        wavelength_m1 = vacuum_wavelength_nm / self.n1
        wavelength_m2 = vacuum_wavelength_nm / self.n2
        wavelength_m3 = vacuum_wavelength_nm / self.n3
        
        speed_m1_formatted = f"{1/self.n1:.2f}c"
        speed_m2_formatted = f"{1/self.n2:.2f}c"
        speed_m3_formatted = f"{1/self.n3:.2f}c"
        
        # Initialize or update labels
        if not hasattr(self, 'grid_value_labels') or len(self.grid_value_labels) == 0:
            # First time creation
            self.grid_value_labels = []
            
            # Create grid line labels
            for x in np.arange(0, 3001, 500):
                physical_nm = x / (500/500)
                label = pg.TextItem(f"{int(physical_nm)} nm", anchor=(0, 0.5), color='white')
                label.setPos(x+10, -1.9)
                self.addItem(label)
                self.grid_value_labels.append(label)
                
            # Create wavelength info labels
            wavelength_info = pg.TextItem("", anchor=(0, 0), color='yellow')
            wavelength_info.setPos(50, 1.5)
            self.addItem(wavelength_info)
            self.grid_value_labels.append(wavelength_info)
            
            # Create medium labels
            for i, pos in enumerate([
                self.boundary1/2,
                self.boundary1 + (self.boundary2-self.boundary1)/2,
                self.boundary2 + (3000-self.boundary2)/2
            ]):
                label = pg.TextItem("", anchor=(0.5, 0), color='yellow')
                label.setPos(pos, 1.5)
                self.addItem(label)
                self.grid_value_labels.append(label)
        
        # Update text content of existing labels
        wavelength_idx = len(self.grid_value_labels) - 4  # Index of wavelength info label
        self.grid_value_labels[wavelength_idx].setText(f"λ₀ = {int(vacuum_wavelength_nm)} nm")
        
        # Update medium labels
        self.grid_value_labels[wavelength_idx+1].setText(f"λ₁ = {int(wavelength_m1)} nm | v₁ = {speed_m1_formatted}")
        self.grid_value_labels[wavelength_idx+2].setText(f"λ₂ = {int(wavelength_m2)} nm | v₂ = {speed_m2_formatted}")
        self.grid_value_labels[wavelength_idx+3].setText(f"λ₃ = {int(wavelength_m3)} nm | v₃ = {speed_m3_formatted}")

    def update_wavelength(self, value):
        """Update the wavelength and colors"""
        self.frequency = value
        self.wavelength = 299792.458 / value


        # Update wave color based on wavelength
        if not self.white_light:
            r, g, b = frequency_to_rgb(value)
            wave_color = QColor(r, g, b)
            
            # Update wave curve color
            self.wave_curve.setPen(pg.mkPen(wave_color, width=4))
            
            # Update ray colors if ray mode is enabled
            if self.show_ray_mode:
                ray_color = pg.mkPen(wave_color, width=4, style=Qt.DashLine)
                self.ray_incident.setPen(ray_color)
                self.ray_refracted1.setPen(ray_color)
                self.ray_refracted2.setPen(ray_color)
        
        # Update the wavelength scale
        self.update_wavelength_scale()
        
        # Update interference waves if they're visible
        if self.show_interference:
            wave1, wave2 = self.calculate_interference_waves(self.wavelength)
            self.interference_wave1.setData(self.x, wave1)
            self.interference_wave2.setData(self.x, wave2)
        
        # Update the wave data - use frequency directly
        wave = self.calculate_wave(self.frequency)
        self.wave_curve.setData(self.x, wave)

        # Update ray visualization if enabled
        if self.show_ray_mode:
            self.update_ray_lines()
            
    def update_amplitude(self, value):
        """Update the amplitude of the wave"""
        self.amplitude = value
        


        # Store current y-range before updating
        y_range = self.getPlotItem().getViewBox().viewRange()[1]
        
        # Update the wave with new amplitude
        wave = self.calculate_wave(self.wavelength)
        self.wave_curve.setData(self.x, wave)
        
        # Restore the previous y-range to prevent auto-zooming
        self.getPlotItem().getViewBox().setYRange(y_range[0], y_range[1], padding=0)

    def update_speed(self, value):
        self.speed = value
        


    def update_n1(self, value):
        """Update refractive index of medium 1"""
        self.n1 = value
        
        # Only update the medium label text, not the entire plot
        for label_pair in self.medium_labels:
            if label_pair[2] == 1:  
                label_pair[0].setText(f"n₁: {self.n1:.4f}")
                break
    

        # Just update the wave data
        if not self.white_light:
            wave = self.calculate_wave(self.frequency)
            self.wave_curve.setData(self.x, wave)
        
        # Update ray visualization if enabled
        if self.show_ray_mode:
            self.update_ray_lines()
        
    def update_n2(self, value):
        """Update refractive index of medium 2"""
        self.n2 = value
        for label_pair in self.medium_labels:
            if label_pair[2] == 2:  # Medium 2
                label_pair[0].setText(f"n₂: {self.n2:.4f}")
                break
        
            # Just update the wave data
        if not self.white_light:
            wave = self.calculate_wave(self.frequency)
            self.wave_curve.setData(self.x, wave)
        
        # Update ray visualization if enabled
        if self.show_ray_mode:
            self.update_ray_lines()            

    def update_n3(self, value):
        """Update refractive index of medium 3"""
        self.n3 = value
# Only update the medium label text, not the entire plot
        for label_pair in self.medium_labels:
            if label_pair[2] == 3:  # Medium 3
                label_pair[0].setText(f"n₃: {self.n3:.4f}")
                break
        
        # Just update the wave data
        if not self.white_light:
            wave = self.calculate_wave(self.frequency)
            self.wave_curve.setData(self.x, wave)
        
        # Update ray visualization if enabled
        if self.show_ray_mode:
            self.update_ray_lines()
        
        
    def toggle_prism_mode(self, enabled):
        """Toggle between normal and prism simulation mode"""
        self.prism_mode = enabled
        self.update_plot()
        
    def toggle_white_light(self, enabled):
        """Toggle white light mode"""
        self.white_light = enabled
        if enabled:
            # Create the curves for visible spectrum
            self.create_white_light_curves()

            # Hide the single wave curve
            if hasattr(self, 'wave_curve'):
                self.wave_curve.setVisible(False)
                
            # Enable superposition by default in white light mode
            self.superposition_enabled = True
            
            # Create superposition wave if it doesn't exist
            if self.superposition_wave is None:
                self.superposition_wave = self.plot(self.x, np.zeros_like(self.x), 
                                                   pen=pg.mkPen('w', width=4))
            
            # Show superposition wave, hide individual components
            self.superposition_wave.setVisible(True)
            for curve_item in self.wave_curves:
                if isinstance(curve_item, tuple) and len(curve_item) == 2:
                    curve, _ = curve_item
                    if hasattr(curve, 'setVisible'):
                        curve.setVisible(False)
                        
            # Update the superposition wave
            self.update_superposition_wave()
            
            # Update the superposition button text if it exists
            if hasattr(self, 'wave_widget') and hasattr(self.wave_widget, 'superposition_button'):
                self.wave_widget.superposition_button.setText("Show Components")
                self.wave_widget.superposition_button.setChecked(False)
        else:
            # Show the single wave curve
            if hasattr(self, 'wave_curve'):
                self.wave_curve.setVisible(True)
                
            # Hide all white light curves
            for curve_item in self.wave_curves:
                if isinstance(curve_item, tuple) and len(curve_item) == 2:
                    curve, _ = curve_item
                    if hasattr(curve, 'setVisible'):
                        curve.setVisible(False)
            
            # Hide superposition wave
            if self.superposition_wave:
                self.superposition_wave.setVisible(False)
                
            # Reset the superposition button text if it exists
            if hasattr(self, 'wave_widget') and hasattr(self.wave_widget, 'superposition_button'):
                self.wave_widget.superposition_button.setText("Superposition")


    def create_white_light_curves(self):
        """Create curves for visible spectrum wavelengths"""
        # Clear any existing curves
        for curve_item in self.wave_curves:
            if isinstance(curve_item, tuple) and len(curve_item) == 2:
                curve, _ = curve_item
                if hasattr(curve, 'setVisible'):
                    self.removeItem(curve)  # Remove the curve from the plot
        
        self.wave_curves = []
        
        # Disable updates while creating curves to improve performance
        self.setUpdatesEnabled(False)
        
        # Define frequency range for visible spectrum (THz)
        min_freq = 400  # ~750 nm (red)
        max_freq = 790  # ~380 nm (violet)
        
        # Adaptive resolution based on system performance
        # Check if we're running smoothly
        if hasattr(self, 'skip_frames') and self.skip_frames > 0:
            # Lower resolution if we're skipping frames
            num_frequencies = 5
            sample_rate = 30  # Take every 30th point
        else:
            # Higher resolution if performance is good
            num_frequencies = 7
            sample_rate = 20  # Take every 20th point
        
        # Use key colors across the spectrum
        frequencies = np.linspace(min_freq, max_freq, num_frequencies)
        
        # Use reduced number of points for white light component waves
        reduced_x = self.x[::sample_rate]
        
        for freq in frequencies:
            # Convert frequency to wavelength
            wavelength_nm = 299792.458 / freq
            
            # Calculate wave with reduced points
            wave = self.calculate_wave(freq)
            reduced_wave = wave[::sample_rate]
            
            # Get color for this wavelength
            r, g, b = wavelength_to_rgb(wavelength_nm)
            wave_color = QColor(r, g, b)
            
            # Create curve with this color and reduced points
            # Use thicker lines to make them more visible despite fewer points
            curve = self.plot(reduced_x, reduced_wave, pen=pg.mkPen(wave_color, width=5))
            curve.setVisible(False)  # Initially hidden
            
            # Store curve and frequency
            self.wave_curves.append((curve, freq))
        
        # Re-enable updates after all curves are created
        self.setUpdatesEnabled(True)

    def update_superposition_wave(self):
        """Calculate and update the superposition wave from all individual waves"""
        if not self.superposition_enabled or not self.white_light or self.superposition_wave is None:
            return
            
        # Create cache key based on current state
        cache_key = (self.time, self.n1, self.n2, self.n3)
        
        # Check if we have this calculation cached
        if hasattr(self, '_superposition_cache') and cache_key in self._superposition_cache:
            superposition = self._superposition_cache[cache_key]
        else:
            # Start with zeros
            superposition = np.zeros_like(self.x)
            
            # Add all individual waves
            for curve, freq in self.wave_curves:
                wave = self.calculate_wave(freq)
                superposition += wave
                
            # Scale the superposition to keep it within reasonable amplitude
            if len(self.wave_curves) > 0:
                superposition = superposition / len(self.wave_curves)
            
            # Cache the result
            if not hasattr(self, '_superposition_cache'):
                self._superposition_cache = {}
            self._superposition_cache[cache_key] = superposition
            
            # Limit cache size
            if len(self._superposition_cache) > 50:
                oldest_keys = list(self._superposition_cache.keys())[:10]
                for key in oldest_keys:
                    self._superposition_cache.pop(key)
            
        # Update the superposition wave
        self.superposition_wave.setData(self.x, superposition)

    def update_ray_lines(self):
        """Update the ray lines to visualize refraction angles"""
        try:
            # Create cache key based on current state
            cache_key = (self.angle_of_incidence, self.n1, self.n2, self.n3, self.ray_target_y)
            
            # Check if we have this calculation cached
            if hasattr(self, '_ray_cache') and cache_key in self._ray_cache:
                ray_data = self._ray_cache[cache_key]
                
                # Apply cached data
                self.ray_incident.setData(ray_data['incident_x'], ray_data['incident_y'])
                self.ray_refracted1.setData(ray_data['refracted1_x'], ray_data['refracted1_y'])
                self.ray_refracted2.setData(ray_data['refracted2_x'], ray_data['refracted2_y'])
                self.ray_incident.setVisible(True)
                self.ray_refracted1.setVisible(ray_data['refracted1_visible'])
                self.ray_refracted2.setVisible(ray_data['refracted2_visible'])
                self.reflection_marker1.setVisible(ray_data['reflection1_visible'])
                self.reflection_marker2.setVisible(ray_data['reflection2_visible'])
                if ray_data['reflection1_visible']:
                    self.reflection_marker1.setData([ray_data['reflection1_x']], [ray_data['reflection1_y']])
                if ray_data['reflection2_visible']:
                    self.reflection_marker2.setData([ray_data['reflection2_x']], [ray_data['reflection2_y']])
                
                # Update angle labels
                self.angle1_label.setHtml(ray_data['angle1_html'])
                self.angle2_label.setHtml(ray_data['angle2_html'])
                self.angle3_label.setHtml(ray_data['angle3_html'])
                
                return
        

            angle_refraction1_deg = 0
            angle_refraction2_deg = 0


            label_html_style = """<div style="font-family: Arial; font-size: 16pt; font-weight: bold; color: white;">"""
            # Hide reflection markers initially
            self.reflection_marker1.setVisible(False)
            self.reflection_marker2.setVisible(False)
            
            # Get the refractive indices of the media
            n1 = self.n1
            n2 = self.n2
            n3 = self.n3
            
            # Get the color from wavelength for all ray segments
            if not self.white_light:
                # Get RGB tuple
                r, g, b = wavelength_to_rgb(self.wavelength)
                # Convert to Qt color
                ray_color = pg.mkPen(QColor(r, g, b), width=4, style=Qt.DashLine)
                self.ray_incident.setPen(ray_color)
                self.ray_refracted1.setPen(ray_color)
                self.ray_refracted2.setPen(ray_color)
            else:
                # Default colors for white light mode
                self.ray_incident.setPen(pg.mkPen('#ffffff', width=4, style=Qt.DashLine))
                self.ray_refracted1.setPen(pg.mkPen('#00aaff', width=4, style=Qt.DashLine))
                self.ray_refracted2.setPen(pg.mkPen('#22ff22', width=4, style=Qt.DashLine))
            
            # Use consistent scale factor for better visualization
            # A smaller scale factor makes the angles less steep visually
            scale_factor = 0.008
            
            # Calculate incident angle in radians
            angle_incidence_rad = np.radians(self.angle_of_incidence)
            self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
            # Set target y-position for the ray to hit at the first boundary
            target_y = self.ray_target_y
            boundary1_y = target_y
            
            # Calculate where to start the ray to hit the boundary
            # Make the incident ray longer by using a much larger value
            # Use the entire plot width for medium 1
            incident_ray_length = self.boundary1  # Use the entire medium 1 width
            start_x = 0  # Start from left edge
            start_y = boundary1_y - np.tan(angle_incidence_rad) * incident_ray_length * scale_factor
            
            # Calculate the critical angles if applicable
            critical_angle1 = None
            critical_angle1_deg = None
            if n1 > n2:
                critical_angle1 = np.arcsin(n2/n1)
                critical_angle1_deg = np.degrees(critical_angle1)
                print(f"Critical angle from medium 1 to 2: {critical_angle1_deg:.2f}°")
            
            critical_angle2 = None
            critical_angle2_deg = None
            if n2 > n3:
                critical_angle2 = np.arcsin(n3/n2)
                critical_angle2_deg = np.degrees(critical_angle2)
                print(f"Critical angle from medium 2 to 3: {critical_angle2_deg:.2f}°")
            
            # Calculate sin(theta2) using Snell's law: n1*sin(theta1) = n2*sin(theta2)
            sin_theta2 = n1 * np.sin(angle_incidence_rad) / n2
            print(f"sin(theta2) = {sin_theta2:.6f}")
            
            # First set the incident ray - this is always displayed
            self.ray_incident.setData([start_x, self.boundary1], [start_y, boundary1_y])
            self.ray_incident.setVisible(True)
            
            # Check for TIR at first boundary
            if abs(sin_theta2) >= 1.0:
                # TIR is occurring - show reflection
                print(f"TIR DETECTED at boundary 1! |sin(theta2)| = {abs(sin_theta2):.6f} > 1.0")
                if critical_angle1_deg:
                    print(f"Incident angle: {self.angle_of_incidence}° > Critical angle: {critical_angle1_deg:.2f}°")
                
                # Update angle labels
                self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
                self.angle2_label.setHtml(f"{label_html_style}θ₂: {angle_refraction1_deg:.1f}°</div>")
                self.angle3_label.setHtml(f"{label_html_style}θ₃: {angle_refraction2_deg:.1f}°</div>")
                
                # Calculate the reflection angle (mirror of incident angle across normal)
                # Law of reflection: angle of incidence = angle of reflection (from normal)
                # In our coordinate system, the reflection goes up and to the left
                # Use same length for reflected ray as incident ray for symmetry
                reflected_end_x = start_x  # Same x-distance from boundary
                reflected_end_y = -start_y  # Mirrored vertically (y-axis is the normal)
                
                # Set the reflected ray
                self.ray_refracted1.setData([self.boundary1, reflected_end_x], [boundary1_y, reflected_end_y])
                self.ray_refracted1.setVisible(True)
                
                # Add reflection marker at the boundary
                self.reflection_marker1.setData([self.boundary1], [boundary1_y])
                self.reflection_marker1.setVisible(True)
                
                # Hide the second refracted ray
                self.ray_refracted2.setVisible(False)
                
                return
            
            # No TIR at first boundary - safe to calculate refraction angle
            angle_refraction1_rad = np.arcsin(sin_theta2)
            angle_refraction1_deg = np.degrees(angle_refraction1_rad)
            print(f"Refraction angle at boundary 1: {angle_refraction1_deg:.2f}°")
            
            # Calculate first refracted ray
            # Determine distance to second boundary
            boundary_distance = self.boundary2 - self.boundary1
            
            # Calculate endpoint for first refracted ray
            refracted1_end_x = self.boundary2
            refracted1_end_y = boundary1_y + np.tan(angle_refraction1_rad) * boundary_distance * scale_factor
            
            # Set first refracted ray
            self.ray_refracted1.setData([self.boundary1, refracted1_end_x], [boundary1_y, refracted1_end_y])
            self.ray_refracted1.setVisible(True)
            
            # Calculate angle of incidence at second boundary (relative to normal)
            # This is the same as angle_refraction1_rad for our coordinate system
            angle_incidence2_rad = angle_refraction1_rad
            angle_incidence2_deg = angle_refraction1_deg
            
            # Check for TIR at second boundary
            # Calculate sin(theta3) using Snell's law
            sin_theta3 = n2 * np.sin(angle_incidence2_rad) / n3
            print(f"sin(theta3) = {sin_theta3:.6f}")
            
            # Check if incident angle at second boundary exceeds critical angle
            # Need to compare absolute values for negative angles
            tir_at_boundary2 = False
            if critical_angle2 is not None:
                # TIR occurs when |angle of incidence| > critical angle
                if abs(angle_incidence2_rad) > critical_angle2:
                    tir_at_boundary2 = True
                    print(f"TIR at boundary 2! |{angle_incidence2_deg:.2f}°| > {np.degrees(critical_angle2):.2f}°")
            
            # Also check if |sin(theta3)| > 1, which is another way to detect TIR
            if abs(sin_theta3) >= 1.0:
                tir_at_boundary2 = True
                print(f"TIR DETECTED at boundary 2! |sin(theta3)| = {abs(sin_theta3):.6f} > 1.0")
            
            if tir_at_boundary2:
                # TIR at second boundary
                # Update angle labels
                self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
                self.angle2_label.setHtml(f"{label_html_style}θ₂: {angle_refraction1_deg:.1f}°</div>")
                self.angle3_label.setHtml(f"{label_html_style}θ₃: TIR</div>")
                
                # Calculate reflection at second boundary
                # Reflection angle equals angle of incidence across normal
                reflection_angle2_rad = -angle_incidence2_rad # Mirror across normal
                
                # Use similar distance for reflected ray
                reflected2_end_x = self.boundary1  # Back to first boundary for clear visualization
                reflected2_delta_x = self.boundary2 - reflected2_end_x
                reflected2_end_y = refracted1_end_y + np.tan(reflection_angle2_rad) * reflected2_delta_x * scale_factor
                
                # Set second ray segment (reflected ray)
                self.ray_refracted2.setData([self.boundary2, reflected2_end_x], [refracted1_end_y, reflected2_end_y])
                self.ray_refracted2.setVisible(True)
                
                # Add reflection marker at second boundary
                self.reflection_marker2.setData([self.boundary2], [refracted1_end_y])
                self.reflection_marker2.setVisible(True)
                
                return
                
            # No TIR at second boundary - calculate normal refraction
            try:
                angle_refraction2_rad = np.arcsin(sin_theta3)
                angle_refraction2_deg = np.degrees(angle_refraction2_rad)
                print(f"Refraction angle at boundary 2: {angle_refraction2_deg:.2f}°")
                
                # Calculate second refracted ray
                # Use a fixed length for visualization clarity
                refracted2_length = 1000
                refracted2_end_x = self.boundary2 + refracted2_length
                refracted2_end_y = refracted1_end_y + np.tan(angle_refraction2_rad) * refracted2_length * scale_factor
                
                # Set second refracted ray
                self.ray_refracted2.setData([self.boundary2, refracted2_end_x], [refracted1_end_y, refracted2_end_y])
                self.ray_refracted2.setVisible(True)
                
                # Update angle labels
                self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
                self.angle2_label.setHtml(f"{label_html_style}θ₂: {angle_refraction1_deg:.1f}°</div>")
                self.angle3_label.setHtml(f"{label_html_style}θ₃: {angle_refraction2_deg:.1f}°</div>")

                if not hasattr(self, '_ray_cache'):
                    self._ray_cache = {}

                # Store all the calculated data
                ray_data = {
                    'incident_x': [start_x, self.boundary1],
                    'incident_y': [start_y, boundary1_y],
                    'refracted1_x': [self.boundary1, refracted1_end_x],
                    'refracted1_y': [boundary1_y, refracted1_end_y],
                    'refracted2_x': [self.boundary2, refracted2_end_x],
                    'refracted2_y': [refracted1_end_y, refracted2_end_y],
                    'refracted1_visible': self.ray_refracted1.isVisible(),
                    'refracted2_visible': self.ray_refracted2.isVisible(),
                    'reflection1_visible': self.reflection_marker1.isVisible(),
                    'reflection2_visible': self.reflection_marker2.isVisible(),
                    'reflection1_x': self.boundary1,
                    'reflection1_y': boundary1_y,
                    'reflection2_x': self.boundary2,
                    'reflection2_y': refracted1_end_y,
                    'angle1_html': self.angle1_label.toHtml(),
                    'angle2_html': self.angle2_label.toHtml(),
                    'angle3_html': self.angle3_label.toHtml()
                }
                
                self._ray_cache[cache_key] = ray_data
                
                # Limit cache size
                if len(self._ray_cache) > 100:
                    oldest_keys = list(self._ray_cache.keys())[:20]
                    for key in oldest_keys:
                        self._ray_cache.pop(key)

            except Exception as e:
                # If arcsin calculation fails, it should be TIR - handle accordingly
                print(f"Exception in calculating refraction angle: {e}")
                
                # Mark as TIR
                self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
                self.angle2_label.setHtml(f"{label_html_style}θ₂: {angle_refraction1_deg:.1f}°</div>")
                self.angle3_label.setHtml(f"{label_html_style}θ₃: TIR (error)</div>")
                
                # Calculate reflection instead
                reflection_angle2_rad = -angle_incidence2_rad
                
                # Calculate reflected ray
                reflected2_end_x = self.boundary1
                reflected2_delta_x = self.boundary2 - reflected2_end_x
                reflected2_end_y = refracted1_end_y + np.tan(reflection_angle2_rad) * reflected2_delta_x * scale_factor
                
                # Set second ray segment (reflected ray)
                self.ray_refracted2.setData([self.boundary2, reflected2_end_x], [refracted1_end_y, reflected2_end_y])
                self.ray_refracted2.setVisible(True)
                
                # Add reflection marker at second boundary
                self.reflection_marker2.setData([self.boundary2], [refracted1_end_y])
                self.reflection_marker2.setVisible(True)
                
        except Exception as e:
            # Top-level exception handling to prevent crashes
            print(f"Error in ray calculation: {type(e).__name__}: {e}")
            # Provide a simple ray visualization
            self.ray_incident.setData([0, self.boundary1], [0, 0])
            self.ray_refracted1.setData([self.boundary1, self.boundary2], [0, 0])
            self.ray_refracted2.setData([self.boundary2, 3000], [0, 0])
            self.ray_incident.setVisible(True)
            self.ray_refracted1.setVisible(True)
            self.ray_refracted2.setVisible(True)
            # Set error labels
            self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
            self.angle2_label.setHtml(f"{label_html_style}θ₂: Error</div>")
            self.angle3_label.setHtml(f"{label_html_style}θ₃: Error</div>")

    def update_angle(self, value):
        """Update the angle of incidence"""
        self.angle_of_incidence = value
        
        # Use cached ray calculations if available
        cache_key = (self.angle_of_incidence, self.n1, self.n2, self.n3, self.ray_target_y)
        if hasattr(self, '_ray_cache') and cache_key in self._ray_cache:
            ray_data = self._ray_cache[cache_key]
            self.angle1_label.setHtml(ray_data['angle1_html'])
            self.angle2_label.setHtml(ray_data['angle2_html'])
            self.angle3_label.setHtml(ray_data['angle3_html'])
        else:

            # Update angle labels even when ray mode is not enabled
            label_html_style = """<div style="font-family: Arial; font-size: 16pt; font-weight: bold; color: white;">"""
            self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
            
            # Calculate refraction angles using Snell's law
            try:
                # First boundary: n1*sin(θ1) = n2*sin(θ2)
                sin_theta2 = self.n1 * np.sin(np.radians(self.angle_of_incidence)) / self.n2
                
                if abs(sin_theta2) <= 1.0:  # Check for total internal reflection
                    angle_refraction1_deg = np.degrees(np.arcsin(sin_theta2))
                    self.angle2_label.setHtml(f"{label_html_style}θ₂: {angle_refraction1_deg:.1f}°</div>")
                    
                    # Second boundary: n2*sin(θ2) = n3*sin(θ3)
                    sin_theta3 = self.n2 * np.sin(np.radians(angle_refraction1_deg)) / self.n3
                    
                    if abs(sin_theta3) <= 1.0:  # Check for total internal reflection
                        angle_refraction2_deg = np.degrees(np.arcsin(sin_theta3))
                        self.angle3_label.setHtml(f"{label_html_style}θ₃: {angle_refraction2_deg:.1f}°</div>")
                    else:
                        self.angle3_label.setHtml(f"{label_html_style}θ₃: TIR</div>")
                else:
                    self.angle2_label.setHtml(f"{label_html_style}θ₂: TIR</div>")
                    self.angle3_label.setHtml(f"{label_html_style}θ₃: N/A</div>")
            except Exception as e:
                print(f"Error updating angle labels: {e}")
            
            # Update ray lines if ray mode is enabled
            if self.show_ray_mode:
                self.update_ray_lines()

    def set_ray_target_y(self, y_value):
        """Set the Y-coordinate where the ray hits the first boundary"""
        self.ray_target_y = y_value
        if self.show_ray_mode:
            self.update_ray_lines()

    def cleanup(self):
        """Clean up resources to prevent memory leaks"""
        if hasattr(self, 'timer') and self.timer.isActive():
            self.timer.stop()
        
        # Clear all cached data
        if hasattr(self, '_wave_calc_cache'):
            self._wave_calc_cache = {}
        if hasattr(self, '_superposition_cache'):
            self._superposition_cache = {}
        if hasattr(self, '_interference_cache'):
            self._interference_cache = {}
        if hasattr(self, '_ray_cache'):
            self._ray_cache = {}

class ModernCheckBox(QCheckBox):
    def __init__(self, text="", parent=None):
        super().__init__(text, parent)
        self.setStyleSheet("""
            QCheckBox {
                spacing: 8px;
                font-size: 14px;
            }
            
            QCheckBox::indicator {
                width: 20px;
                height: 20px;
                border-radius: 10px;
                border: 2px solid #555555;
                background-color: #333337;
            }
            
            QCheckBox::indicator:unchecked:hover {
                border: 2px solid #007ACC;
            }
            
            QCheckBox::indicator:checked {
                background-color: #007ACC;
                border: 2px solid #007ACC;
                image: url(data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNCIgaGVpZ2h0PSIxNCIgdmlld0JveD0iMCAwIDI0IDI0IiBmaWxsPSJub25lIiBzdHJva2U9IiNmZmZmZmYiIHN0cm9rZS13aWR0aD0iMyIgc3Ryb2tlLWxpbmVjYXA9InJvdW5kIiBzdHJva2UtbGluZWpvaW49InJvdW5kIj48cG9seWxpbmUgcG9pbnRzPSIyMCA2IDkgMTcgNCAxMiI+PC9wb2x5bGluZT48L3N2Zz4=);
            }
            
            QCheckBox::indicator:checked:hover {
                background-color: #1C97EA;
                border: 2px solid #1C97EA;
            }
        """)

class ColoredSlider(QSlider):
    def __init__(self, parent=None):
        super().__init__(Qt.Horizontal, parent)
        self.setStyleSheet("""
            QSlider::groove:horizontal {
                border: none;
                height: 10px;
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0.000 #ff0000,    /* 400 THz - red */
                    stop:0.200 #ff8800,    /* 460 THz - orange */
                    stop:0.350 #ffff00,    /* 510 THz - yellow */
                    stop:0.500 #00ff00,    /* 540 THz - green */
                    stop:0.650 #0088ff,    /* 580 THz - light blue */
                    stop:0.800 #4400ff,    /* 680 THz - deep blue */
                    stop:1.000 #8800ff);   /* 790 THz - violet */
                margin: 2px 0;
                border-radius: 5px;
            }
            QSlider::handle:horizontal {
                background: white;
                border: none;
                width: 18px;
                margin: -4px 0;
                border-radius: 9px;
            }
        """)

class BlueSlider(QSlider):
    def __init__(self, parent=None):
        super().__init__(Qt.Horizontal, parent)
        self.setStyleSheet("""
            QSlider::groove:horizontal {
                border: none;
                height: 10px;
                background: #333337;
                margin: 2px 0;
                border-radius: 5px;
            }
            QSlider::sub-page:horizontal {
                background: #0088ff;
                border: none;
                height: 10px;
                margin: 2px 0;
                border-radius: 5px;
            }
            QSlider::handle:horizontal {
                background: white;
                border: none;
                width: 18px;
                margin: -4px 0;
                border-radius: 9px;
            }
        """)

class LightSimulationApp(QMainWindow):
    def __init__(self):
        super().__init__()
        
        # Set up the main window
        self.setWindowTitle("Light Wave Simulation")

        central_widget = QWidget()
        central_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setCentralWidget(central_widget)
        self.main_layout = QVBoxLayout(central_widget)
        self.main_layout.setContentsMargins(0, 0, 0, 0)  # Remove margins
        self.main_layout.setSpacing(0)

                # Medium presets (refractive indices at ~550nm wavelength)
        self.medium_presets = {
            'Air': {'n': 1.0003, 'color': '#22222259'},  # Subtle gray with 35% opacity
            'Water': {'n': 1.33, 'color': '#1E90FF59'},  # Soft blue with 35% opacity
            'Glass (Crown)': {'n': 1.52, 'color': '#87CEEB59'},  # Sky blue with 35% opacity
            'Glass (Flint)': {'n': 1.62, 'color': '#9370DB59'},  # Medium purple with 35% opacity
            'Diamond': {'n': 2.42, 'color': '#E0E0E059'},  # Light silver with 35% opacity
            'Acrylic': {'n': 1.49, 'color': '#98FB9859'},  # Pale green with 35% opacity
            'Glycerine': {'n': 1.47, 'color': '#F0E68C59'},  # Khaki with 35% opacity
            'Ethanol': {'n': 1.36, 'color': '#D8BFD859'},  # Thistle with 35% opacity
            'Quartz': {'n': 1.54, 'color': '#DEB88759'},  # Burlywood with 35% opacity
            'Sapphire': {'n': 1.77, 'color': '#4682B459'}  # Steel blue with 35% opacity
        }

        
        # Preset scenarios
        self.scenario_materials = {
            'Air → Water → Glass': ('Air', 'Water', 'Glass (Crown)'),
            'Air → Glass → Water': ('Air', 'Glass (Crown)', 'Water'),
            'Water → Air → Glass': ('Water', 'Air', 'Glass (Crown)'),
            'Air → Diamond → Glass': ('Air', 'Diamond', 'Glass (Crown)'),
            'Glass → Air → Water': ('Glass (Crown)', 'Air', 'Water')
        }
        
        
        # Set dark mode stylesheet
        self.setStyleSheet("""
            QMainWindow, QWidget {
                background-color: #2D2D30;
                color: #FFFFFF;
            }
            QGroupBox {
                border: 1px solid #3F3F46;
                border-radius: 5px;
                margin-top: 10px;
                font-weight: bold;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
            QPushButton {
                background-color: #007ACC;
                border: none;
                border-radius: 3px;
                padding: 5px;
                color: white;
            }
            QPushButton:hover {
                background-color: #1C97EA;
            }
            QComboBox, QLineEdit, QSpinBox, QDoubleSpinBox {
                background-color: #333337;
                border: 1px solid #3F3F46;
                border-radius: 3px;
                padding: 2px;
            }
            QSlider::groove:horizontal {
                border: 1px solid #999999;
                height: 8px;
                background: #333337;
                margin: 2px 0;
                border-radius: 2px;
            }
            QSlider::handle:horizontal {
                background: #007ACC;
                border: 1px solid #5c5c5c;
                width: 18px;
                margin: -2px 0;
                border-radius: 3px;
            }
        """)
        
  
 
        
        # Create wave widget directly (no tabs)
        self.wave_widget = WaveSimulationWidget()
        self.wave_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.main_layout.addWidget(self.wave_widget, 1)
        
        # Create a modern play/pause button
        self.play_pause_button = PlayPauseButton()
        self.play_pause_button.clicked.connect(self.toggle_animation)

        # Add the button to a centered container at the bottom of the plot
        button_container = QWidget()
        button_layout = QHBoxLayout(button_container)
        button_layout.addStretch(1)
        button_layout.addWidget(self.play_pause_button)
        button_layout.addStretch(1)
        button_layout.setContentsMargins(0, 5, 0, 5)

  
        self.main_layout.addWidget(button_container)
        
 
        self.wave_widget.setContextMenuPolicy(Qt.NoContextMenu)
        # Share medium presets with the wave widget
        self.wave_widget.medium_presets = self.medium_presets
        
        # Create controls container
        controls_container = QWidget()
        controls_container.setMaximumHeight(200)  # Limit height
        controls_layout = QHBoxLayout()
        controls_layout.setContentsMargins(10, 10, 10, 10)  # Add some padding
        controls_layout.setSpacing(10)  # Add spacing between controls
        controls_container.setLayout(controls_layout)
        self.main_layout.addWidget(controls_container)
        
        # Add wave controls
        self.setup_wave_controls(controls_layout)


        self.wave_widget.medium1_color = self.medium_presets['Air']['color']
        self.wave_widget.medium2_color = self.medium_presets['Water']['color']
        self.wave_widget.medium3_color = self.medium_presets['Glass (Crown)']['color']

        
        self.showMaximized()

    def resizeEvent(self, event):
        """Handle window resize events"""
        super().resizeEvent(event)
        # Reposition the pause button at the bottom left of the wave widget
        if hasattr(self, 'pause_container') and hasattr(self, 'wave_widget'):
            self.pause_container.move(20, self.wave_widget.height() - 100)

    def toggle_animation(self):
        paused = self.play_pause_button.isChecked()
        self.wave_widget.toggle_pause(paused)
        self.play_pause_button.update_icon(paused)
     
    def setup_wave_controls(self, layout):
        # Left controls group (wave properties)
        # Create a group for wave controls
        wave_group = QGroupBox("Wave Controls")
        wave_layout = QVBoxLayout()
        wave_group.setLayout(wave_layout)
        layout.addWidget(wave_group)
        
        
        
        # Wavelength slider with color gradient
        wavelength_layout = QHBoxLayout()
        frequency_label = QLabel("Frequency (THz):")
        self.frequency_slider = ColoredSlider()
        self.frequency_slider.setMinimum(400)
        self.frequency_slider.setMaximum(790)
        self.frequency_slider.setValue(545)  # ~550 nm
        self.frequency_slider.setTickPosition(QSlider.TicksBelow)
        self.frequency_slider.setTickInterval(50)
        
        
        self.wavelength_value = QLabel("550 nm")
        self.frequency_value_label = QLabel("545 THz")

        wavelength_layout.addWidget(frequency_label)
        wavelength_layout.addWidget(self.frequency_slider)  # Use frequency_slider instead of wavelength_slider
        wavelength_layout.addWidget(self.frequency_value_label)
        wave_layout.addLayout(wavelength_layout)
        
        # Amplitude slider
        amplitude_layout = QHBoxLayout()
        amplitude_label = QLabel("Amplitude:")
        self.amplitude_slider = BlueSlider()
        self.amplitude_slider.setMinimum(1)
        self.amplitude_slider.setMaximum(10)
        self.amplitude_slider.setValue(5)
        self.amplitude_value = QLabel("5")
        
        amplitude_layout.addWidget(amplitude_label)
        amplitude_layout.addWidget(self.amplitude_slider)
        amplitude_layout.addWidget(self.amplitude_value)
        wave_layout.addLayout(amplitude_layout)
        
        # Speed slider
        speed_layout = QHBoxLayout()
        speed_label = QLabel("Wave Speed:")
        self.speed_slider = BlueSlider()
        self.speed_slider.setMinimum(1)
        self.speed_slider.setMaximum(10)
        self.speed_slider.setValue(4)
        self.speed_value = QLabel("4")
        
        speed_layout.addWidget(speed_label)
        speed_layout.addWidget(self.speed_slider)
        speed_layout.addWidget(self.speed_value)
        wave_layout.addLayout(speed_layout)

        # Angle of incidence slider
        angle_layout = QHBoxLayout()
        angle_label = QLabel("Angle of Incidence:")
        self.angle_slider = BlueSlider()
        self.angle_slider.setMinimum(-90)
        self.angle_slider.setMaximum(90)
        self.angle_slider.setValue(0)
        self.angle_value = QLabel("0°")
        
        angle_layout.addWidget(angle_label)
        angle_layout.addWidget(self.angle_slider)
        angle_layout.addWidget(self.angle_value)
        wave_layout.addLayout(angle_layout)

        # Ray target Y position slider
        ray_target_layout = QHBoxLayout()
        ray_target_label = QLabel("Ray Target Y:")
        self.ray_target_slider = BlueSlider()
        self.ray_target_slider.setMinimum(-100)
        self.ray_target_slider.setMaximum(100)
        self.ray_target_slider.setValue(0)
        self.ray_target_value = QLabel("0.0")
        
        ray_target_layout.addWidget(ray_target_label)
        ray_target_layout.addWidget(self.ray_target_slider)
        ray_target_layout.addWidget(self.ray_target_value)
        wave_layout.addLayout(ray_target_layout)

        mode_layout = QHBoxLayout()
        self.white_light_check = ModernCheckBox("White Light")
        self.interference_check = ModernCheckBox("Show Interference")
        self.interference_check.setChecked(False)
        self.ray_mode_check = ModernCheckBox("Show Ray Path")
        self.superposition_check = ModernCheckBox("Show Components")
        self.superposition_check.setToolTip("Toggle between showing individual color components or a single superposition wave")


        mode_layout.addWidget(self.white_light_check)
        mode_layout.addWidget(self.interference_check)
        mode_layout.addWidget(self.ray_mode_check)
        mode_layout.addWidget(self.superposition_check)
        wave_layout.addLayout(mode_layout)

        self.superposition_check.setEnabled(False)
        self.white_light_check.stateChanged.connect(self.update_superposition_enabled)
        self.superposition_check.toggled.connect(self.toggle_superposition)

        wave_group.setLayout(wave_layout)
        layout.addWidget(wave_group)
        
        # Middle group (medium 1)
        medium1_group = QGroupBox("Medium 1")
        medium1_layout = QVBoxLayout()
        medium1_group.setLayout(medium1_layout)
        layout.addWidget(medium1_group)
        
        # Medium 1 selection
        self.medium1_combo = QComboBox()
        for medium in sorted(self.medium_presets.keys()):
            self.medium1_combo.addItem(medium)
        medium1_layout.addWidget(self.medium1_combo)
        

        # Medium 1 n slider
        n1_layout = QHBoxLayout()
        n1_label = QLabel("n₁:")
        self.n1_slider = BlueSlider()
        self.n1_slider.setMinimum(100)
        self.n1_slider.setMaximum(300)
        self.n1_slider.setValue(100)
        self.n1_value = QLabel("1.0003")
        
        n1_layout.addWidget(n1_label)
        n1_layout.addWidget(self.n1_slider)
        n1_layout.addWidget(self.n1_value)
        medium1_layout.addLayout(n1_layout)
        self.medium1_combo.setCurrentText(self.wave_widget.medium1_name)

        # Middle group (medium 2)
        medium2_group = QGroupBox("Medium 2")
        medium2_layout = QVBoxLayout()
        medium2_group.setLayout(medium2_layout)
        layout.addWidget(medium2_group)
        
        # Medium 2 selection
        self.medium2_combo = QComboBox()
        for medium in sorted(self.medium_presets.keys()):
            self.medium2_combo.addItem(medium)
        medium2_layout.addWidget(self.medium2_combo)
        

        # Medium 2 n slider
        n2_layout = QHBoxLayout()
        n2_label = QLabel("n₂:")
        self.n2_slider = BlueSlider()
        self.n2_slider.setMinimum(100)
        self.n2_slider.setMaximum(300)
        self.n2_slider.setValue(133)
        self.n2_value = QLabel("1.33")
        
        n2_layout.addWidget(n2_label)
        n2_layout.addWidget(self.n2_slider)
        n2_layout.addWidget(self.n2_value)
        medium2_layout.addLayout(n2_layout)
        self.medium2_combo.setCurrentText(self.wave_widget.medium2_name)


        # Right group (medium 3)
        medium3_group = QGroupBox("Medium 3")
        medium3_layout = QVBoxLayout()
        medium3_group.setLayout(medium3_layout)
        layout.addWidget(medium3_group)
        

        # Medium 3 selection
        self.medium3_combo = QComboBox()
        for medium in sorted(self.medium_presets.keys()):
            self.medium3_combo.addItem(medium)
        medium3_layout.addWidget(self.medium3_combo)
        
        # Medium 3 n slider
        n3_layout = QHBoxLayout()
        n3_label = QLabel("n₃:")
        self.n3_slider = BlueSlider()
        self.n3_slider.setMinimum(100)
        self.n3_slider.setMaximum(300)
        self.n3_slider.setValue(152)
        self.n3_value = QLabel("1.52")
        
        n3_layout.addWidget(n3_label)
        n3_layout.addWidget(self.n3_slider)
        n3_layout.addWidget(self.n3_value)
        medium3_layout.addLayout(n3_layout)
        self.medium3_combo.setCurrentText(self.wave_widget.medium3_name)
        # Scenario presets
        scenario_group = QGroupBox("Presets")
        scenario_layout = QVBoxLayout()
        scenario_group.setLayout(scenario_layout)
        layout.addWidget(scenario_group)
        
        self.scenario_combo = QComboBox()
        for scenario in self.scenario_materials.keys():
            self.scenario_combo.addItem(scenario)
        scenario_layout.addWidget(self.scenario_combo)
        
        apply_button = QPushButton("Apply Preset")
        scenario_layout.addWidget(apply_button)
        
        # Connect signals
        self.frequency_slider.valueChanged.connect(self.update_frequency) 
        self.amplitude_slider.valueChanged.connect(self.update_amplitude)
        self.speed_slider.valueChanged.connect(self.update_speed)
        self.angle_slider.valueChanged.connect(self.update_angle)
        self.ray_target_slider.valueChanged.connect(self.update_ray_target)

        self.n1_slider.valueChanged.connect(self.update_n1)
        self.n2_slider.valueChanged.connect(self.update_n2)
        self.n3_slider.valueChanged.connect(self.update_n3)
        
        self.medium1_combo.currentTextChanged.connect(self.update_medium1)
        self.medium2_combo.currentTextChanged.connect(self.update_medium2)
        self.medium3_combo.currentTextChanged.connect(self.update_medium3)
        
        self.white_light_check.stateChanged.connect(self.toggle_white_light)
        
        apply_button.clicked.connect(self.apply_scenario)

        self.interference_check.stateChanged.connect(self.toggle_interference)
        self.ray_mode_check.stateChanged.connect(self.toggle_ray_mode)
        self.superposition_check.stateChanged.connect(self.toggle_superposition)

    def toggle_interference(self, state):
        """Toggle interference visualization"""
        enabled = state == Qt.Checked
        
        # Toggle interference in the wave widget
        self.wave_widget.toggle_interference(enabled)
        
        # Force an immediate update to ensure proper display
        if enabled:
            # Update the wave widget immediately
            self.wave_widget.update_animation()
        
    def update_frequency(self, value):
        """Update the frequency value and pass to wave widget"""
        self.frequency_value_label.setText(f"{value} THz")
        # Update displayed value
        wavelength_nm = 299792.458 / value
        self.wavelength_value.setText(f"{wavelength_nm:.1f} nm")
        
        # Update widget
        self.wave_widget.frequency = value
        self.wave_widget.wavelength = wavelength_nm

        r, g, b = frequency_to_rgb(value)
        wave_color = QColor(r, g, b)
        # We're not updating the color indicator anymore since it's hidden

        self.wave_widget.update_wavelength(value)
        
    def update_amplitude(self, value):
        """Update amplitude value"""
        self.amplitude_value.setText(str(value))
        self.wave_widget.update_amplitude(value)
        
    def update_speed(self, value):
        """Update speed value"""
        self.speed_value.setText(str(value))
        self.wave_widget.update_speed(value)
        
    def update_n1(self, value):
        """Update refractive index of medium 1"""
        n = value / 100
        self.n1_value.setText(f"{n:.4f}")
        self.wave_widget.update_n1(n)
  
        
    def update_n2(self, value):
        """Update refractive index of medium 2"""
        n = value / 100
        self.n2_value.setText(f"{n:.4f}")
        self.wave_widget.update_n2(n)

        
    def update_n3(self, value):
        """Update refractive index of medium 3"""
        n = value / 100
        self.n3_value.setText(f"{n:.4f}")
        self.wave_widget.update_n3(n)


 
        
    def update_medium1(self, medium_name):
        """Update medium 1 selection"""
        if medium_name in self.medium_presets:
            # Update the refractive index
            n1 = self.medium_presets[medium_name]['n']
            self.n1_slider.setValue(int(n1 * 100))
            
            # Update the wave widget's n1 value
            self.wave_widget.update_n1(n1)
            
            # Set the medium name and color in the wave widget
            self.wave_widget.medium1_name = medium_name
            self.wave_widget.medium1_color = self.medium_presets[medium_name]['color']
            
            # Update the wave widget's medium color
            self.wave_widget.update_plot()
            
    def update_medium2(self, medium_name):
        """Update medium 2 selection"""
        if medium_name in self.medium_presets:
            # Update the refractive index
            n2 = self.medium_presets[medium_name]['n']
            self.n2_slider.setValue(int(n2 * 100))
            
            # Update the wave widget's n2 value
            self.wave_widget.update_n2(n2)
            
            # Set the medium name and color in the wave widget
            self.wave_widget.medium2_name = medium_name
            self.wave_widget.medium2_color = self.medium_presets[medium_name]['color']
            
            # Update the wave widget's medium color
            self.wave_widget.update_plot()
            
    def update_medium3(self, medium_name):
        """Update medium 3 selection"""
        if medium_name in self.medium_presets:
            # Update the refractive index
            n3 = self.medium_presets[medium_name]['n']
            self.n3_slider.setValue(int(n3 * 100))
            
            # Update the wave widget's n3 value
            self.wave_widget.update_n3(n3)
            
            # Set the medium name and color in the wave widget
            self.wave_widget.medium3_name = medium_name
            self.wave_widget.medium3_color = self.medium_presets[medium_name]['color']
            
            # Update the wave widget's medium color
            self.wave_widget.update_plot()
        
    def toggle_white_light(self, state):
        """Toggle white light mode"""
        enabled = state == Qt.Checked
        self.wave_widget.toggle_white_light(enabled)

    def update_superposition_enabled(self, state):
        """Enable superposition checkbox only when white light is enabled"""
        self.superposition_check.setEnabled(state == Qt.Checked)
        if state != Qt.Checked:
            self.superposition_check.setChecked(False)
            
    def toggle_superposition(self, state):
        """Toggle superposition visualization in white light mode"""
        enabled = state == Qt.Checked
        self.wave_widget.toggle_superposition(enabled)    
        
    def apply_scenario(self):
        """Apply the selected scenario preset"""
        scenario_name = self.scenario_combo.currentText()
        if scenario_name in self.scenario_materials and hasattr(self, 'medium_presets'):
            medium1, medium2, medium3 = self.scenario_materials[scenario_name]
            
            # Update medium selections
            self.medium1_combo.setCurrentText(medium1)
            self.medium2_combo.setCurrentText(medium2)
            self.medium3_combo.setCurrentText(medium3)
                        
            
            # Force a plot update
            self.wave_widget.update_plot()

            # Update refractive indices
            if medium1 in self.medium_presets:
                n1 = self.medium_presets[medium1]['n']
                self.n1_slider.setValue(int(n1 * 100))
                self.n1_value.setText(f"{n1:.4f}")
            
            if medium2 in self.medium_presets:
                n2 = self.medium_presets[medium2]['n']
                self.n2_slider.setValue(int(n2 * 100))
                self.n2_value.setText(f"{n2:.4f}")
            
            if medium3 in self.medium_presets:
                n3 = self.medium_presets[medium3]['n']
                self.n3_slider.setValue(int(n3 * 100))
                self.n3_value.setText(f"{n3:.4f}")
        
    def update_angle(self, value):
        """Update the angle value and pass to wave widget"""
        self.angle_value.setText(f"{value}°")
        
        # Update the wave widget's angle
        self.wave_widget.update_angle(value)
        
        # Force an update of the wave visualization
        self.wave_widget.update_plot()
        # Update ray visualization
        if self.wave_widget.show_ray_mode:
            self.wave_widget.update_ray_lines()
            
    def update_ray_target(self, value):
        """Update the ray target Y position"""
        # Scale the slider value to a smaller range for better control
        y_value = value / 100.0
        # Update label
        self.ray_target_value.setText(f"{y_value:.1f}")
        # Update target position in wave widget
        self.wave_widget.set_ray_target_y(y_value)
            
    def toggle_ray_mode(self, state):
        """Toggle ray mode"""
        enabled = state == Qt.Checked
        self.wave_widget.toggle_ray_mode(enabled)

    def closeEvent(self, event):
        """Handle window close event to ensure proper cleanup"""
        # Stop all timers
        if hasattr(self.wave_widget, 'timer') and self.wave_widget.timer.isActive():
            self.wave_widget.timer.stop()
        
        # Clear all plot items to release memory
        if hasattr(self.wave_widget, 'clear'):
            self.wave_widget.clear()
        
        # Delete wave curves explicitly
        if hasattr(self.wave_widget, 'wave_curves'):
            for curve_item in self.wave_widget.wave_curves:
                if isinstance(curve_item, tuple) and len(curve_item) == 2:
                    curve, _ = curve_item
                    if hasattr(self.wave_widget, 'removeItem'):
                        self.wave_widget.removeItem(curve)
        
        # Accept the close event
        event.accept()
        
        # Force garbage collection
        import gc
        gc.collect()

# Run the application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = LightSimulationApp()
    window.show()
    app.exec_()
    # Ensure cleanup happens when application exits
    if hasattr(window, 'wave_widget') and hasattr(window.wave_widget, 'cleanup'):
        window.wave_widget.cleanup()