import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                            QLabel, QSlider, QPushButton, QComboBox, QGroupBox, QFrame,
                            QCheckBox, QTabWidget, QRadioButton, QButtonGroup)
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
        self.medium1_color = '#3A3A3A80'  
        self.medium2_color = '#0D47A1C0'  
        self.medium3_color = '#4488AACC'  
        
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
        
        # Create and add the medium rectangles with proper fill
        # Medium 1 (Air)
        medium1_rect = pg.FillBetweenItem(
            pg.PlotCurveItem([0, self.boundary1], [2, 2]),
            pg.PlotCurveItem([0, self.boundary1], [-2, -2]),
            brush=pg.mkBrush(self.medium1_color))
        self.addItem(medium1_rect)
        self.medium_rects.append(medium1_rect)

        
        # Medium 2 (Water)
        medium2_rect = pg.FillBetweenItem(
            pg.PlotCurveItem([self.boundary1, self.boundary2], [2, 2]),
            pg.PlotCurveItem([self.boundary1, self.boundary2], [-2, -2]),
            brush=pg.mkBrush(self.medium2_color))
        self.addItem(medium2_rect)
        self.medium_rects.append(medium2_rect)
        
        # Medium 3 (Glass)
        medium3_rect = pg.FillBetweenItem(
            pg.PlotCurveItem([self.boundary2, 3000], [2, 2]),
            pg.PlotCurveItem([self.boundary2, 3000], [-2, -2]),
            brush=pg.mkBrush(self.medium3_color))
        self.addItem(medium3_rect)
        self.medium_rects.append(medium3_rect)
        
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


        # Set up the animation timer
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_animation)
        self.timer.start(16)  # 50ms interval (20 fps)
        
    def calculate_interference_waves(self, wavelength):
        """Calculate the component waves that create interference"""
        if not self.show_interference:
            return np.zeros_like(self.x), np.zeros_like(self.x)
            
        cache_key = (
            wavelength,
            self.time,
            self.frequency,
            self.n1,
            self.n2,
            self.n3,
            self.amplitude,
            self.speed
        )

        if hasattr(self, '_interference_cache') and self._interference_cache.get('key') == cache_key:
            return self._interference_cache.get('wave1'), self._interference_cache.get('wave2')


        # Get the actual refracted wave
        actual_wave = self.calculate_wave(self.frequency)
        
        # Wave 1 (red) - original wave as if it was coming from vacuum (n=1.0)
        wave1 = np.zeros_like(self.x)
        
        # Calculate vacuum wavelength
        vacuum_wavelength = wavelength
        # Calculate wave number for vacuum (n=1.0)
        k_vacuum = (2 * np.pi * 1.0) / vacuum_wavelength
        
        # Calculate phase with time component to ensure animation
        phase_vacuum = k_vacuum * self.x - self.speed * self.time
        
        # For all media, calculate the wave as if it was in vacuum (n=1.0)
        # This represents the original wave without any medium effects
        wave1 = self.amplitude * self.visualization_scale * np.sin(phase_vacuum)
        
        # Wave 2 (green) - the difference between actual wave and vacuum wave
        # This represents the effect of the media on the wave
        wave2 = actual_wave - wave1
        
        if not hasattr(self, '_interference_cache'):
            self._interference_cache = {}
        self._interference_cache['key'] = cache_key
        self._interference_cache['wave1'] = wave1
        self._interference_cache['wave2'] = wave2

        return wave1, wave2
        
    def calculate_wave(self, frequency):
        """Calculate the wave based on current parameters for a specific wavelength"""

        cache_key = (
            frequency, 
            self.time,
            self.n1, 
            self.n2, 
            self.n3, 
            self.angle_of_incidence,
            self.boundary1, 
            self.boundary2,
            self.prism_mode
        )

        if hasattr(self, '_wave_cache') and self._wave_cache.get('key') == cache_key:
            return self._wave_cache.get('wave')

        wave = np.zeros_like(self.x)
        wavelength = 299792.458 / frequency
        
        if self.prism_mode:
            wl_nm = wavelength
            
            # Apply Cauchy's equation for dispersion: n(λ) = A + B/λ² + C/λ⁴
            # Using simplified coefficients for glass-like materials
            n2_wl = self.n2 + 0.0006 * (550 / wl_nm) ** 2
            n3_wl = self.n3 + 0.0008 * (550 / wl_nm) ** 2
        else:
            # No dispersion in normal mode
            n2_wl = self.n2
            n3_wl = self.n3
        
        # Angle of incidence in radians
        angle_incidence = np.radians(self.angle_of_incidence)
        
        # Add safety checks for total internal reflection
        sin_refraction1 = np.sin(angle_incidence) * self.n1 / n2_wl
        if abs(sin_refraction1) > 1:
            # Total internal reflection at first boundary
            angle_refraction1 = np.pi/2  # Set to 90 degrees
            angle_refraction2 = np.pi/2  # Set to 90 degrees
            total_reflection = True
        else:
            angle_refraction1 = np.arcsin(sin_refraction1)
            sin_refraction2 = np.sin(angle_refraction1) * n2_wl / n3_wl
            if abs(sin_refraction2) > 1:
                # Total internal reflection at second boundary
                angle_refraction2 = np.pi/2  # Set to 90 degrees
                total_reflection = True
            else:
                angle_refraction2 = np.arcsin(sin_refraction2)
                total_reflection = False
                
        # Update angle labels if they exist
        if hasattr(self, 'angle1_label') and hasattr(self, 'angle2_label') and hasattr(self, 'angle3_label'):
            label_html_style = """<div style="font-family: Arial; font-size: 16pt; font-weight: bold; color: white;">"""
            self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
            self.angle2_label.setHtml(f"{label_html_style}θ₂: {np.degrees(angle_refraction1):.1f}°</div>")
            self.angle3_label.setHtml(f"{label_html_style}θ₃: {np.degrees(angle_refraction2):.1f}°</div>")

            
        # Wave parameters for each medium
        k1 = (2 * np.pi * self.n1) / wavelength  # Wave number in medium 1
        k2 = (2 * np.pi * n2_wl) / wavelength    # Wave number in medium 2
        k3 = (2 * np.pi * n3_wl) / wavelength    # Wave number in medium 3
        
        # Create y-coordinates for vertical displacement
        y = np.zeros_like(self.x)

        omega = self.speed  # Angular frequency
        
        # Calculate wave in medium 1 with angular propagation
        mask1 = self.x <= self.boundary1
        x1 = self.x[mask1]
        phase1 = k1 * x1 - self.speed * self.time
        wave[mask1] = self.amplitude * self.visualization_scale * np.sin(phase1)
        
        # Calculate wave in medium 2 with angular propagation
        mask2 = (self.x > self.boundary1) & (self.x <= self.boundary2)
        x2 = self.x[mask2] - self.boundary1
        
        # Handle differently based on whether we have total internal reflection
        if 'total_reflection' in locals() and total_reflection and sin_refraction1 > 1:
            # For total internal reflection at first boundary, create reflected wave
            phase2 = k1 * x2 + k1 * self.boundary1 - self.speed * self.time
        else:
            # Normal refraction
            phase2 = k2 * x2 + k1 * self.boundary1 - self.speed * self.time
                  
        wave[mask2] = self.amplitude * self.visualization_scale * np.sin(phase2)
        
        # Calculate wave in medium 3 with angular propagation
        mask3 = self.x > self.boundary2
        x3 = self.x[mask3] - self.boundary2
        
        
        # Handle differently based on whether we have total internal reflection
        if 'total_reflection' in locals() and total_reflection:
            if sin_refraction1 > 1:
                # Total reflection at first boundary means no wave in medium 3
                wave[mask3] = 0
            elif sin_refraction2 > 1:
                # Total reflection at second boundary
                phase3 = k2 * x3 + k2 * (self.boundary2 - self.boundary1) + k1 * self.boundary1 - self.speed * self.time
                wave[mask3] = self.amplitude * self.visualization_scale * np.sin(phase3)
        else:
            # Normal refraction to medium 3
            phase3 = k3 * x3 + k2 * (self.boundary2 - self.boundary1) + k1 * self.boundary1 - self.speed * self.time
            wave[mask3] = self.amplitude * self.visualization_scale * np.sin(phase3)
        
        if not hasattr(self, '_wave_cache'):
            self._wave_cache = {}
        self._wave_cache['key'] = cache_key
        self._wave_cache['wave'] = wave

        return wave
        

    def update_animation(self):
        """Update the animation for each timer tick"""
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
                'frequency': self.frequency
            }
            needs_update = True
        
        # Check if any relevant state has changed
        if (abs(self._prev_state['time'] - self.time) > 0.009 or
            self._prev_state['white_light'] != self.white_light or
            self._prev_state['show_interference'] != self.show_interference or
            self._prev_state['frequency'] != self.frequency):
            needs_update = True
        
        if needs_update:
            # Update the wave data
            if self.white_light:
                # Update multiple waves with different wavelengths
                for curve, freq in self.wave_curves:
                    wave = self.calculate_wave(freq)
                    curve.setData(self.x, wave)
            else:
                # Update single wave
                wave = self.calculate_wave(self.frequency)
                self.wave_curve.setData(self.x, wave)
            
            # Update interference waves if enabled
            if self.show_interference:
                wave1, wave2 = self.calculate_interference_waves(self.wavelength)
                self.interference_wave1.setData(self.x, wave1)
                self.interference_wave2.setData(self.x, wave2)
                
                # Ensure they're visible
                if not self.interference_wave1.isVisible():
                    self.interference_wave1.setVisible(True)
                if not self.interference_wave2.isVisible():
                    self.interference_wave2.setVisible(True)
            else:
                # Ensure they're hidden when not enabled
                if self.interference_wave1.isVisible():
                    self.interference_wave1.setVisible(False)
                if self.interference_wave2.isVisible():
                    self.interference_wave2.setVisible(False)
            
            # Update previous state
            self._prev_state = {
                'time': self.time,
                'white_light': self.white_light,
                'show_interference': self.show_interference,
                'frequency': self.frequency
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

        self.interference_wave1.setVisible(enabled)
        self.interference_wave2.setVisible(enabled)
        
        if enabled:
            # Force a significant time update to ensure waves are animated
            old_time = self.time
            self.time += 2.5  # Add a significant time offset
            
            # Force recalculation with the new time
            wave1, wave2 = self.calculate_interference_waves(self.wavelength)
            
            # Update the interference waves
            self.interference_wave1.setData(self.x, wave1)
            self.interference_wave2.setData(self.x, wave2)
            
            
            # Reset time to original plus a small increment to keep animation flowing
            self.time = old_time + 0.1
        else:
            # Hide interference waves when disabled
            self.interference_wave1.setVisible(False)
            self.interference_wave2.setVisible(False)

    def update_plot(self):
        """Update the plot with current parameters"""
        try:
            # Store current y range to restore after updates
            y_range = self.getViewBox().viewRange()[1]
            
            # Clear existing medium rectangles and labels
            for rect in self.medium_rects:
                self.removeItem(rect)
            self.medium_rects = []
            
            for label_pair in self.medium_labels:
                for label in label_pair[0:2]:
                    self.removeItem(label)
            self.medium_labels = []
            
            # Create medium rectangles with proper fill
            # Medium 1 (Air)
        # Force color to be a QColor object
            medium1_color = QColor(self.medium1_color)
            medium1_rect = pg.FillBetweenItem(
                pg.PlotCurveItem([0, self.boundary1], [2, 2]),
                pg.PlotCurveItem([0, self.boundary1], [-2, -2]),
                brush=pg.mkBrush(medium1_color))
            self.addItem(medium1_rect)
            self.medium_rects.append(medium1_rect)
            
            # Medium 2 (Water)
 
            # Force color to be a QColor object
            medium2_color = QColor(self.medium2_color)
            medium2_rect = pg.FillBetweenItem(
                pg.PlotCurveItem([self.boundary1, self.boundary2], [2, 2]),
                pg.PlotCurveItem([self.boundary1, self.boundary2], [-2, -2]),
                brush=pg.mkBrush(medium2_color))
            self.addItem(medium2_rect)
            self.medium_rects.append(medium2_rect)
            
            # Medium 3 (Glass)
            # Force color to be a QColor object
            medium3_color = QColor(self.medium3_color)
            medium3_rect = pg.FillBetweenItem(
                pg.PlotCurveItem([self.boundary2, 3000], [2, 2]),
                pg.PlotCurveItem([self.boundary2, 3000], [-2, -2]),
                brush=pg.mkBrush(medium3_color))
            self.addItem(medium3_rect)
            self.medium_rects.append(medium3_rect)
        

            # Create wavelength scale bar and labels
            self.scale_bar = None
            self.scale_label = None
            self.scale_title = None
            
            # Create initial scale bar (will be updated with wavelength changes)
            self.update_wavelength_scale()
            
            self.medium_rects.append(self.boundary1_line)
            self.medium_rects.append(self.boundary2_line)
            
            # Add boundary lines
            boundary1_line = pg.InfiniteLine(pos=self.boundary1, angle=90, pen=pg.mkPen('w', width=2, style=Qt.DashLine))
            boundary2_line = pg.InfiniteLine(pos=self.boundary2, angle=90, pen=pg.mkPen('w', width=2, style=Qt.DashLine))
            self.addItem(boundary1_line)
            self.addItem(boundary2_line)
            self.medium_rects.append(boundary1_line)
            self.medium_rects.append(boundary2_line)
            
             # Bring grid lines to front by removing and re-adding them
            for line in self.grid_lines_x + self.grid_lines_y:
                self.removeItem(line)
                self.addItem(line)

            grid_pen = pg.mkPen(color=(255, 255, 255, 80), width=0.5, style=Qt.DotLine)
            self.showGrid(x=True, y=True, alpha=0.5)
            # Update medium labels
            # Medium 1
            medium1_color_label = pg.TextItem(f"n₁: {self.n1:.4f}", anchor=(0.5, 0), color='white')
            medium1_color_label.setPos(self.boundary1 / 2, -1.2)
            self.addItem(medium1_color_label)
            
            medium1_name_label = pg.TextItem(f"{self.medium1_name}", anchor=(0.5, 1), color='white')
            medium1_name_label.setPos(self.boundary1 / 2, 1.8)
            self.addItem(medium1_name_label)
            
            self.medium_labels.append((medium1_color_label, medium1_name_label, 1))
            
            # Medium 2
            medium2_color_label = pg.TextItem(f"n₂: {self.n2:.4f}", anchor=(0.5, 0), color='white')
            medium2_color_label.setPos(self.boundary1 + (self.boundary2 - self.boundary1) / 2, -1.2)
            self.addItem(medium2_color_label)
            
            medium2_name_label = pg.TextItem(f"{self.medium2_name}", anchor=(0.5, 1), color='white')
            medium2_name_label.setPos(self.boundary1 + (self.boundary2 - self.boundary1) / 2, 1.8)
            self.addItem(medium2_name_label)
            
            self.medium_labels.append((medium2_color_label, medium2_name_label, 2))
            
            
            # Medium 3
            medium3_color_label = pg.TextItem(f"n₃: {self.n3:.4f}", anchor=(0.5, 0), color='white')
            medium3_color_label.setPos(self.boundary2 + (3000 - self.boundary2) / 2, -1.2)
            self.addItem(medium3_color_label)
            
            medium3_name_label = pg.TextItem(f"{self.medium3_name}", anchor=(0.5, 1), color='white')
            medium3_name_label.setPos(self.boundary2 + (3000 - self.boundary2) / 2, 1.8)
            self.addItem(medium3_name_label)
            
            self.medium_labels.append((medium3_color_label, medium3_name_label, 3))
            
            # Update the wave curve with current parameters
            if not self.white_light:
                # Get color from wavelength
                r, g, b = wavelength_to_rgb(299792.458 / self.frequency)
                wave_color = QColor(r, g, b)
                self.wave_curve.setPen(pg.mkPen(wave_color, width=4))
                
                # Update wave data
                if hasattr(self, 'wave_curve'):
                    self.wave_curve.setPen(pg.mkPen(wave_color, width=4))
                    wave = self.calculate_wave(self.frequency)
                    self.wave_curve.setData(self.x, wave)
            else:
                # Update white light curves
                if hasattr(self, 'wave_curves'):
                    for curve, wl in self.wave_curves:
                        r, g, b = wavelength_to_rgb(wl)
                        curve.setPen(pg.mkPen(QColor(r, g, b), width=2))
                        wave = self.calculate_wave(wl)
                        curve.setData(self.x, wave)
            
            # Update ray lines if ray mode is enabled
            if self.show_ray_mode:
                self.update_ray_lines()
                
            # Restore y range to prevent auto-scaling
            self.setYRange(y_range[0], y_range[1], padding=0)
            
        except Exception as e:
            print(f"Error in update_plot: {str(e)}")
            import traceback
            traceback.print_exc()
            return

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
        # Remove existing scale elements if they exist

        scale_key = (self.frequency, self.n1, self.n2, self.n3)
    
        # Skip update if nothing has changed
        if hasattr(self, '_scale_key') and self._scale_key == scale_key:
            return
        
        # Store new key
        self._scale_key = scale_key

        if hasattr(self, 'scale_bar') and self.scale_bar is not None:
            self.removeItem(self.scale_bar)
        if hasattr(self, 'scale_label') and self.scale_label is not None:
            self.removeItem(self.scale_label)
        if hasattr(self, 'scale_title') and self.scale_title is not None:
            self.removeItem(self.scale_title)
            
        # Remove existing grid value labels if they exist
        if hasattr(self, 'grid_value_labels'):
            for label in self.grid_value_labels:
                self.removeItem(label)
        
        # Initialize grid value labels list
        self.grid_value_labels = []

        vacuum_wavelength_nm = 299792.458 / self.frequency
        wavelength_m1 = vacuum_wavelength_nm / self.n1
        wavelength_m2 = vacuum_wavelength_nm / self.n2
        wavelength_m3 = vacuum_wavelength_nm / self.n3

        c = 299792458  
        speed_m1 = c / self.n1  
        speed_m2 = c / self.n2  
        speed_m3 = c / self.n3

        speed_m1_formatted = f"{1/self.n1:.2f}c"
        speed_m2_formatted = f"{1/self.n2:.2f}c"
        speed_m3_formatted = f"{1/self.n3:.2f}c"

        scale_factor = 500 / 500
        # Add wavelength values to vertical grid lines
        for x in np.arange(0, 3001, 500):
            # Calculate corresponding wavelength value
            physical_nm = x / scale_factor  # Scale to make one grid unit = wavelength/500
            
            # Create label with wavelength value
            label = pg.TextItem(f"{int(physical_nm)} nm", anchor=(0, 0.5), color='white')
            label.setPos(x+10, -1.9)  # Position below the plot
            self.addItem(label)
            self.grid_value_labels.append(label)

        # Add vacuum wavelength reference
        wavelength_info = pg.TextItem(f"λ₀ = {int(vacuum_wavelength_nm)} nm", 
                                     anchor=(0, 0), color='yellow')
        wavelength_info.setPos(50, 1.5)  # Position at top-left
        self.addItem(wavelength_info)
        self.grid_value_labels.append(wavelength_info)
        
         # Add medium-specific wavelength and speed labels
        m1_label = pg.TextItem(f"λ₁ = {int(wavelength_m1)} nm | v₁ = {speed_m1_formatted}", 
                              anchor=(0.5, 0), color='yellow')
        m1_label.setPos(self.boundary1/2, 1.5)
        self.addItem(m1_label)
        self.grid_value_labels.append(m1_label)
        
        m2_label = pg.TextItem(f"λ₂ = {int(wavelength_m2)} nm | v₂ = {speed_m2_formatted}", 
                              anchor=(0.5, 0), color='yellow')
        m2_label.setPos(self.boundary1 + (self.boundary2-self.boundary1)/2, 1.5)
        self.addItem(m2_label)
        self.grid_value_labels.append(m2_label)
        
        m3_label = pg.TextItem(f"λ₃ = {int(wavelength_m3)} nm | v₃ = {speed_m3_formatted}", 
                              anchor=(0.5, 0), color='yellow')
        m3_label.setPos(self.boundary2 + (3000-self.boundary2)/2, 1.5)
        self.addItem(m3_label)
        self.grid_value_labels.append(m3_label)

    def update_wavelength(self, value):
        """Update the wavelength and colors"""
        self.frequency = value
        self.wavelength = 299792.458 / value

        if hasattr(self, '_wave_cache'):
            del self._wave_cache
        if hasattr(self, '_interference_cache'):
            del self._interference_cache
        if hasattr(self, '_scale_key'):
            del self._scale_key

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
        
        if hasattr(self, '_wave_cache'):
            del self._wave_cache
        if hasattr(self, '_interference_cache'):
            del self._interference_cache

        # Store current y-range before updating
        y_range = self.getPlotItem().getViewBox().viewRange()[1]
        
        # Update the wave with new amplitude
        wave = self.calculate_wave(self.wavelength)
        self.wave_curve.setData(self.x, wave)
        
        # Restore the previous y-range to prevent auto-zooming
        self.getPlotItem().getViewBox().setYRange(y_range[0], y_range[1], padding=0)

    def update_speed(self, value):
        self.speed = value
        
        if hasattr(self, '_wave_cache'):
            del self._wave_cache
        if hasattr(self, '_interference_cache'):
            del self._interference_cache

    def update_n1(self, value):
        """Update refractive index of medium 1"""
        self.n1 = value
        # Invalidate caches when parameters change
        if hasattr(self, '_wave_cache'):
            del self._wave_cache
        if hasattr(self, '_interference_cache'):
            del self._interference_cache
        if hasattr(self, '_scale_key'):
            del self._scale_key
        self.update_plot()
        
    def update_n2(self, value):
        """Update refractive index of medium 2"""
        self.n2 = value
        # Invalidate caches when parameters change
        if hasattr(self, '_wave_cache'):
            del self._wave_cache
        if hasattr(self, '_interference_cache'):
            del self._interference_cache
        if hasattr(self, '_scale_key'):
            del self._scale_key
        self.update_plot()
        
    def update_n3(self, value):
        """Update refractive index of medium 3"""
        self.n3 = value
        # Invalidate caches when parameters change
        if hasattr(self, '_wave_cache'):
            del self._wave_cache
        if hasattr(self, '_interference_cache'):
            del self._interference_cache
        if hasattr(self, '_scale_key'):
            del self._scale_key
        self.update_plot()
        
    def update_medium1(self, medium_name):
        if medium_name in self.medium_presets:
            # Update the refractive index
            n1 = self.medium_presets[medium_name]['n']
            self.n1_slider.setValue(int(n1 * 100))
            
            # Update the wave widget's n1 value and name
            self.wave_widget.update_n1(n1)
            self.wave_widget.medium1_name = medium_name
            
            # Update the medium color
            self.wave_widget.medium1_color = self.medium_presets[medium_name]['color']
        
            self.wave_widget.update_plot()
            
    def update_medium2(self, medium_name):

        if medium_name in self.medium_presets:
            # Update the refractive index
            n2 = self.medium_presets[medium_name]['n']
            self.n2_slider.setValue(int(n2 * 100))
            
            # Update the wave widget's n2 value and name
            self.wave_widget.update_n2(n2)
            self.wave_widget.medium2_name = medium_name
            
            # Update the medium color
            self.wave_widget.medium2_color = self.medium_presets[medium_name]['color']

            self.wave_widget.update_plot()
            
    def update_medium3(self, medium_name):
   
        if medium_name in self.medium_presets:
            # Update the refractive index
            n3 = self.medium_presets[medium_name]['n']
            self.n3_slider.setValue(int(n3 * 100))
            
            # Update the wave widget's n3 value and name
            self.wave_widget.update_n3(n3)
            self.wave_widget.medium3_name = medium_name
            
            # Update the medium color
            self.wave_widget.medium3_color = self.medium_presets[medium_name]['color']
            
            self.wave_widget.update_plot()
        
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

            if hasattr(self, 'wave_curve'):
                self.wave_curve.setVisible(False)
            # Hide the original single wavelength curve when in white light mode
            for curve_item in self.wave_curves:
                if isinstance(curve_item, tuple) and len(curve_item) == 2:
                    curve, _ = curve_item
                    if hasattr(curve, 'setVisible'):
                        curve.setVisible(True)
        else:
            # Show the original single wavelength curve when not in white light mode
            if hasattr(self, 'wave_curve'):
                self.wave_curve.setVisible(True)
        
             # Hide all white light curves
            for curve_item in self.wave_curves:
                if isinstance(curve_item, tuple) and len(curve_item) == 2:
                    curve, _ = curve_item
                    if hasattr(curve, 'setVisible'):
                        curve.setVisible(False)

    def create_white_light_curves(self):
        """Create curves for multiple wavelengths to simulate white light"""
        for curve_item in self.wave_curves:
            if isinstance(curve_item, tuple) and len(curve_item) == 2:
                curve, _ = curve_item
                self.removeItem(curve)
      
        self.wave_curves = []
        # Clear existing curves if any
        for freq in self.prism_frequencies:
            # Calculate wavelength from frequency
            wl = 299792.458 / freq
                
        # Calculate wave for this frequency
            wave = self.calculate_wave(freq)
            
            # Get color for this wavelength
            r, g, b = wavelength_to_rgb(wl)
            wave_color = QColor(r, g, b)
            
            # Create curve with this color
            curve = self.plot(self.x, wave, pen=pg.mkPen(wave_color, width=4))
            
             # Store curve and frequency
            self.wave_curves.append((curve, freq))

    def update_ray_lines(self):
        """Update the ray lines to visualize refraction angles"""
        try:
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
        if self.show_ray_mode:
            self.update_ray_lines()

    def set_ray_target_y(self, y_value):
        """Set the Y-coordinate where the ray hits the first boundary"""
        self.ray_target_y = y_value
        if self.show_ray_mode:
            self.update_ray_lines()

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
        self.setGeometry(100, 100, 1200, 800)
        

        # Medium presets (refractive indices at ~550nm wavelength)
        self.medium_presets = {
            'Air': {'n': 1.0003, 'color': '#3A3A3A80'},  # Very dark gray with transparency
            'Water': {'n': 1.33, 'color': '#0D47A1C0'},  # Rich blue with higher opacity
            'Glass (Crown)': {'n': 1.52, 'color': '#4488AACC'},  # Deep cyan with higher opacity
            'Glass (Flint)': {'n': 1.62, 'color': '#55AABBCC'},  # Dark cyan with higher opacity
            'Diamond': {'n': 2.42, 'color': '#331a37C0'},  # Dark teal with higher opacity
            'Acrylic': {'n': 1.49, 'color': '#1565C0C0'},  # Bright blue with higher opacity
            'Glycerine': {'n': 1.47, 'color': '#2f2016C0'},  # Light blue with higher opacity
            'Ethanol': {'n': 1.36, 'color': '#66CCFFCC'},  # Lighter blue with higher opacity
            'Quartz': {'n': 1.54, 'color': '#200e0eC0'},  # Medium cyan with higher opacity
            'Sapphire': {'n': 1.77, 'color': '#311B92C0'}  # Deep purple with higher opacity
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
        
        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        self.main_layout = QVBoxLayout(central_widget)
        self.main_layout.setContentsMargins(0, 0, 0, 0)  # Remove margins
        self.main_layout.setSpacing(0)  # Remove spacing
        
        # Create wave widget directly (no tabs)
        self.wave_widget = WaveSimulationWidget()
        self.main_layout.addWidget(self.wave_widget)
        
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

        # Add interference mode checkbox next to other mode controls
        mode_layout = QHBoxLayout()
        self.white_light_check = QCheckBox("White Light")
        self.interference_check = QCheckBox("Show Interference")
        self.interference_check.setChecked(False)
        self.ray_mode_check = QCheckBox("Show Ray Path")
        
        mode_layout.addWidget(self.white_light_check)
        mode_layout.addWidget(self.interference_check)
        mode_layout.addWidget(self.ray_mode_check)
        wave_layout.addLayout(mode_layout)

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
        self.medium3_combo.setCurrentText('Glass (Crown)')
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
        """Update the angle of incidence"""
        # Update label
        self.angle_value.setText(f"{value}°")
        # Update angle in wave widget
        self.wave_widget.angle_of_incidence = value
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



# Run the application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    mainWindow = LightSimulationApp()
    mainWindow.show()
    sys.exit(app.exec_())