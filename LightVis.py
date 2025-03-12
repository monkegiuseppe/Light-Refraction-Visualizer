import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                            QLabel, QSlider, QPushButton, QComboBox, QGroupBox, QFrame,
                            QCheckBox, QTabWidget, QRadioButton, QButtonGroup)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QPainter,QFont, QColor, QPen, QBrush
import pyqtgraph as pg

# Define wavelength to RGB color mapping
def wavelength_to_rgb(wavelength):
    """Convert wavelength in nm to RGB color"""
    # Map simulation wavelength to approximate visible spectrum (380-750nm)
    # Assuming wavelength 10-100 in simulation maps to 380-750nm
    nm = wavelength
    
    # Based on algorithm by Dan Bruton (www.physics.sfasu.edu/astro/color/spectra.html)
    if 380 <= nm < 440:
        r = -(nm - 440) / (440 - 380)
        g = 0.0
        b = 1.0
    elif 440 <= nm < 490:
        r = 0.0
        g = (nm - 440) / (490 - 440)
        b = 1.0
    elif 490 <= nm < 510:
        r = 0.0
        g = 1.0
        b = -(nm - 510) / (510 - 490)
    elif 510 <= nm < 580:
        r = (nm - 510) / (580 - 510)
        g = 1.0
        b = 0.0
    elif 580 <= nm < 645:
        r = 1.0
        g = -(nm - 645) / (645 - 580)
        b = 0.0
    elif 645 <= nm <= 750:
        r = 1.0
        g = 0.0
        b = 0.0
    else:
        r, g, b = 0.5, 0.5, 0.5  # Outside visible spectrum
    
    # Attenuate brightness at the edges of the visible spectrum
    if 380 <= nm < 420:
        factor = 0.3 + 0.7 * (nm - 380) / (420 - 380)
    elif 420 <= nm < 700:
        factor = 1.0
    elif 700 <= nm <= 750:
        factor = 0.3 + 0.7 * (750 - nm) / (750 - 700)
    else:
        factor = 0.0
        
    r = round(255 * r * factor)
    g = round(255 * g * factor)
    b = round(255 * b * factor)
    
    return (r, g, b)

class WaveSimulationWidget(pg.PlotWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent

        # Create font
        plot_font = QFont("Foto")

        # Set up the plot
        self.setAntialiasing(True)
        self.setBackground('k')
        self.setLabel('left', 'Amplitude', color='w')
        self.setLabel('bottom', 'Position', color='w')
        
        # Ensure background stays black regardless of parent styling
        self.setStyleSheet("background-color: black; border: none;")

        # Add grid
        self.showGrid(x=True, y=True, alpha=0.3)
        self.getPlotItem().getAxis('left').setPen('w')
        self.getPlotItem().getAxis('bottom').setPen('w')

        # Hide the auto-range button
        self.getPlotItem().hideButtons()
        
        # Remove the border and axis rectangle
        self.getPlotItem().getViewBox().setBackgroundColor(None)
        self.getPlotItem().getViewBox().setBorder(None)
        
        # Hide the axis box and remove borders
        for axis in ['left', 'bottom', 'top', 'right']:
            if self.getPlotItem().getAxis(axis):
                self.getPlotItem().getAxis(axis).setStyle(showValues=True, tickLength=5)
                if axis in ['top', 'right']:
                    self.getPlotItem().showAxis(axis, False)
        
        # Remove the plot border
        self.getPlotItem().setContentsMargins(0, 0, 0, 0)
        self.getPlotItem().layout.setContentsMargins(0, 0, 0, 0)
        self.getPlotItem().vb.setContentsMargins(0, 0, 0, 0)

        # Set grid line spacing
        self.getPlotItem().getAxis('left').setTicks([[(i, str(i)) for i in range(-2, 3, 1)]])
        self.getPlotItem().getAxis('bottom').setTicks([[(i, str(i)) for i in range(0, 3001, 500)]])

        # Set up the x-axis
        self.x = np.linspace(0, 3000, 3000)

        # Define boundaries first
        self.boundary1 = 1000
        self.boundary2 = 2000

        # Initial parameters
        self.wavelength = 550
        self.amplitude = 5
        self.speed = 2
        self.visualization_scale = 0.1  # Add visualization scaling factor
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
        
        # Create angle labels (initialize them here to avoid errors)
        self.angle1_label = pg.TextItem("θ₁: 0°", anchor=(0, 0), color='white')
        self.angle2_label = pg.TextItem("θ₂: 0°", anchor=(0, 0), color='white')
        self.angle3_label = pg.TextItem("θ₃: 0°", anchor=(0, 0), color='white')
        self.addItem(self.angle1_label)
        self.addItem(self.angle2_label)
        self.addItem(self.angle3_label)
        self.angle1_label.setPos(500, -5)
        self.angle2_label.setPos(1500, -5)
        self.angle3_label.setPos(2500, -5)

        # Add interference wave curves with specific colors
        self.interference_wave1 = self.plot(self.x, np.zeros_like(self.x), 
                                          pen=pg.mkPen('r', width=2, style=Qt.DashLine))  # Red for original wave (dashed)
        self.interference_wave2 = self.plot(self.x, np.zeros_like(self.x), 
                                          pen=pg.mkPen('g', width=2, style=Qt.DashLine))  # Green for difference wave (dashed)
        
        # Initialize interference waves with proper data
        # Set a non-zero time to ensure waves are not straight lines
        self.time = 2.5  # Use a larger time value for initialization
        wave1, wave2 = self.calculate_interference_waves(self.wavelength)
        self.interference_wave1.setData(self.x, wave1)
        self.interference_wave2.setData(self.x, wave2)
        self.time = 0  # Reset time after initialization
        
        # Hide interference waves initially
        self.interference_wave1.setVisible(False)
        self.interference_wave2.setVisible(False)
        self.show_interference = False
        
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

        # Define boundaries
        self.boundary1 = 1000
        self.boundary2 = 2000
        
        # Add angle labels at boundaries 
        self.angle1_label = pg.TextItem("θ₁: 0°", anchor=(0.5, 0), color='#ffffff')
        self.angle1_label.setPos(self.boundary1 - 200, -1.5)  # Moved further left
        # Set a larger, bolder font for better visibility
        font = QFont()
        font.setPointSize(12)
        font.setBold(True)
        self.angle1_label.setFont(font)
        # Add a semi-transparent background to the text for better visibility
        self.angle1_label.fill = pg.mkBrush(0, 0, 0, 150)  # Semi-transparent black

        self.angle2_label = pg.TextItem("θ₂: 0°", anchor=(0.5, 0), color='#ffffff')
        self.angle2_label.setPos(self.boundary1 + 200, -1.5)  # Moved further right
        self.angle2_label.setFont(font)
        self.angle2_label.fill = pg.mkBrush(0, 0, 0, 150)  # Semi-transparent black

        self.angle3_label = pg.TextItem("θ₃: 0°", anchor=(0.5, 0), color='#ffffff')
        self.angle3_label.setPos(self.boundary2 + 200, -1.5)  # Moved further right
        self.angle3_label.setFont(font)
        self.angle3_label.fill = pg.mkBrush(0, 0, 0, 150)  # Semi-transparent black

        self.addItem(self.angle1_label)
        self.addItem(self.angle2_label)
        self.addItem(self.angle3_label)

        # Simulation mode
        self.prism_mode = False
        self.white_light = False
        
    
        
        
        # Set fixed Y-axis range to show waves better
        self.setXRange(0, 3000, padding=0)
        self.setYRange(-2, 2, padding=0)

        # Disable mouse interaction
        self.setMouseEnabled(x=False, y=False)
        # Disable auto-range
        self.enableAutoRange(axis='x', enable=False)
        self.enableAutoRange(axis='y', enable=False)

        # Medium colors
        self.medium_colors = {
            'Air': (230, 230, 255, 50),        # Very light blue
            'Water': (153, 204, 255, 80),      # Light blue
            'Glass (Crown)': (204, 230, 230, 100),  # Light cyan
            'Glass (Flint)': (179, 204, 204, 100),  # Darker cyan
            'Diamond': (242, 242, 255, 130),   # White/blue tint
            'Acrylic': (230, 230, 179, 80),    # Light yellow
            'Glycerine': (230, 204, 230, 80),  # Light purple
            'Ethanol': (204, 204, 230, 80),    # Light purple-blue
            'Quartz': (255, 255, 230, 80),     # Very light yellow
            'Sapphire': (179, 179, 230, 100)   # Medium blue
        }
        
        # Current medium selections
        self.current_medium1 = 'Air'
        self.current_medium2 = 'Water'
        self.current_medium3 = 'Glass (Crown)'
        
        # Create the medium regions (transparent colored rectangles)
        self.medium1_region = pg.LinearRegionItem([0, self.boundary1], movable=False, 
                                                brush=QBrush(QColor(*self.medium_colors['Air'])))
        self.medium2_region = pg.LinearRegionItem([self.boundary1, self.boundary2], movable=False, 
                                                brush=QBrush(QColor(*self.medium_colors['Water'])))
        self.medium3_region = pg.LinearRegionItem([self.boundary2, 3000], movable=False, 
                                                brush=QBrush(QColor(*self.medium_colors['Glass (Crown)'])))
        
        # Add regions to plot
        self.addItem(self.medium1_region)
        self.addItem(self.medium2_region)
        self.addItem(self.medium3_region)
        
        # Add boundary lines
        self.boundary1_line = pg.InfiniteLine(pos=self.boundary1, angle=90, pen=pg.mkPen('w', width=2, style=Qt.DashLine))
        self.boundary2_line = pg.InfiniteLine(pos=self.boundary2, angle=90, pen=pg.mkPen('w', width=2, style=Qt.DashLine))
        self.addItem(self.boundary1_line)
        self.addItem(self.boundary2_line)
        
        # Add medium labels
        self.medium1_label = pg.TextItem(f'{self.current_medium1} (n₁ = {self.n1:.4f})', anchor=(0.5, 0), color='w')
        self.medium1_label.setPos(self.boundary1/2, 40)
        self.medium1_label.setFont(QFont("Foto", 12, QFont.Bold))
        
        self.medium2_label = pg.TextItem(f'{self.current_medium2} (n₂ = {self.n2:.4f})', anchor=(0.5, 0), color='w')
        self.medium2_label.setPos(self.boundary1 + (self.boundary2-self.boundary1)/2, 40)
        self.medium2_label.setFont(QFont("Foto", 12, QFont.Bold))
        
        self.medium3_label = pg.TextItem(f'{self.current_medium3} (n₃ = {self.n3:.4f})', anchor=(0.5, 0), color='w')
        self.medium3_label.setPos(self.boundary2 + (3000-self.boundary2)/2, 40)
        self.medium3_label.setFont(QFont("Foto", 12, QFont.Bold))
        
        self.addItem(self.medium1_label)
        self.addItem(self.medium2_label)
        self.addItem(self.medium3_label)

        
       
        
        # For white light/prism mode - create multiple curves for different wavelengths
        self.wave_curves = []
        self.prism_wavelengths = [400, 450, 500, 550, 600, 650, 700]  # Different wavelengths
        
        for wl in self.prism_wavelengths:
            # Calculate wave for this wavelength
            wave = self.calculate_wave(wl)
            
            # Get color for this wavelength
            color = wavelength_to_rgb(wl)
            
            # Create curve with this color
            curve = self.plot(self.x, wave, pen=pg.mkPen(color, width=4))
            
            # Hide initially
            curve.setVisible(False)
            
            # Add to list
            self.wave_curves.append((wl, curve))
        
        # Set up the animation timer
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_animation)
        self.timer.start(10)  # 50ms interval (20 fps)

         # Create the wave curve for single wavelength
        wave = self.calculate_wave(self.wavelength)
        # Use wavelength_to_rgb to set the initial color based on self.wavelength
        r, g, b = wavelength_to_rgb(self.wavelength)
        wave_color = QColor(r, g, b)
        self.wave_curve = self.plot(self.x, wave, pen=pg.mkPen(wave_color, width=4))

        # Medium setup
        self.boundary1 = 1000
        self.boundary2 = 2000
        
        # Initialize medium properties
        self.n1 = 1.0
        self.n2 = 1.0
        self.n3 = 1.0
        
        # Initialize medium colors with default values
        self.medium1_color = '#222233'  # Dark blue-gray
        self.medium2_color = '#223344'  # Medium blue-gray
        self.medium3_color = '#334455'  # Light blue-gray
        
        # Initialize containers for medium visualization
        self.medium_rects = []  # Store rectangle items for medium backgrounds
        self.medium_labels = []  # Store labels for media
        
        # Initialize medium presets
        self.medium_presets = {
            'Air': {'n': 1.0003, 'color': '#222233'},
            'Water': {'n': 1.33, 'color': '#2244AA'},
            'Glass (Crown)': {'n': 1.52, 'color': '#445566'},
            'Glass (Flint)': {'n': 1.62, 'color': '#556677'},
            'Diamond': {'n': 2.42, 'color': '#778899'},
            'Acrylic': {'n': 1.49, 'color': '#223355'},
            'Glycerine': {'n': 1.47, 'color': '#336699'},
            'Ethanol': {'n': 1.36, 'color': '#3377AA'},
            'Quartz': {'n': 1.54, 'color': '#667788'},
            'Sapphire': {'n': 1.77, 'color': '#5566AA'}
        }
    
    def calculate_interference_waves(self, wavelength):
        """Calculate the component waves that create interference"""
        # Get the actual refracted wave
        actual_wave = self.calculate_wave(wavelength)
        
        # Wave 1 (red) - original wave continuing straight without refraction
        wave1 = np.zeros_like(self.x)
        k1 = (2 * np.pi * self.n1) / wavelength
        
        # Calculate phase with time component to ensure animation
        phase1 = k1 * self.x - self.speed * self.time
        
        # In medium 1, the actual wave and the straight wave are identical
        # (when angle of incidence is 0)
        mask_medium1 = self.x <= self.boundary1
        
        # For medium 1, if angle is 0, there's no difference between actual and straight
        if abs(self.angle_of_incidence) < 0.001:  # Practically zero
            wave1[mask_medium1] = actual_wave[mask_medium1]
        else:
            # If angle is not 0, calculate the straight wave in medium 1
            wave1[mask_medium1] = self.amplitude * self.visualization_scale * np.sin(phase1[mask_medium1])
        
        # For medium 2 and 3, calculate the straight wave as if it continued from medium 1
        mask_medium23 = self.x > self.boundary1
        wave1[mask_medium23] = self.amplitude * self.visualization_scale * np.sin(phase1[mask_medium23])
        
        # Wave 2 (green) - the difference between actual wave and original wave
        wave2 = actual_wave - wave1
        
        return wave1, wave2
        
    def calculate_wave(self, wavelength):
        """Calculate the wave based on current parameters for a specific wavelength"""
        wave = np.zeros_like(self.x)
        
        # We need to adjust refractive indices for wavelength in prism mode
        # This simulates dispersion - different wavelengths refract differently
        if self.prism_mode:
            # Calculate approximate wavelength in nm (for dispersion calculation)
            # Map wavelength 10-100 to 380-750nm (visible spectrum)
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
            self.angle1_label.setText(f"θ₁: {self.angle_of_incidence}°")
            self.angle2_label.setText(f"θ₂: {np.degrees(angle_refraction1):.1f}°")
            self.angle3_label.setText(f"θ₃: {np.degrees(angle_refraction2):.1f}°")

            
        # Wave parameters for each medium
        k1 = (2 * np.pi * self.n1) / wavelength  # Wave number in medium 1
        k2 = (2 * np.pi * n2_wl) / wavelength    # Wave number in medium 2
        k3 = (2 * np.pi * n3_wl) / wavelength    # Wave number in medium 3
        
        # Create y-coordinates for vertical displacement
        y = np.linspace(-1, 1, len(self.x))

        omega = self.speed  # Angular frequency
        
        # Calculate wave in medium 1 with angular propagation
        mask1 = self.x <= self.boundary1
        x1 = self.x[mask1]
        phase1 = k1 * (x1 * np.cos(angle_incidence) + y[mask1] * np.sin(angle_incidence)) - self.speed * self.time
        wave[mask1] = self.amplitude * self.visualization_scale * np.sin(phase1)
        
        # Calculate wave in medium 2 with angular propagation
        mask2 = (self.x > self.boundary1) & (self.x <= self.boundary2)
        x2 = self.x[mask2] - self.boundary1
        
        # Handle differently based on whether we have total internal reflection
        if 'total_reflection' in locals() and total_reflection and sin_refraction1 > 1:
            # For total internal reflection at first boundary, create reflected wave
            reflected_angle = np.pi - angle_incidence
            phase2 = (k1 * (x2 * np.cos(reflected_angle) + y[mask2] * np.sin(reflected_angle)) + 
                  k1 * self.boundary1 * np.cos(angle_incidence)) - self.speed * self.time
        else:
            # Normal refraction
            phase2 = (k2 * (x2 * np.cos(angle_refraction1) + y[mask2] * np.sin(angle_refraction1)) + 
                  k1 * self.boundary1 * np.cos(angle_incidence)) - self.speed * self.time
                  
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
                reflected_angle = np.pi - angle_refraction1
                phase3 = (k2 * (x3 * np.cos(reflected_angle) + y[mask3] * np.sin(reflected_angle)) + 
                      k2 * (self.boundary2 - self.boundary1) * np.cos(angle_refraction1) +
                      k1 * self.boundary1 * np.cos(angle_incidence)) - self.speed * self.time
                wave[mask3] = self.amplitude * self.visualization_scale * np.sin(phase3)
        else:
            # Normal refraction to medium 3
            phase3 = (k3 * (x3 * np.cos(angle_refraction2) + y[mask3] * np.sin(angle_refraction2)) + 
                  k2 * (self.boundary2 - self.boundary1) * np.cos(angle_refraction1) +
                  k1 * self.boundary1 * np.cos(angle_incidence)) - self.speed * self.time
            wave[mask3] = self.amplitude * self.visualization_scale * np.sin(phase3)
        
        return wave
        
    def update_animation(self):
        """Update the animation for each timer tick"""
        self.time += 0.01
    
        if self.white_light:
            # Update multiple waves with different wavelengths
            for wl, curve in self.wave_curves:
                wave = self.calculate_wave(wl)
                curve.setData(self.x, wave)
        else:
            # Update single wave
            wave = self.calculate_wave(self.wavelength)
            self.wave_curve.setData(self.x, wave)
        
        # Update interference waves if enabled
        if self.show_interference:
            # Always recalculate interference waves to ensure they're animated properly
            wave1, wave2 = self.calculate_interference_waves(self.wavelength)
            self.interference_wave1.setData(self.x, wave1)
            self.interference_wave2.setData(self.x, wave2)
            
            # Ensure they're visible (in case they were toggled on but not showing)
            if not self.interference_wave1.isVisible() or not self.interference_wave2.isVisible():
                self.interference_wave1.setVisible(True)
                self.interference_wave2.setVisible(True)

    def toggle_interference(self, enabled):
        """Toggle visibility of interference waves"""
        self.show_interference = enabled
        
        if enabled:
            # Force a significant time update to ensure waves are animated
            old_time = self.time
            self.time += 2.5  # Add a significant time offset
            
            # Force recalculation with the new time
            wave1, wave2 = self.calculate_interference_waves(self.wavelength)
            
            # Update the interference waves
            self.interference_wave1.setData(self.x, wave1)
            self.interference_wave2.setData(self.x, wave2)
            
            # Make them visible
            self.interference_wave1.setVisible(True)
            self.interference_wave2.setVisible(True)
            
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
            
            # Create medium rectangles using PyQtGraph's built-in rectangle
            # Medium 1
            medium1_rect = pg.PlotCurveItem(
                x=[0, 0, self.boundary1, self.boundary1, 0],
                y=[-10, 10, 10, -10, -10],
                fillLevel=0,
                fillBrush=pg.mkBrush(self.medium1_color),
                pen=pg.mkPen(None))
            self.addItem(medium1_rect)
            self.medium_rects.append(medium1_rect)
            
            # Medium 2
            medium2_rect = pg.PlotCurveItem(
                x=[self.boundary1, self.boundary1, self.boundary2, self.boundary2, self.boundary1],
                y=[-10, 10, 10, -10, -10],
                fillLevel=0,
                fillBrush=pg.mkBrush(self.medium2_color),
                pen=pg.mkPen(None))
            self.addItem(medium2_rect)
            self.medium_rects.append(medium2_rect)
            
            # Medium 3
            medium3_rect = pg.PlotCurveItem(
                x=[self.boundary2, self.boundary2, 10, 10, self.boundary2],
                y=[-10, 10, 10, -10, -10],
                fillLevel=0,
                fillBrush=pg.mkBrush(self.medium3_color),
                pen=pg.mkPen(None))
            self.addItem(medium3_rect)
            self.medium_rects.append(medium3_rect)
            
            # Add boundary lines
            boundary1_line = pg.InfiniteLine(pos=self.boundary1, angle=90, pen=pg.mkPen('w', width=2, style=Qt.DashLine))
            boundary2_line = pg.InfiniteLine(pos=self.boundary2, angle=90, pen=pg.mkPen('w', width=2, style=Qt.DashLine))
            self.addItem(boundary1_line)
            self.addItem(boundary2_line)
            self.medium_rects.append(boundary1_line)
            self.medium_rects.append(boundary2_line)
            
            # Update medium labels
            for i, (color, name, index) in enumerate(self.medium_labels):
                self.removeItem(color)
                self.removeItem(name)
            self.medium_labels = []
            
            # Add medium labels
            # Medium 1
            medium1_color_label = pg.TextItem(f"n₁: {self.n1:.4f}", anchor=(0.5, 0), color='white')
            medium1_color_label.setPos(self.boundary1 / 2, -5)
            self.addItem(medium1_color_label)
            
            medium1_name_label = pg.TextItem(f"Medium 1", anchor=(0.5, 1), color='white')
            medium1_name_label.setPos(self.boundary1 / 2, 5)
            self.addItem(medium1_name_label)
            
            self.medium_labels.append((medium1_color_label, medium1_name_label, 1))
            
            # Medium 2
            medium2_color_label = pg.TextItem(f"n₂: {self.n2:.4f}", anchor=(0.5, 0), color='white')
            medium2_color_label.setPos(self.boundary1 + (self.boundary2 - self.boundary1) / 2, -5)
            self.addItem(medium2_color_label)
            
            medium2_name_label = pg.TextItem(f"Medium 2", anchor=(0.5, 1), color='white')
            medium2_name_label.setPos(self.boundary1 + (self.boundary2 - self.boundary1) / 2, 5)
            self.addItem(medium2_name_label)
            
            self.medium_labels.append((medium2_color_label, medium2_name_label, 2))
            
            # Medium 3
            medium3_color_label = pg.TextItem(f"n₃: {self.n3:.4f}", anchor=(0.5, 0), color='white')
            medium3_color_label.setPos(self.boundary2 + (3000 - self.boundary2) / 2, -5)
            self.addItem(medium3_color_label)
            
            medium3_name_label = pg.TextItem(f"Medium 3", anchor=(0.5, 1), color='white')
            medium3_name_label.setPos(self.boundary2 + (3000 - self.boundary2) / 2, 5)
            self.addItem(medium3_name_label)
            
            self.medium_labels.append((medium3_color_label, medium3_name_label, 3))
            
            # Update the wave curve with current parameters
            if not self.white_light:
                # Get color from wavelength
                r, g, b = wavelength_to_rgb(self.wavelength)
                self.wave_curve.setPen(pg.mkPen(QColor(r, g, b), width=4))
                
                # Update wave data
                wave = self.calculate_wave(self.wavelength)
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
            print(f"Error in update_plot: {e}")
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

    def update_wavelength(self, value):
        """Update the wavelength and colors"""
        self.wavelength = value
        
        # Update wave color based on wavelength
        if not self.white_light:
            r, g, b = wavelength_to_rgb(self.wavelength)
            wave_color = QColor(r, g, b)
            
            # Update wave curve color
            self.wave_curve.setPen(pg.mkPen(wave_color, width=4))
            
            # Update ray colors if ray mode is enabled
            if self.show_ray_mode:
                ray_color = pg.mkPen(wave_color, width=4, style=Qt.DashLine)
                self.ray_incident.setPen(ray_color)
                self.ray_refracted1.setPen(ray_color)
                self.ray_refracted2.setPen(ray_color)
        
        # Update the wave data
        self.update_plot()
        
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
        self.n1 = value
        self.update_plot()
        
    def update_n2(self, value):
        self.n2 = value
        self.update_plot()
        
    def update_n3(self, value):
        self.n3 = value
        self.update_plot()
        
    def update_medium1(self, medium_name):
        """Update medium 1 selection"""
        if hasattr(self, 'medium_presets') and medium_name in self.medium_presets:
            # Update the refractive index directly
            self.n1 = self.medium_presets[medium_name]['n']
            
            # Update medium color in the visualization
            if medium_name in self.medium_colors:
                self.medium1_region.setBrush(QBrush(QColor(*self.medium_colors[medium_name])))
                
            # Update medium label
            self.medium1_label.setText(f'{medium_name} (n₁ = {self.n1:.4f})')
            
            # Update the visualization
            self.update_plot()
        
    def update_medium2(self, medium_name):
        """Update medium 2 selection"""
        if hasattr(self, 'medium_presets') and medium_name in self.medium_presets:
            # Update the refractive index directly
            self.n2 = self.medium_presets[medium_name]['n']
            
            # Update medium color in the visualization
            if medium_name in self.medium_colors:
                self.medium2_region.setBrush(QBrush(QColor(*self.medium_colors[medium_name])))
                
            # Update medium label
            self.medium2_label.setText(f'{medium_name} (n₂ = {self.n2:.4f})')
            
            # Update the visualization
            self.update_plot()
        
    def update_medium3(self, medium_name):
        """Update medium 3 selection"""
        if hasattr(self, 'medium_presets') and medium_name in self.medium_presets:
            # Update the refractive index directly
            self.n3 = self.medium_presets[medium_name]['n']
            
            # Update medium color in the visualization
            if medium_name in self.medium_colors:
                self.medium3_region.setBrush(QBrush(QColor(*self.medium_colors[medium_name])))
                
            # Update medium label
            self.medium3_label.setText(f'{medium_name} (n₃ = {self.n3:.4f})')
            
            # Update the visualization
            self.update_plot()
        
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
            
            # Hide the original single wavelength curve when in white light mode
            self.wave_curve.setVisible(False)
        else:
            # Show the original single wavelength curve when not in white light mode
            self.wave_curve.setVisible(True)
        
        # Update visibility of wave curves
        for wl, curve in self.wave_curves:
            curve.setVisible(enabled)

    def create_white_light_curves(self):
        """Create curves for multiple wavelengths to simulate white light"""
        # Clear existing curves if any
        if hasattr(self, 'wave_curves'):
            for wl, curve in self.wave_curves:
                self.removeItem(curve)
                
        self.wave_curves = []
        
        # Create curves for different wavelengths in the visible spectrum
        for wl in self.prism_wavelengths:
            # Calculate wave for this wavelength
            wave = self.calculate_wave(wl)
            
            # Get color for this wavelength
            r, g, b = wavelength_to_rgb(wl)
            wave_color = QColor(r, g, b)
            
            # Create curve with this color
            curve = self.plot(self.x, wave, pen=pg.mkPen(wave_color, width=3))
            
            # Store wavelength and curve
            self.wave_curves.append((wl, curve))

    def update_ray_lines(self):
        """Update the ray lines to visualize refraction angles"""
        try:
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
                self.angle1_label.setText(f"θ₁: {self.angle_of_incidence}°")
                self.angle2_label.setText("θ₂: TIR")
                self.angle3_label.setText("θ₃: N/A")
                
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
                self.angle1_label.setText(f"θ₁: {self.angle_of_incidence}°")
                self.angle2_label.setText(f"θ₂: {angle_refraction1_deg:.1f}°")
                self.angle3_label.setText("θ₃: TIR")
                
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
            # Safely calculate refraction angle, avoiding arcsin errors
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
                self.angle1_label.setText(f"θ₁: {self.angle_of_incidence}°")
                self.angle2_label.setText(f"θ₂: {angle_refraction1_deg:.1f}°")
                self.angle3_label.setText(f"θ₃: {angle_refraction2_deg:.1f}°")
            except Exception as e:
                # If arcsin calculation fails, it should be TIR - handle accordingly
                print(f"Exception in calculating refraction angle: {e}")
                
                # Mark as TIR
                self.angle1_label.setText(f"θ₁: {self.angle_of_incidence}°")
                self.angle2_label.setText(f"θ₂: {angle_refraction1_deg:.1f}°")
                self.angle3_label.setText("θ₃: TIR (error)")
                
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
            self.angle1_label.setText(f"θ₁: {self.angle_of_incidence}°")
            self.angle2_label.setText("θ₂: Error")
            self.angle3_label.setText("θ₃: Error")

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
                border: 1px solid #999999;
                height: 10px;
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0.000 #380082,    /* 380nm - violet */
                    stop:0.161 #4400ff,    /* 440nm - deep blue */
                    stop:0.297 #0088ff,    /* 490nm - light blue */
                    stop:0.351 #00ff00,    /* 510nm - green */
                    stop:0.540 #ffff00,    /* 580nm - yellow */
                    stop:0.716 #ff8800,    /* 620nm - orange */
                    stop:0.716 #ff0000,    /* 645nm - red */
                    stop:1.000 #820000);   /* 750nm - deep red */
                margin: 2px 0;
            }
            QSlider::handle:horizontal {
                background: white;
                border: 1px solid #565a5e;
                width: 10px;
                margin: -4px 0;
                border-radius: 3px;
            }
        """)


class BlueSlider(QSlider):
    def __init__(self, parent=None):
        super().__init__(Qt.Horizontal, parent)
        self.setStyleSheet("""
            QSlider::groove:horizontal {
                border: 1px solid #999999;
                height: 10px;
                background: white;
                margin: 2px 0;
            }
            QSlider::sub-page:horizontal {
                background: #0088ff;
                border: 1px solid #999999;
                height: 10px;
                margin: 2px 0;
            }
            QSlider::handle:horizontal {
                background: white;
                border: 1px solid #565a5e;
                width: 10px;
                margin: -4px 0;
                border-radius: 3px;
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
            'Air': {'n': 1.0003, 'color': '#222233'},
            'Water': {'n': 1.33, 'color': '#2244AA'},
            'Glass (Crown)': {'n': 1.52, 'color': '#445566'},
            'Glass (Flint)': {'n': 1.62, 'color': '#556677'},
            'Diamond': {'n': 2.42, 'color': '#778899'},
            'Acrylic': {'n': 1.49, 'color': '#223355'},
            'Glycerine': {'n': 1.47, 'color': '#336699'},
            'Ethanol': {'n': 1.36, 'color': '#3377AA'},
            'Quartz': {'n': 1.54, 'color': '#667788'},
            'Sapphire': {'n': 1.77, 'color': '#5566AA'}
        }
        
        # Preset scenarios
        self.scenario_materials = {
            'Air → Water → Glass': ('Air', 'Water', 'Glass (Crown)'),
            'Air → Glass → Water': ('Air', 'Glass (Crown)', 'Water'),
            'Water → Air → Glass': ('Water', 'Air', 'Glass (Crown)'),
            'Air → Diamond → Glass': ('Air', 'Diamond', 'Glass (Crown)'),
            'Glass → Air → Water': ('Glass (Crown)', 'Air', 'Water')
        }
        
        # Create a color indicator for wavelength (we'll hide this)
        self.wavelength_color_indicator = QFrame()
        self.wavelength_color_indicator.setFixedSize(0, 0)  # Set to zero size to hide it
        self.wavelength_color_indicator.setFrameShape(QFrame.NoFrame)
        self.wavelength_color_indicator.setStyleSheet("background-color: transparent;")
        
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

    def setup_wave_controls(self, layout):
        # Left controls group (wave properties)
        wave_group = QGroupBox("Wavelength")
        wave_layout = QVBoxLayout()
        
        wavelength_slider_layout = QHBoxLayout()
        wavelength_slider_layout.addWidget(QLabel("Wavelength:"))
        self.wavelength_slider = ColoredSlider()
        self.wavelength_slider.setOrientation(Qt.Horizontal)
        self.wavelength_slider.setMinimum(380)
        self.wavelength_slider.setMaximum(750)
        self.wavelength_slider.setValue(550)
        self.wavelength_slider.setTickPosition(QSlider.TicksBelow)
        self.wavelength_slider.setTickInterval(50)
        self.wavelength_slider.valueChanged.connect(self.update_wavelength)
        wavelength_slider_layout.addWidget(self.wavelength_slider)
        
        # Simplified wavelength value display (no color indicator)
        wavelength_value_layout = QHBoxLayout()
        self.wavelength_value_label = QLabel("550 nm")
        wavelength_value_layout.addWidget(self.wavelength_value_label)
        wavelength_value_layout.addStretch()
        
        wave_layout.addLayout(wavelength_slider_layout)
        wave_layout.addLayout(wavelength_value_layout)
        wave_group.setLayout(wave_layout)
        layout.addWidget(wave_group)
        
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
        self.speed_slider.setValue(2)
        self.speed_value = QLabel("2")
        
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
        self.angle_slider.setValue(0)  # 0 degrees is straight line
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
        self.interference_check.setChecked(False)  # Set unchecked by default
        self.ray_mode_check = QCheckBox("Show Ray Path")
        
        mode_layout.addWidget(self.white_light_check)
        mode_layout.addWidget(self.interference_check)
        mode_layout.addWidget(self.ray_mode_check)
        wave_layout.addLayout(mode_layout)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
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
        self.medium3_combo.setCurrentText("Glass (Crown)")
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
        self.wavelength_slider.valueChanged.connect(self.update_wavelength)
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
        
    def update_wavelength(self, value):
        """Update wavelength value"""
        # Update displayed value
        nm_value = value
        self.wavelength_value_label.setText(f"{nm_value} nm")
        
        # Update widget
        self.wave_widget.update_wavelength(value)
        
        # We're not updating the color indicator anymore since it's hidden
        
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
        if hasattr(self, 'medium_presets') and medium_name in self.medium_presets:
            # Update the refractive index directly
            self.n1 = self.medium_presets[medium_name]['n']
            
            # Update medium color in the visualization
            if medium_name in self.medium_colors:
                self.medium1_region.setBrush(QBrush(QColor(*self.medium_colors[medium_name])))
                
            # Update medium label
            self.medium1_label.setText(f'{medium_name} (n₁ = {self.n1:.4f})')
            
            # Update the visualization
            self.update_plot()
        
    def update_medium2(self, medium_name):
        """Update medium 2 selection"""
        if hasattr(self, 'medium_presets') and medium_name in self.medium_presets:
            # Update the refractive index directly
            self.n2 = self.medium_presets[medium_name]['n']
            
            # Update medium color in the visualization
            if medium_name in self.medium_colors:
                self.medium2_region.setBrush(QBrush(QColor(*self.medium_colors[medium_name])))
                
            # Update medium label
            self.medium2_label.setText(f'{medium_name} (n₂ = {self.n2:.4f})')
            
            # Update the visualization
            self.update_plot()
        
    def update_medium3(self, medium_name):
        """Update medium 3 selection"""
        if hasattr(self, 'medium_presets') and medium_name in self.medium_presets:
            # Update the refractive index directly
            self.n3 = self.medium_presets[medium_name]['n']
            
            # Update medium color in the visualization
            if medium_name in self.medium_colors:
                self.medium3_region.setBrush(QBrush(QColor(*self.medium_colors[medium_name])))
                
            # Update medium label
            self.medium3_label.setText(f'{medium_name} (n₃ = {self.n3:.4f})')
            
            # Update the visualization
            self.update_plot()
        
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