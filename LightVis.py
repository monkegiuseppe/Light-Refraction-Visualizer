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
        
       
       

        
        # Add grid
        self.showGrid(x=True, y=True, alpha=0.3)
        self.getPlotItem().getAxis('left').setPen('w')
        self.getPlotItem().getAxis('bottom').setPen('w')

        # Hide the auto-range button
        self.getPlotItem().hideButtons()
        

        # Set grid line spacing
        self.getPlotItem().getAxis('left').setTicks([[(i, str(i)) for i in range(-2, 3, 1)]])
        self.getPlotItem().getAxis('bottom').setTicks([[(i, str(i)) for i in range(0, 3001, 500)]])

        # Set up the x-axis
        self.x = np.linspace(0, 3000, 3000)

        # Initial parameters
        self.wavelength = 550
        self.amplitude = 5
        self.speed = 2
        self.visualization_scale = 0.1  # Add visualization scaling factor
        self.n1 = 1.0003  # Air
        self.n2 = 1.33    # Water
        self.n3 = 1.52    # Glass (Crown)
        self.time = 0
        self.angle_of_incidence = 0

       

        # Add interference wave curves
        self.interference_wave1 = self.plot(self.x, np.zeros_like(self.x), 
                                          pen=pg.mkPen('r', width=2, style=Qt.DashLine))
        self.interference_wave2 = self.plot(self.x, np.zeros_like(self.x), 
                                          pen=pg.mkPen('g', width=2, style=Qt.DashLine))
        
        # Hide interference waves initially
        self.interference_wave1.setVisible(False)
        self.interference_wave2.setVisible(False)
        self.show_interference = False
        

        # Define boundaries
        self.boundary1 = 1000
        self.boundary2 = 2000
        
        # Add angle labels at boundaries 
        self.angle1_label = pg.TextItem("θ₁: 0°", anchor=(0.5, 0), color='w')
        self.angle1_label.setPos(self.boundary1 - 100, -1.5)

        self.angle2_label = pg.TextItem("θ₂: 0°", anchor=(0.5, 0), color='w')
        self.angle2_label.setPos(self.boundary1 + 100, -1.5)

        self.angle3_label = pg.TextItem("θ₃: 0°", anchor=(0.5, 0), color='w')
        self.angle3_label.setPos(self.boundary2 + 100, -1.5)

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
        self.wave_curve = self.plot(self.x, wave, pen=pg.mkPen('b', width=4))

    
    def calculate_interference_waves(self, wavelength):
        """Calculate the component waves that create interference"""
        wave = np.zeros_like(self.x)
        
        # Calculate incident wave continuing straight (wave1)
        angle_incidence = np.radians(self.angle_of_incidence)
        k1 = (2 * np.pi * self.n1) / wavelength
        
        # Wave 1 - incident wave continuing straight
        wave1 = np.zeros_like(self.x)
        phase1 = k1 * self.x - self.speed * self.time
        wave1 = self.amplitude * self.visualization_scale * np.sin(phase1)
        
        # Wave 2 - additional wave that creates interference
        wave2 = np.zeros_like(self.x)
        wave2 = self.calculate_wave(wavelength) - wave1
        
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
        angle_refraction1 = np.arcsin(np.sin(angle_incidence) * self.n1 / n2_wl)
        angle_refraction2 = np.arcsin(np.sin(angle_refraction1) * n2_wl / n3_wl)
        
        # Snell's law to calculate refraction angles
        angle_refraction1 = np.arcsin(np.sin(angle_incidence) / self.n2)
        angle_refraction2 = np.arcsin(np.sin(angle_refraction1) * self.n2 / self.n3)


        # Update angle labels
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
        phase2 = (k2 * (x2 * np.cos(angle_refraction1) + y[mask2] * np.sin(angle_refraction1)) + 
              k1 * self.boundary1 * np.cos(angle_incidence)) - self.speed * self.time
        wave[mask2] = self.amplitude * self.visualization_scale * np.sin(phase2)

        # Calculate wave in medium 3 with angular propagation
        mask3 = self.x > self.boundary2
        x3 = self.x[mask3] - self.boundary2
        phase3 = (k3 * (x3 * np.cos(angle_refraction2) + y[mask3] * np.sin(angle_refraction2)) +
              k1 * self.boundary1 * np.cos(angle_incidence) +
              k2 * (self.boundary2 - self.boundary1) * np.cos(angle_refraction1)) - self.speed * self.time
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
            
            if self.show_interference:
                wave1, wave2 = self.calculate_interference_waves(wl)
                self.interference_wave1.setData(self.x, wave1)
                self.interference_wave2.setData(self.x, wave2)
                
            
        else:
        # Update single wave
            wave = self.calculate_wave(self.wavelength)
            self.wave_curve.setData(self.x, wave)
        
        if self.show_interference:
            wave1, wave2 = self.calculate_interference_waves(self.wavelength)
            self.interference_wave1.setData(self.x, wave1)
            self.interference_wave2.setData(self.x, wave2)
            
        


    def toggle_interference(self, enabled):
        """Toggle visibility of interference waves"""
        self.show_interference = enabled
        self.interference_wave1.setVisible(enabled)
        self.interference_wave2.setVisible(enabled)
        
    def update_plot(self):
        """Update the plot with current parameters"""
        # Update medium labels
        self.medium1_label.setText(f"{self.current_medium1} (n₁ = {self.n1:.4f})")
        self.medium2_label.setText(f"{self.current_medium2} (n₂ = {self.n2:.4f})")
        self.medium3_label.setText(f"{self.current_medium3} (n₃ = {self.n3:.4f})")
        
        # Update medium colors
        self.medium1_region.setBrush(QBrush(QColor(*self.medium_colors[self.current_medium1])))
        self.medium2_region.setBrush(QBrush(QColor(*self.medium_colors[self.current_medium2])))
        self.medium3_region.setBrush(QBrush(QColor(*self.medium_colors[self.current_medium3])))

        
        # Update wavelength color for single wave mode
        if not self.white_light:
            color = wavelength_to_rgb(self.wavelength)
            self.wave_curve.setPen(pg.mkPen(color, width=4))

        # Update wave curve
        wave = self.calculate_wave(self.wavelength)
        self.wave_curve.setData(self.x, wave)

    def update_wavelength(self, value):
        self.wavelength = value
        self.update_plot()
        
        
        
    def update_amplitude(self, value):
        self.amplitude = value
        # Let the zoom control handle the Y range adjustment if it exists
        if hasattr(self, 'parent') and hasattr(self.parent, 'zoom_slider'):
            zoom_factor = self.parent.zoom_slider.value() / 10.0
            y_range = self.amplitude * self.visualization_scale * 1.5 / zoom_factor
            self.setYRange(-y_range, y_range)
        else:
            # Default behavior with visualization scaling
            self.setYRange(-1, 1, padding=0)

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
        self.current_medium1 = medium_name
        self.update_plot()
        
    def update_medium2(self, medium_name):
        self.current_medium2 = medium_name
        self.update_plot()
        
    def update_medium3(self, medium_name):
        self.current_medium3 = medium_name
        self.update_plot()
        
    def toggle_prism_mode(self, enabled):
        """Toggle between normal and prism simulation mode"""
        self.prism_mode = enabled
        self.update_plot()
        
    def toggle_white_light(self, enabled):
        """Toggle between single wavelength and white light (multiple wavelengths)"""
        self.white_light = enabled
        
        # Show/hide appropriate wave curves
        self.wave_curve.setVisible(not enabled)
        
        for _, curve in self.wave_curves:
            curve.setVisible(enabled)
        
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



class LightSimulationApp(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("Light Wave Refraction Simulation")
        self.setGeometry(100, 100, 1200, 800)
        
        # Set dark mode stylesheet
        self.setStyleSheet("""
        QMainWindow {
            background-color: #2b2b2b;
            font-family: 'Foto';
        }
        QWidget {
            background-color: #2b2b2b;
            color: #ffffff;
            font-family: 'Foto';
        }
        QGroupBox {
            border: 1px solid #ffffff;
            margin-top: 10px;
            padding-top: 10px;
            font-family: 'Foto';
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            subcontrol-position: top left;
            padding: 0 3px;
            color: #ffffff;
            font-family: 'Foto';
        }
        QLabel {
            color: #ffffff;
            font-family: 'Foto';
        }
        QSlider {
            background-color: transparent;
        }
        QSlider::handle {
            background-color: #ffffff;
        }
        QComboBox {
            color: #ffffff;
            background-color: #3c3c3c;
            selection-background-color: #404040;
            border: 1px solid #ffffff;
            font-family: 'Foto';
        }
        QPushButton {
            color: #ffffff;
            background-color: #3c3c3c;
            border: 1px solid #ffffff;
            padding: 5px;
            font-family: 'Foto';
        }
        QPushButton:hover {
            background-color: #404040;
        }
        QCheckBox {
            color: #ffffff;
            font-family: 'Foto';
        }
    """)

        # Medium presets (refractive indices at ~550nm wavelength)
        self.medium_presets = {
            'Air': 1.0003,
            'Water': 1.33,
            'Glass (Crown)': 1.52,
            'Glass (Flint)': 1.62,
            'Diamond': 2.42,
            'Acrylic': 1.49,
            'Glycerine': 1.47,
            'Ethanol': 1.36,
            'Quartz': 1.54,
            'Sapphire': 1.77
        }
        
        # Preset scenarios
        self.scenario_materials = {
            'Air → Water → Glass': ('Air', 'Water', 'Glass (Crown)'),
            'Air → Glass → Water': ('Air', 'Glass (Crown)', 'Water'),
            'Water → Air → Glass': ('Water', 'Air', 'Glass (Crown)'),
            'Air → Diamond → Glass': ('Air', 'Diamond', 'Glass (Crown)'),
            'Glass → Air → Water': ('Glass (Crown)', 'Air', 'Water')
        }
        
        # Create the main widget and layout
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        
        # Main vertical layout
        self.main_layout = QVBoxLayout()
        self.main_layout.setSpacing(0)
        self.main_widget.setLayout(self.main_layout)
        
        # Create the wave visualization
        self.wave_widget = WaveSimulationWidget()
        self.main_layout.addWidget(self.wave_widget, stretch=1)

        # Create and add controls
        controls_container = QWidget()
        controls_container.setMaximumHeight(200)
        controls_container.setMinimumHeight(200)
        controls_layout = QHBoxLayout()
        controls_container.setLayout(controls_layout)


        
        
        # Add wave controls
        self.setup_wave_controls(controls_layout)
        
        self.main_layout.addWidget(controls_container)
    
    def setup_wave_controls(self,controls_layout):
        """Set up the control panel for wave simulation"""
        
        
        # Left controls group (wave properties)
        wave_group = QGroupBox("Wave Properties")
        wave_layout = QVBoxLayout()
        wave_group.setLayout(wave_layout)
        controls_layout.addWidget(wave_group)
        
        # Wavelength slider
        wavelength_layout = QHBoxLayout()
        wavelength_label = QLabel("Wavelength:")
        self.wavelength_slider = ColoredSlider()  # Use the custom slider
        self.wavelength_slider.setMinimum(380)
        self.wavelength_slider.setMaximum(750)
        self.wavelength_slider.setValue(550)
        self.wavelength_value = QLabel("550 nm")
        
        wavelength_layout.addWidget(wavelength_label)
        wavelength_layout.addWidget(self.wavelength_slider)
        wavelength_layout.addWidget(self.wavelength_value)
        wave_layout.addLayout(wavelength_layout)
        
        # Amplitude slider
        amplitude_layout = QHBoxLayout()
        amplitude_label = QLabel("Amplitude:")
        self.amplitude_slider = QSlider(Qt.Horizontal)
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
        self.speed_slider = QSlider(Qt.Horizontal)
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
        self.angle_slider = QSlider(Qt.Horizontal)
        self.angle_slider.setMinimum(0)
        self.angle_slider.setMaximum(90)
        self.angle_slider.setValue(0)
        self.angle_value = QLabel("0°")
    
        angle_layout.addWidget(angle_label)
        angle_layout.addWidget(self.angle_slider)
        angle_layout.addWidget(self.angle_value)
        wave_layout.addLayout(angle_layout)


        # Add interference mode checkbox next to other mode controls
        mode_layout = QHBoxLayout()
        self.white_light_check = QCheckBox("White Light")
        self.interference_check = QCheckBox("Show Interference")

        mode_layout.addWidget(self.white_light_check)
        mode_layout.addWidget(self.interference_check)
        wave_layout.addLayout(mode_layout)
        
        
        # Middle group (medium 1)
        medium1_group = QGroupBox("Medium 1")
        medium1_layout = QVBoxLayout()
        medium1_group.setLayout(medium1_layout)
        controls_layout.addWidget(medium1_group)
        
        # Medium 1 selection
        self.medium1_combo = QComboBox()
        for medium in sorted(self.medium_presets.keys()):
            self.medium1_combo.addItem(medium)
        medium1_layout.addWidget(self.medium1_combo)
        
        # Medium 1 n slider
        n1_layout = QHBoxLayout()
        n1_label = QLabel("n₁:")
        self.n1_slider = QSlider(Qt.Horizontal)
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
        controls_layout.addWidget(medium2_group)
        
        # Medium 2 selection
        self.medium2_combo = QComboBox()
        for medium in sorted(self.medium_presets.keys()):
            self.medium2_combo.addItem(medium)
        self.medium2_combo.setCurrentText("Water")
        medium2_layout.addWidget(self.medium2_combo)
        
        # Medium 2 n slider
        n2_layout = QHBoxLayout()
        n2_label = QLabel("n₂:")
        self.n2_slider = QSlider(Qt.Horizontal)
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
        controls_layout.addWidget(medium3_group)
        
        # Medium 3 selection
        self.medium3_combo = QComboBox()
        for medium in sorted(self.medium_presets.keys()):
            self.medium3_combo.addItem(medium)
        self.medium3_combo.setCurrentText("Glass (Crown)")
        medium3_layout.addWidget(self.medium3_combo)
        
        # Medium 3 n slider
        n3_layout = QHBoxLayout()
        n3_label = QLabel("n₃:")
        self.n3_slider = QSlider(Qt.Horizontal)
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
        controls_layout.addWidget(scenario_group)
        
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

        self.n1_slider.valueChanged.connect(self.update_n1)
        self.n2_slider.valueChanged.connect(self.update_n2)
        self.n3_slider.valueChanged.connect(self.update_n3)
        
        self.medium1_combo.currentTextChanged.connect(self.update_medium1)
        self.medium2_combo.currentTextChanged.connect(self.update_medium2)
        self.medium3_combo.currentTextChanged.connect(self.update_medium3)
        
        
        self.white_light_check.stateChanged.connect(self.toggle_white_light)
        
        apply_button.clicked.connect(self.apply_scenario)

        self.interference_check.stateChanged.connect(self.toggle_interference)
        


    def toggle_interference(self, state):
        """Toggle interference visualization"""
        enabled = state == Qt.Checked
        self.wave_widget.toggle_interference(enabled)
    

    
    
    def on_tab_changed(self, index):
        """Handle tab change event"""
        # Pause/resume appropriate timers based on visible tab
        if index == 0:  # Wave tab
            self.wave_widget.timer.start(50)
        else:
            self.wave_widget.timer.stop()
            
    # --- Wave Control Event Handlers ---
    def update_wavelength(self, value):
        """Update wavelength value"""
        self.wavelength_value.setText(f"{value} nm")
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
        n = self.medium_presets[medium_name]
        self.n1_slider.setValue(int(n * 100))
        self.wave_widget.n1 = n
        self.wave_widget.update_medium1(medium_name)
        
    def update_medium2(self, medium_name):
        """Update medium 2 selection"""
        n = self.medium_presets[medium_name]
        self.n2_slider.setValue(int(n * 100))
        self.wave_widget.n2 = n
        self.wave_widget.update_medium2(medium_name)
        
    def update_medium3(self, medium_name):
        """Update medium 3 selection"""
        n = self.medium_presets[medium_name]
        self.n3_slider.setValue(int(n * 100))
        self.wave_widget.n3 = n
        self.wave_widget.update_medium3(medium_name)
    
        
    

    def update_zoom(self, value):
    #Update vertical zoom factor#
        zoom_factor = value / 10.0  # Convert to a scale where 10 = 1.0x
        self.zoom_value.setText(f"{zoom_factor:.1f}x")
    # Adjust the y-range of the plot widget
        amplitude = self.wave_widget.amplitude
        y_range = amplitude * 1.5 / zoom_factor
        self.wave_widget.setYRange(-y_range, y_range)
        
            
    def toggle_white_light(self, state):
        """Toggle white light mode"""
        enabled = state == Qt.Checked
        self.wave_widget.toggle_white_light(enabled)
        
    def apply_scenario(self):
        """Apply selected scenario preset"""
        scenario = self.scenario_combo.currentText()
        medium1, medium2, medium3 = self.scenario_materials[scenario]
        
        self.medium1_combo.setCurrentText(medium1)
        self.medium2_combo.setCurrentText(medium2)
        self.medium3_combo.setCurrentText(medium3)

        n1 = self.medium_presets[medium1]
        n2 = self.medium_presets[medium2]
        n3 = self.medium_presets[medium3]

    def update_angle(self, value):
        """Update angle of incidence"""
        self.angle_value.setText(f"{value}°")
        self.wave_widget.angle_of_incidence = value
        self.wave_widget.update_plot()


# Run the application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    mainWindow = LightSimulationApp()
    mainWindow.show()
    sys.exit(app.exec_())