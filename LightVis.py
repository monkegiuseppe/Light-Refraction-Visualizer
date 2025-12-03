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
            self.setToolTip("Play Animation")
            self.setText("▶")
        else:
            self.setToolTip("Pause Animation")
            self.setText("||")

        self.setStyleSheet(self.styleSheet() + """
            QPushButton {
                color: white;
                font-size: 18px;
                font-weight: bold;
            }
        """)

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

def wavelength_to_rgb(wavelength):
    """Convert wavelength (in nm) to RGB color values"""
    if wavelength < 380:
        wavelength = 380
    elif wavelength > 750:
        wavelength = 750
    
    if 380 <= wavelength < 440:
        r = ((440 - wavelength) / (440 - 380)) * 0.8
        g = 0.0
        b = 1.0
    elif 440 <= wavelength < 490:
        r = 0.0
        g = ((wavelength - 440) / (490 - 440))
        b = 1.0
    elif 490 <= wavelength < 510:
        r = 0.0
        g = 1.0
        b = ((510 - wavelength) / (510 - 490))
    elif 510 <= wavelength < 580:
        r = ((wavelength - 510) / (580 - 510))
        g = 1.0
        b = 0.0
    elif 580 <= wavelength < 645:
        r = 1.0
        g = ((645 - wavelength) / (645 - 580))
        b = 0.0
    else: 
        r = 1.0
        g = 0.0
        b = 0.0
    
    r = int(r * 255)
    g = int(g * 255)
    b = int(b * 255)
    
    return r, g, b

def frequency_to_rgb(frequency_THz):
    """Convert frequency in THz to RGB color"""
    wavelength_nm = 299792.458 / frequency_THz
    return wavelength_to_rgb(wavelength_nm)

class WaveSimulationWidget(pg.PlotWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.angle_of_incidence = 0
        self.time = 0.0
        self.visualization_scale = 0.1
        self.medium1_color = '#22222259'  
        self.medium2_color = '#1E90FF59'  
        self.medium3_color = '#88DDFF80'  
        
        self.medium1_name = 'Air'
        self.medium2_name = 'Water'
        self.medium3_name = 'Glass (Crown)'

        self.medium_presets = {}
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
        self.show_reflections = False

        self.prism_frequencies = np.linspace(400, 790, 10)
        self.prism_wavelengths = [299792.458 / freq for freq in self.prism_frequencies]
        self.n1 = 1.0  
        self.n2 = 1.33  
        self.n3 = 1.5  
        
        self.boundary1 = 1000  # First boundary position
        self.boundary2 = 2000  # Second boundary position
        
        self.setBackground(None)  
        self.getPlotItem().setTitle("")  
        self.getPlotItem().showAxis('left', True)
        self.getPlotItem().getAxis('left').setLabel('Electric Field', units='V/m') 
        self.getPlotItem().hideAxis('bottom')  

        grid_pen = pg.mkPen(color=(255, 255, 255, 80), width=0.5, style=Qt.DotLine)
        
        self.getPlotItem().getAxis('left').setPen(grid_pen)
        self.getPlotItem().getAxis('bottom').setPen(grid_pen)
        
        
        self.setMouseEnabled(x=False, y=False)
        self.setMenuEnabled(False)
        
        self.setXRange(0, 3000, padding=0)
        self.setYRange(-2, 2, padding=0)
        
        self.medium_rects = []
        self.medium_labels = []
        
        self.x = np.linspace(0, 3000, 6000)
 
        
        self.wave_curves = []   
        

        self.reflected_wave1 = self.plot(self.x, np.zeros_like(self.x), 
                                        pen=pg.mkPen('#00FFFF', width=2, style=Qt.DashLine))
        self.reflected_wave2 = self.plot(self.x, np.zeros_like(self.x), 
                                        pen=pg.mkPen('#FFFF00', width=2, style=Qt.DashLine))
        
        self.reflected_wave1.setVisible(False)
        self.reflected_wave2.setVisible(False)

        self.setAntialiasing(True)

        self.grid_lines_x = []
        self.grid_lines_y = []
        
        for y in np.arange(-2, 3, 1):
            line = pg.InfiniteLine(pos=y, angle=0, pen=pg.mkPen(color=(200, 200, 200, 150), width=0.5, style=Qt.DotLine))
            self.addItem(line)
            self.grid_lines_y.append(line)
        
        
        for x in np.arange(0, 3001, 500):
            line = pg.InfiniteLine(pos=x, angle=90, pen=pg.mkPen(color=(200, 200, 200, 150), width=0.5, style=Qt.DotLine))
            self.addItem(line)
            self.grid_lines_x.append(line)

        self.getPlotItem().setContentsMargins(0, 0, 0, 0)
        self.getPlotItem().layout.setContentsMargins(0, 0, 0, 0)
        self.getPlotItem().vb.setContentsMargins(0, 0, 0, 0)
        self.getPlotItem().getViewBox().setBackgroundColor(None)
        self.getPlotItem().getViewBox().setBorder(None)
        
        for axis in ['left', 'bottom', 'top', 'right']:
            if self.getPlotItem().getAxis(axis):
                self.getPlotItem().getAxis(axis).setStyle(showValues=True, tickLength=5)
                if axis in ['top', 'right']:
                    self.getPlotItem().showAxis(axis, False)

        self.getPlotItem().getAxis('left').setTicks([[(i, f"{int(i*50)} V/m") for i in range(-2, 3, 1)]])
        self.getPlotItem().getAxis('bottom').setTicks([[(i, str(i)) for i in range(0, 3001, 500)]])

        self.x = np.linspace(0, 3000, 3000)
        self.boundary1_line = pg.InfiniteLine(pos=self.boundary1, angle=90, pen=pg.mkPen('w', width=2, style=Qt.DashLine))
        self.boundary2_line = pg.InfiniteLine(pos=self.boundary2, angle=90, pen=pg.mkPen('w', width=2, style=Qt.DashLine))
        self.addItem(self.boundary1_line)
        self.addItem(self.boundary2_line)
        
        self.medium_rects = []  # Initialize as empty list
        self.medium_rects.append(self.boundary1_line)
        self.medium_rects.append(self.boundary2_line)
        

        self.wavelength = 550
        self.amplitude = 5
        self.speed = 4
        self.n1 = 1.0003  # Air
        self.n2 = 1.33    # Water
        self.n3 = 1.52    # Glass (Crown)
        self.time = 0
        self.angle_of_incidence = 0  # Default angle is 0 (straight line)
        self.ray_target_y = 0  # Y-coordinate where ray hits first boundary
        
        self.show_interference = False
        self.show_ray_mode = False
        self.white_light = False
        self.prism_mode = False

        
        label_html_style = """<div style="font-family: Arial; font-size: 16pt; font-weight: bold; color: white;">"""
        self.angle1_label = pg.TextItem(html=f"{label_html_style}θ₁: 0°</div>", anchor=(0.5, 0.5))
        self.angle2_label = pg.TextItem(html=f"{label_html_style}θ₂: 0°</div>", anchor=(0.5, 0.5))
        self.angle3_label = pg.TextItem(html=f"{label_html_style}θ₃: 0°</div>", anchor=(0.5, 0.5))
        self.addItem(self.angle1_label)
        self.addItem(self.angle2_label)
        self.addItem(self.angle3_label)
        
        self.angle1_label.setPos(self.boundary1/2, -1.5)
        self.angle2_label.setPos(self.boundary1 + (self.boundary2-self.boundary1)/2, -1.5)
        self.angle3_label.setPos(self.boundary2 + (3000-self.boundary2)/2, -1.5)

        
        self.interference_wave1 = self.plot(self.x, np.zeros_like(self.x), 
                                          pen=pg.mkPen('r', width=2, style=Qt.DashLine))
        self.interference_wave2 = self.plot(self.x, np.zeros_like(self.x), 
                                          pen=pg.mkPen('g', width=2, style=Qt.DashLine))
        
        wave = self.calculate_wave(self.wavelength)
        r, g, b = wavelength_to_rgb(self.wavelength)
        wave_color = QColor(r, g, b)
        self.wave_curve = self.plot(self.x, wave, pen=pg.mkPen(wave_color, width=4))

        
        reduced_x, wave1, wave2 = self.calculate_interference_waves(self.wavelength)
        self.interference_wave1.setData(reduced_x, wave1)
        self.interference_wave2.setData(reduced_x, wave2)
        
        self.interference_wave1.setVisible(False)
        self.interference_wave2.setVisible(False)
        self.update_plot()
        
        
        self.ray_incident = pg.PlotDataItem([], [], pen=pg.mkPen('#ffffff', width=4, style=Qt.DashLine))
        self.ray_refracted1 = pg.PlotDataItem([], [], pen=pg.mkPen('#00aaff', width=4, style=Qt.DashLine))
        self.ray_refracted2 = pg.PlotDataItem([], [], pen=pg.mkPen('#22ff22', width=4, style=Qt.DashLine))
        
        self.reflection_marker1 = pg.ScatterPlotItem(size=12, brush='y', symbol='o')
        self.reflection_marker2 = pg.ScatterPlotItem(size=12, brush='y', symbol='o')
        
        
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

 
        
        self.last_update_time = time.time()
        self.frame_time = 1/60  # Target 60 FPS
        self.skip_frames = 0
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_animation)
        self.timer.start(16)  # 16ms interval (60 fps)

        self.update_wavelength_scale()
        

    def calculate_interference_waves(self, wavelength):
        """Calculate the component waves that create interference"""
        reduced_x = self.x[::3]
        if not self.show_interference:
            return reduced_x, np.zeros_like(reduced_x), np.zeros_like(reduced_x)
        
        actual_wave = self.calculate_wave(self.frequency)
        reduced_actual = actual_wave[::3]
        
        
        vacuum_wavelength = 299792.458 / self.frequency
        
        k_vacuum = 2 * np.pi / vacuum_wavelength
        
        omega_vacuum = 2 * np.pi * self.frequency * (self.speed / 50) / 100
        phase_vacuum = omega_vacuum * self.time
        
        
        reduced_wave1 = self.amplitude * self.visualization_scale * 1.5 * np.sin(k_vacuum * reduced_x - phase_vacuum)
        
        
        reduced_wave2 = reduced_actual - reduced_wave1
        
        return reduced_x, reduced_wave1, reduced_wave2
        
    def calculate_superposition_interference(self):
        """Calculate interference waves based on the superposition of all wavelengths"""

        cache_key = (self.time, self.n1, self.n2, self.n3)
        if hasattr(self, '_interference_cache') and cache_key in self._interference_cache:
            return self._interference_cache[cache_key]
        reduced_x = self.x[::3]
        
        superposition = np.zeros_like(reduced_x)
        
        for curve, freq in self.wave_curves:
            wave = self.calculate_wave(freq)
            reduced_wave = wave[::3]  # Take every 3rd point
            superposition += reduced_wave
            
        if len(self.wave_curves) > 0:
            superposition = superposition / len(self.wave_curves)
        
        
        vacuum_superposition = np.zeros_like(reduced_x)
        
        
        for curve, freq in self.wave_curves:
            
            vacuum_wavelength = 299792.458 / freq
            
            
            k_vacuum = 2 * np.pi / vacuum_wavelength
            
            
            omega_vacuum = 2 * np.pi * freq * (self.speed / 50) / 100
            phase_vacuum = omega_vacuum * self.time
            
            
            vacuum_wave = self.amplitude * self.visualization_scale * 1.5 * np.sin(k_vacuum * reduced_x - phase_vacuum)
            
            
            vacuum_superposition += vacuum_wave
        
        
        if len(self.wave_curves) > 0:
            vacuum_superposition = vacuum_superposition / len(self.wave_curves)
        
        wave1 = vacuum_superposition
        
        wave2 = superposition - vacuum_superposition
            
        
        if not hasattr(self, '_interference_cache'):
            self._interference_cache = {}
        self._interference_cache[cache_key] = (reduced_x, wave1, wave2)
        
        
        if len(self._interference_cache) > 50:
            oldest_keys = list(self._interference_cache.keys())[:10]
            for key in oldest_keys:
                self._interference_cache.pop(key)

        return reduced_x, wave1, wave2

    def calculate_wave(self, frequency, additional_phase=0):
        """Calculate wave values for the given frequency"""
        
        wavelength_nm = 299792.458 / frequency
        
        cache_key = (frequency, self.time, self.n1, self.n2, self.n3, 
                     self.amplitude, self.angle_of_incidence, additional_phase, self.show_reflections)
        
        
        if hasattr(self, '_wave_calc_cache') and cache_key in self._wave_calc_cache:
            return self._wave_calc_cache[cache_key]
        
        
        if not hasattr(self, '_wave_calc_cache'):
            self._wave_calc_cache = {}
        
        
        wavelength_m1 = wavelength_nm / self.n1
        wavelength_m2 = wavelength_nm / self.n2
        wavelength_m3 = wavelength_nm / self.n3
        
        
        k1 = 2 * np.pi / wavelength_m1
        k2 = 2 * np.pi / wavelength_m2
        k3 = 2 * np.pi / wavelength_m3
        
        
        omega = 2 * np.pi * frequency * (self.speed / 50) / 100
        
        phase = omega * self.time + additional_phase
        
        y = np.zeros_like(self.x)
        
        
        scaled_amplitude = self.amplitude * self.visualization_scale * 1.5 

        
        angle_rad = np.radians(self.angle_of_incidence)
        coeffs = self.calculate_reflection_coefficients(angle_rad)

        A1 = scaled_amplitude 
        A2 = A1 * coeffs['A_t1'] 
        A3 = A2 * coeffs['A_t2'] 

        mask1 = self.x < self.boundary1
        mask2 = (self.x >= self.boundary1) & (self.x < self.boundary2)
        mask3 = self.x >= self.boundary2
        
        
        phase_shift2 = (k1 - k2) * self.boundary1
        phase_shift3 = (k1 - k2) * self.boundary1 + (k2 - k3) * self.boundary2
        
        
        y[mask1] = A1 * np.sin(k1 * self.x[mask1] - phase)
        y[mask2] = A2 * np.sin(k2 * self.x[mask2] - phase + phase_shift2)
        y[mask3] = A3 * np.sin(k3 * self.x[mask3] - phase + phase_shift3)
        
        
        self._wave_calc_cache[cache_key] = y
        
        
        if len(self._wave_calc_cache) > 100:
            keys = list(self._wave_calc_cache.keys())
            oldest_keys = keys[:10]
            for key in oldest_keys:
                self._wave_calc_cache.pop(key)
        
        return y

   

    def calculate_reflected_waves(self):
        """Calculate reflected waves at each boundary"""
        if not self.show_reflections:
            return
            
        
        angle_rad = np.radians(self.angle_of_incidence)
        coeffs = self.calculate_reflection_coefficients(angle_rad)
        
        
        wavelength_nm = 299792.458 / self.frequency
        wavelength_m1 = wavelength_nm / self.n1
        wavelength_m2 = wavelength_nm / self.n2
        
        
        k1 = 2 * np.pi / wavelength_m1
        k2 = 2 * np.pi / wavelength_m2
        
        
        omega = 2 * np.pi * self.frequency * (self.speed / 50) / 100
        phase = omega * self.time
        
        
        scaled_amplitude = self.amplitude * self.visualization_scale * 1.5
        
        
        reflected1 = np.zeros_like(self.x)
        mask1 = self.x < self.boundary1
        
        
        A_r1 = scaled_amplitude * abs(coeffs['r1'])
        
        
        phase_shift_r1 = np.pi if coeffs['r1'] < 0 else 0
        
        
        reflected1[mask1] = A_r1 * np.sin(k1 * (2 * self.boundary1 - self.x[mask1]) - phase + phase_shift_r1)
        
        
        reflected2 = np.zeros_like(self.x)
        mask2 = (self.x >= self.boundary1) & (self.x < self.boundary2)
        
        
        A_r2 = scaled_amplitude * coeffs['A_t1'] * abs(coeffs['r2'])
        
        
        phase_shift_r2 = np.pi if coeffs['r2'] < 0 else 0
        
        
        phase_shift2 = (k1 - k2) * self.boundary1
        
        
        reflected2[mask2] = A_r2 * np.sin(k2 * (2 * self.boundary2 - self.x[mask2]) - phase + phase_shift2 + phase_shift_r2)
        
        
        self.reflected_wave1.setData(self.x, reflected1)
        self.reflected_wave2.setData(self.x, reflected2)

    def toggle_main_wave(self, enabled):
        """Toggle visibility of the main wave"""
        if not self.white_light:
            self.wave_curve.setVisible(enabled)
        else:
            
            if self.superposition_enabled:
                if self.superposition_wave:
                    self.superposition_wave.setVisible(enabled)
            else:
                
                for curve_item in self.wave_curves:
                    if isinstance(curve_item, tuple) and len(curve_item) == 2:
                        curve, _ = curve_item
                        if hasattr(curve, 'setVisible'):
                            curve.setVisible(enabled)
        self.update_plot()

    def toggle_reflections(self, enabled):
        """Toggle visibility of reflection waves"""
        self.show_reflections = enabled
        self.reflected_wave1.setVisible(enabled)
        self.reflected_wave2.setVisible(enabled)
        self.update_plot()

    def calculate_energy_density(self, amplitude, refractive_index):
        """
        Calculate energy density in each medium
        Energy density = (1/2) × ε₀ × n × E²
        Where E is the amplitude (V/m)
        """
        
        epsilon_0 = 8.85e-12
        
        
        e_field = amplitude * 100
        
        
        energy_density = 0.5 * epsilon_0 * refractive_index * (e_field ** 2)
        
        return energy_density

    def calculate_reflection_coefficients(self, angle_incidence_rad):
        """Calculate reflection and transmission coefficients"""
        
        n1 = self.n1
        n2 = self.n2
        
        
        sin_theta_t = (n1 / n2) * np.sin(angle_incidence_rad)
        
        if abs(sin_theta_t) >= 1.0:
            
            r1 = 1.0
            t1 = 0.0
            theta_t1 = 0.0
        else:
            theta_t1 = np.arcsin(sin_theta_t)
            
            r1 = ((n1 * np.cos(angle_incidence_rad) - n2 * np.cos(theta_t1)) / 
                 (n1 * np.cos(angle_incidence_rad) + n2 * np.cos(theta_t1)))
            
            t1 = (2 * n1 * np.cos(angle_incidence_rad)) / (n1 * np.cos(angle_incidence_rad) + n2 * np.cos(theta_t1))
        
        
        n2 = self.n2
        n3 = self.n3
        
        
        angle_incidence2_rad = theta_t1
        
        
        sin_theta_t2 = (n2 / n3) * np.sin(angle_incidence2_rad)
        
        if abs(sin_theta_t2) >= 1.0:
            
            r2 = 1.0
            t2 = 0.0
            theta_t2 = 0.0
        else:
            theta_t2 = np.arcsin(sin_theta_t2)
            
            r2 = ((n2 * np.cos(angle_incidence2_rad) - n3 * np.cos(theta_t2)) / 
                 (n2 * np.cos(angle_incidence2_rad) + n3 * np.cos(theta_t2)))
            
            t2 = (2 * n2 * np.cos(angle_incidence2_rad)) / (n2 * np.cos(angle_incidence2_rad) + n3 * np.cos(theta_t2))
        
        
        if abs(sin_theta_t) < 1.0:
            A_t1 = t1 * np.sqrt((n2 * np.cos(theta_t1)) / (n1 * np.cos(angle_incidence_rad)))
        else:
            A_t1 = 0.0
            
        if abs(sin_theta_t2) < 1.0:
            A_t2 = t2 * np.sqrt((n3 * np.cos(theta_t2)) / (n2 * np.cos(angle_incidence2_rad)))
        else:
            A_t2 = 0.0
            
        
        A_final = A_t1 * A_t2
        
        return {
            'r1': r1,
            'r2': r2,
            't1': t1,
            't2': t2,
            'A_t1': A_t1,
            'A_t2': A_t2
        }

    def update_animation(self):
        """Update the animation for each timer tick"""
        
        current_time = time.time()
        elapsed = current_time - self.last_update_time

        if elapsed > 0.033 and self.skip_frames < 2:  
            self.skip_frames += 1
            return

        self.skip_frames = 0
        self.last_update_time = current_time

        if not self.paused: 
            self.time += 0.01
        
        
        needs_update = False
        
        
        if not hasattr(self, '_prev_state'):
            self._prev_state = {
                'time': self.time - 1,  
                'white_light': self.white_light,
                'show_interference': self.show_interference,
                'show_reflections': self.show_reflections,
                'frequency': self.frequency,
                'n1': self.n1,
                'n2': self.n2,
                'n3': self.n3,
                'superposition_enabled': self.superposition_enabled
            }
            needs_update = True
        
        
        if (abs(self._prev_state['time'] - self.time) > 0.009 or
            self._prev_state['white_light'] != self.white_light or
            self._prev_state['show_interference'] != self.show_interference or
            self._prev_state['show_reflections'] != self.show_reflections or
            self._prev_state['frequency'] != self.frequency or
            self._prev_state['superposition_enabled'] != self.superposition_enabled or
            abs(self._prev_state['n1'] - self.n1) > 0.0001 or
            abs(self._prev_state['n2'] - self.n2) > 0.0001 or
            abs(self._prev_state['n3'] - self.n3) > 0.0001):
            needs_update = True
        
        if needs_update:
            
            if self.white_light:
                
                if (self._prev_state['frequency'] != self.frequency or
                    abs(self._prev_state['n1'] - self.n1) > 0.0001 or
                    abs(self._prev_state['n2'] - self.n2) > 0.0001 or
                    abs(self._prev_state['n3'] - self.n3) > 0.0001):
                    self.update_wavelength_scale()
                
                if self.superposition_enabled:
                    self.update_superposition_wave()
                else:
                    visible_components = False
                    for curve_item in self.wave_curves:
                        if isinstance(curve_item, tuple) and len(curve_item) == 2:
                            curve, _ = curve_item
                            if hasattr(curve, 'setVisible') and curve.isVisible():
                                visible_components = True
                                break
                    
                    if visible_components:
                        
                        self.setUpdatesEnabled(False)  
                        
                        reduced_x = self.x[::20]
                        
                        for curve, freq in self.wave_curves:
                            if curve.isVisible():
                                wave = self.calculate_wave(freq)
                                reduced_wave = wave[::20]  
                                curve.setData(reduced_x, reduced_wave)
                        
                        self.setUpdatesEnabled(True)  
            else:
                
                wave = self.calculate_wave(self.frequency)
                self.wave_curve.setData(self.x, wave)
                
                
                if (self._prev_state['frequency'] != self.frequency or
                    abs(self._prev_state['n1'] - self.n1) > 0.0001 or
                    abs(self._prev_state['n2'] - self.n2) > 0.0001 or
                    abs(self._prev_state['n3'] - self.n3) > 0.0001):
                    self.update_wavelength_scale()
            
            
            if self.show_interference and needs_update:
                if self.white_light:
                    
                    reduced_x, wave1, wave2 = self.calculate_superposition_interference()
                else:
                    reduced_x, wave1, wave2 = self.calculate_interference_waves(self.wavelength)
                
                
                self.interference_wave1.setData(reduced_x, wave1)
                self.interference_wave2.setData(reduced_x, wave2)
            
             
            if self.show_reflections:
                self.calculate_reflected_waves()

            
            self._prev_state = {
                'time': self.time,
                'white_light': self.white_light,
                'show_interference': self.show_interference,
                'show_reflections': self.show_reflections,
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

        
        self.interference_wave1.setVisible(enabled)
        self.interference_wave2.setVisible(enabled)
        
        if enabled:
            
            self.interference_wave1.setPen(pg.mkPen('r', width=2, style=Qt.DashLine))
            self.interference_wave2.setPen(pg.mkPen('g', width=2, style=Qt.DashLine))
            
            
            old_time = self.time
            self.time += 0.1  
            
            
            if self.white_light:
                reduced_x, wave1, wave2 = self.calculate_superposition_interference()
            else:
                
                reduced_x, wave1, wave2 = self.calculate_interference_waves(self.wavelength)
            
            
            self.interference_wave1.setData(reduced_x, wave1)
            self.interference_wave2.setData(reduced_x, wave2)
            
            self.time = old_time + 0.01
        else: 
            self.interference_wave1.setVisible(False)
            self.interference_wave2.setVisible(False)

    def toggle_superposition(self, enabled):
        """Toggle between showing individual waves or a single superposition wave"""
        if self.white_light:
            
            self.superposition_enabled = not enabled
            
            
            if self.superposition_wave is None:
                self.superposition_wave = self.plot(self.x, np.zeros_like(self.x), 
                                                   pen=pg.mkPen('w', width=4))
            
            
            if self.superposition_wave:
                self.superposition_wave.setVisible(not enabled)
            
            
            if enabled:
                
                self.create_white_light_curves()
                
                
                for curve_item in self.wave_curves:
                    if isinstance(curve_item, tuple) and len(curve_item) == 2:
                        curve, _ = curve_item
                        if hasattr(curve, 'setVisible'):
                            curve.setVisible(True)
            else:
                
                for curve_item in self.wave_curves:
                    if isinstance(curve_item, tuple) and len(curve_item) == 2:
                        curve, _ = curve_item
                        if hasattr(curve, 'setVisible'):
                            curve.setVisible(False)
                
                self.update_superposition_wave()
        else:
            
            self.superposition_enabled = enabled
            
            
            if enabled and self.superposition_wave is None:
                self.superposition_wave = self.plot(self.x, np.zeros_like(self.x), 
                                                   pen=pg.mkPen('w', width=4))
            
            
            if self.superposition_wave:
                self.superposition_wave.setVisible(enabled)
            
            
            if enabled:
                self.update_superposition_wave()

    def update_plot(self):
        """Update the plot with current parameters"""
        try:
            self.setUpdatesEnabled(False)
            
            y_range = self.getViewBox().viewRange()[1]
            
            if not hasattr(self, '_plot_initialized'):
                self._create_initial_plot_elements()
                self._plot_initialized = True
            
            self._update_medium_rectangles()
            
            self._update_medium_labels() 
            
            self.update_energy_display()

            if not self.white_light:
                
                r, g, b = wavelength_to_rgb(299792.458 / self.frequency)
                wave_color = QColor(r, g, b)
                self.wave_curve.setPen(pg.mkPen(wave_color, width=4))
                
                
                wave = self.calculate_wave(self.frequency)
                self.wave_curve.setData(self.x, wave)
            else:
                
                all_waves = []
                for curve, freq in self.wave_curves:
                    if curve.isVisible():  
                        wave = self.calculate_wave(freq)
                        all_waves.append((curve, wave))
                
                
                for curve, wave in all_waves:
                    curve.setData(self.x, wave)
            
            if self.show_ray_mode:
                self.update_ray_lines()
                
            self.setYRange(y_range[0], y_range[1], padding=0)
            
            self.setUpdatesEnabled(True)
            
        except Exception as e:
            print(f"Error in update_plot: {str(e)}")
            self.setUpdatesEnabled(True) 
            return
            
    def _create_initial_plot_elements(self):
        """Create initial plot elements that will be reused"""
        
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
        
        self.medium1_rect.setCurves(
            pg.PlotCurveItem([0, self.boundary1], [2, 2]),
            pg.PlotCurveItem([0, self.boundary1], [-2, -2]))
        self.medium1_rect.setBrush(pg.mkBrush(self.medium1_color))
        
        
        self.medium2_rect.setCurves(
            pg.PlotCurveItem([self.boundary1, self.boundary2], [2, 2]),
            pg.PlotCurveItem([self.boundary1, self.boundary2], [-2, -2]))
        self.medium2_rect.setBrush(pg.mkBrush(self.medium2_color))
        
        self.medium3_rect.setCurves(
            pg.PlotCurveItem([self.boundary2, 3000], [2, 2]),
            pg.PlotCurveItem([self.boundary2, 3000], [-2, -2]))
        self.medium3_rect.setBrush(pg.mkBrush(self.medium3_color))
        
        
        self.boundary1_line.setValue(self.boundary1)
        self.boundary2_line.setValue(self.boundary2)
        
    def _update_medium_labels(self):
        
        self.medium1_color_label.setText(f"n₁: {self.n1:.4f}")
        self.medium1_color_label.setPos(self.boundary1 / 2, -1.2)
        self.medium1_name_label.setText(f"{self.medium1_name}")
        self.medium1_name_label.setPos(self.boundary1 / 2, 1.8)
        
        self.medium2_color_label.setText(f"n₂: {self.n2:.4f}")
        self.medium2_color_label.setPos(self.boundary1 + (self.boundary2 - self.boundary1) / 2, -1.2)
        self.medium2_name_label.setText(f"{self.medium2_name}")
        self.medium2_name_label.setPos(self.boundary1 + (self.boundary2 - self.boundary1) / 2, 1.8)
        
        
        self.medium3_color_label.setText(f"n₃: {self.n3:.4f}")
        self.medium3_color_label.setPos(self.boundary2 + (3000 - self.boundary2) / 2, -1.2)
        self.medium3_name_label.setText(f"{self.medium3_name}")
        self.medium3_name_label.setPos(self.boundary2 + (3000 - self.boundary2) / 2, 1.8)

    def toggle_ray_mode(self, enabled):
        self.show_ray_mode = enabled
        self.ray_incident.setVisible(enabled)
        self.ray_refracted1.setVisible(enabled)
        self.ray_refracted2.setVisible(enabled)
        
        
        self.reflection_marker1.setVisible(False)
        self.reflection_marker2.setVisible(False)
        
        if enabled:
            self.update_ray_lines()


    def update_wavelength_scale(self):
        scale_key = (self.frequency, self.n1, self.n2, self.n3)
        
        if hasattr(self, '_scale_key') and self._scale_key == scale_key:
            return
        
        
        self._scale_key = scale_key

        
        vacuum_wavelength_nm = 299792.458 / self.frequency
        wavelength_m1 = vacuum_wavelength_nm / self.n1
        wavelength_m2 = vacuum_wavelength_nm / self.n2
        wavelength_m3 = vacuum_wavelength_nm / self.n3
        
        speed_m1_formatted = f"{1/self.n1:.2f}c"
        speed_m2_formatted = f"{1/self.n2:.2f}c"
        speed_m3_formatted = f"{1/self.n3:.2f}c"
        
        
        if not hasattr(self, 'grid_value_labels') or len(self.grid_value_labels) == 0:
            
            self.grid_value_labels = []
            
            for x in np.arange(0, 3001, 500):
                physical_nm = x / (500/500)
                label = pg.TextItem(f"{int(physical_nm)} nm", anchor=(0, 0.5), color='white')
                label.setPos(x+10, -1.9)
                self.addItem(label)
                self.grid_value_labels.append(label)
                
            wavelength_info = pg.TextItem("", anchor=(0, 0), color='yellow')
            wavelength_info.setPos(50, 1.5)
            self.addItem(wavelength_info)
            self.grid_value_labels.append(wavelength_info)
            
            
            for i, pos in enumerate([
                self.boundary1/2,
                self.boundary1 + (self.boundary2-self.boundary1)/2,
                self.boundary2 + (3000-self.boundary2)/2
            ]):
                label = pg.TextItem("", anchor=(0.5, 0), color='yellow')
                label.setPos(pos, 1.5)
                self.addItem(label)
                self.grid_value_labels.append(label)
        
        
        wavelength_idx = len(self.grid_value_labels) - 4  
        self.grid_value_labels[wavelength_idx].setText(f"λ₀ = {int(vacuum_wavelength_nm)} nm")
        
        
        self.grid_value_labels[wavelength_idx+1].setText(f"λ₁ = {int(wavelength_m1)} nm | v₁ = {speed_m1_formatted}")
        self.grid_value_labels[wavelength_idx+2].setText(f"λ₂ = {int(wavelength_m2)} nm | v₂ = {speed_m2_formatted}")
        self.grid_value_labels[wavelength_idx+3].setText(f"λ₃ = {int(wavelength_m3)} nm | v₃ = {speed_m3_formatted}")

    def update_energy_display(self):
        
        angle_rad = np.radians(self.angle_of_incidence)
        coeffs = self.calculate_reflection_coefficients(angle_rad)
        
        initial_energy = self.calculate_energy_density(self.amplitude, self.n1)
        
        
        transmitted_energy1 = initial_energy * (coeffs['A_t1'] ** 2)
        
        
        transmitted_energy2 = transmitted_energy1 * (coeffs['A_t2'] ** 2)
        
        
        reflected_energy1 = initial_energy * (coeffs['r1'] ** 2)
        
        
        reflected_energy2 = transmitted_energy1 * (coeffs['r2'] ** 2)
        
        
        if initial_energy > 0:
            transmitted_percent1 = (transmitted_energy1 / initial_energy) * 100
            transmitted_percent2 = (transmitted_energy2 / initial_energy) * 100
            reflected_percent1 = (reflected_energy1 / initial_energy) * 100
            reflected_percent2 = (reflected_energy2 / initial_energy) * 100
        else:
            transmitted_percent1 = transmitted_percent2 = reflected_percent1 = reflected_percent2 = 0
        
        
        energy_html = f"""
        <div style='background-color: rgba(0,0,0,0.8); padding: 12px; border-radius: 8px; border: 1px solid #3498db;'>
            <div style='color: #3498db; font-size: 16px; font-weight: bold; margin-bottom: 5px;'>Energy Analysis</div>
            <div style='color: white; font-size: 14px;'>
                <span style='color: #2ecc71;'>Medium 1 → 2:</span> {transmitted_percent1:.1f}% transmitted, <span style='color: #e74c3c;'>{reflected_percent1:.1f}% reflected</span><br>
                <span style='color: #2ecc71;'>Medium 2 → 3:</span> {transmitted_percent2:.1f}% transmitted, <span style='color: #e74c3c;'>{reflected_percent2:.1f}% reflected</span>
            </div>
        </div>
        """
        
        
        if hasattr(self, 'external_energy_label') and self.external_energy_label is not None:
            self.external_energy_label.setText(energy_html)
            
            if hasattr(self, 'energy_label') and self.energy_label is not None:
                self.removeItem(self.energy_label)
                self.energy_label = None

    def update_wavelength(self, value):
        
        self.frequency = value
        self.wavelength = 299792.458 / value


        if not self.white_light:
            r, g, b = frequency_to_rgb(value)
            wave_color = QColor(r, g, b)
            
            self.wave_curve.setPen(pg.mkPen(wave_color, width=4))
            
            
            if self.show_ray_mode:
                ray_color = pg.mkPen(wave_color, width=4, style=Qt.DashLine)
                self.ray_incident.setPen(ray_color)
                self.ray_refracted1.setPen(ray_color)
                self.ray_refracted2.setPen(ray_color)
        
        self.update_wavelength_scale()
        
        if self.show_interference:
            reduced_x, wave1, wave2 = self.calculate_interference_waves(self.wavelength)
            self.interference_wave1.setData(reduced_x, wave1)
            self.interference_wave2.setData(reduced_x, wave2)
        
        
        wave = self.calculate_wave(self.frequency)
        self.wave_curve.setData(self.x, wave)

        if self.show_ray_mode:
            self.update_ray_lines()
            
    def update_amplitude(self, value):
        
        self.amplitude = value
        


        
        y_range = self.getPlotItem().getViewBox().viewRange()[1]
        
        wave = self.calculate_wave(self.wavelength)
        self.wave_curve.setData(self.x, wave)
        
        self.getPlotItem().getViewBox().setYRange(y_range[0], y_range[1], padding=0)

    def update_speed(self, value):
        self.speed = value
        


    def update_n1(self, value):
        
        self.n1 = value
        
        for label_pair in self.medium_labels:
            if label_pair[2] == 1:  
                label_pair[0].setText(f"n₁: {self.n1:.4f}")
                break
    

        if not self.white_light:
            wave = self.calculate_wave(self.frequency)
            self.wave_curve.setData(self.x, wave)
        
        
        if self.show_ray_mode:
            self.update_ray_lines()
        
    def update_n2(self, value):
        
        self.n2 = value
        for label_pair in self.medium_labels:
            if label_pair[2] == 2:  
                label_pair[0].setText(f"n₂: {self.n2:.4f}")
                break
        
        if not self.white_light:
            wave = self.calculate_wave(self.frequency)
            self.wave_curve.setData(self.x, wave)
        
        if self.show_ray_mode:
            self.update_ray_lines()            

    def update_n3(self, value):
        
        self.n3 = value
        for label_pair in self.medium_labels:
            if label_pair[2] == 3:  
                label_pair[0].setText(f"n₃: {self.n3:.4f}")
                break
        
        
        if not self.white_light:
            wave = self.calculate_wave(self.frequency)
            self.wave_curve.setData(self.x, wave)
        
        if self.show_ray_mode:
            self.update_ray_lines()
        
        
    def toggle_prism_mode(self, enabled):
        
        self.prism_mode = enabled
        self.update_plot()
        
    def set_energy_label(self, label):
        
        self.external_energy_label = label

    def toggle_white_light(self, enabled):
        
        self.white_light = enabled
        if enabled:
            
            self.create_white_light_curves()

            
            if hasattr(self, 'wave_curve'):
                self.wave_curve.setVisible(False)
                
            
            self.superposition_enabled = True
            
            if self.superposition_wave is None:
                self.superposition_wave = self.plot(self.x, np.zeros_like(self.x), 
                                                   pen=pg.mkPen('w', width=4))
            
            self.superposition_wave.setVisible(True)
            for curve_item in self.wave_curves:
                if isinstance(curve_item, tuple) and len(curve_item) == 2:
                    curve, _ = curve_item
                    if hasattr(curve, 'setVisible'):
                        curve.setVisible(False)
                        
            self.update_superposition_wave()
            
            if hasattr(self, 'wave_widget') and hasattr(self.wave_widget, 'superposition_button'):
                self.wave_widget.superposition_button.setText("Show Components")
                self.wave_widget.superposition_button.setChecked(False)
        else:
            
            if hasattr(self, 'wave_curve'):
                self.wave_curve.setVisible(True)
                
            
            for curve_item in self.wave_curves:
                if isinstance(curve_item, tuple) and len(curve_item) == 2:
                    curve, _ = curve_item
                    if hasattr(curve, 'setVisible'):
                        curve.setVisible(False)
            
            if self.superposition_wave:
                self.superposition_wave.setVisible(False)
                
            if hasattr(self, 'wave_widget') and hasattr(self.wave_widget, 'superposition_button'):
                self.wave_widget.superposition_button.setText("Superposition")


    def create_white_light_curves(self):
        
        for curve_item in self.wave_curves:
            if isinstance(curve_item, tuple) and len(curve_item) == 2:
                curve, _ = curve_item
                if hasattr(curve, 'setVisible'):
                    self.removeItem(curve)  
        
        self.wave_curves = []
        
        self.setUpdatesEnabled(False)
        
        min_freq = 400  
        max_freq = 790 
        
        
        if hasattr(self, 'skip_frames') and self.skip_frames > 0:
            
            num_frequencies = 5
            sample_rate = 30  
        else:
            
            num_frequencies = 7
            sample_rate = 20  
        
        
        frequencies = np.linspace(min_freq, max_freq, num_frequencies)
        
        
        reduced_x = self.x[::sample_rate]
        
        for freq in frequencies:
            
            wavelength_nm = 299792.458 / freq
            
            
            wave = self.calculate_wave(freq)
            reduced_wave = wave[::sample_rate]
            
            
            r, g, b = wavelength_to_rgb(wavelength_nm)
            wave_color = QColor(r, g, b)
            
            curve = self.plot(reduced_x, reduced_wave, pen=pg.mkPen(wave_color, width=5))
            curve.setVisible(False)  
            
            self.wave_curves.append((curve, freq))
        
        
        self.setUpdatesEnabled(True)

    def update_superposition_wave(self):
        
        if not self.superposition_enabled or not self.white_light or self.superposition_wave is None:
            return
            
        
        cache_key = (self.time, self.n1, self.n2, self.n3)
        
        if hasattr(self, '_superposition_cache') and cache_key in self._superposition_cache:
            superposition = self._superposition_cache[cache_key]
        else:
            
            superposition = np.zeros_like(self.x)
            
            
            for curve, freq in self.wave_curves:
                wave = self.calculate_wave(freq)
                superposition += wave
                
            
            if len(self.wave_curves) > 0:
                superposition = superposition / len(self.wave_curves)
            
            
            if not hasattr(self, '_superposition_cache'):
                self._superposition_cache = {}
            self._superposition_cache[cache_key] = superposition
            
            
            if len(self._superposition_cache) > 50:
                oldest_keys = list(self._superposition_cache.keys())[:10]
                for key in oldest_keys:
                    self._superposition_cache.pop(key)
            
        
        self.superposition_wave.setData(self.x, superposition)

    def update_ray_lines(self):
        
        try:
            
            cache_key = (self.angle_of_incidence, self.n1, self.n2, self.n3, self.ray_target_y)
            
            if hasattr(self, '_ray_cache') and cache_key in self._ray_cache:
                ray_data = self._ray_cache[cache_key]
                
                
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
                
                
                self.angle1_label.setHtml(ray_data['angle1_html'])
                self.angle2_label.setHtml(ray_data['angle2_html'])
                self.angle3_label.setHtml(ray_data['angle3_html'])
                
                return
        

            angle_refraction1_deg = 0
            angle_refraction2_deg = 0


            label_html_style = """<div style="font-family: Arial; font-size: 16pt; font-weight: bold; color: white;">"""
            
            self.reflection_marker1.setVisible(False)
            self.reflection_marker2.setVisible(False)
            
            n1 = self.n1
            n2 = self.n2
            n3 = self.n3
            
            if not self.white_light:
                
                r, g, b = wavelength_to_rgb(self.wavelength)
                
                ray_color = pg.mkPen(QColor(r, g, b), width=4, style=Qt.DashLine)
                self.ray_incident.setPen(ray_color)
                self.ray_refracted1.setPen(ray_color)
                self.ray_refracted2.setPen(ray_color)
            else:
                
                self.ray_incident.setPen(pg.mkPen('#ffffff', width=4, style=Qt.DashLine))
                self.ray_refracted1.setPen(pg.mkPen('#00aaff', width=4, style=Qt.DashLine))
                self.ray_refracted2.setPen(pg.mkPen('#22ff22', width=4, style=Qt.DashLine))
            
            
            scale_factor = 0.008
            
            angle_incidence_rad = np.radians(self.angle_of_incidence)
            self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
            
            target_y = self.ray_target_y
            boundary1_y = target_y
            
            incident_ray_length = self.boundary1  
            start_x = 0  
            start_y = boundary1_y - np.tan(angle_incidence_rad) * incident_ray_length * scale_factor
            
            
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
            
            sin_theta2 = n1 * np.sin(angle_incidence_rad) / n2
            print(f"sin(theta2) = {sin_theta2:.6f}")
            
            self.ray_incident.setData([start_x, self.boundary1], [start_y, boundary1_y])
            self.ray_incident.setVisible(True)
            
            if abs(sin_theta2) >= 1.0:
                
                print(f"TIR DETECTED at boundary 1! |sin(theta2)| = {abs(sin_theta2):.6f} > 1.0")
                if critical_angle1_deg:
                    print(f"Incident angle: {self.angle_of_incidence}° > Critical angle: {critical_angle1_deg:.2f}°")
                
                self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
                self.angle2_label.setHtml(f"{label_html_style}θ₂: {angle_refraction1_deg:.1f}°</div>")
                self.angle3_label.setHtml(f"{label_html_style}θ₃: {angle_refraction2_deg:.1f}°</div>")
                
                
                
                
                reflected_end_x = start_x  
                reflected_end_y = -start_y  
                
                
                self.ray_refracted1.setData([self.boundary1, reflected_end_x], [boundary1_y, reflected_end_y])
                self.ray_refracted1.setVisible(True)
                
                
                self.reflection_marker1.setData([self.boundary1], [boundary1_y])
                self.reflection_marker1.setVisible(True)
                
                
                self.ray_refracted2.setVisible(False)
                
                return
            
            
            angle_refraction1_rad = np.arcsin(sin_theta2)
            angle_refraction1_deg = np.degrees(angle_refraction1_rad)
            print(f"Refraction angle at boundary 1: {angle_refraction1_deg:.2f}°")
            
            
            
            boundary_distance = self.boundary2 - self.boundary1
            
            
            
            refracted1_end_x = self.boundary2
            refracted1_end_y = boundary1_y + np.tan(angle_refraction1_rad) * boundary_distance * scale_factor
            
            
            self.ray_refracted1.setData([self.boundary1, refracted1_end_x], [boundary1_y, refracted1_end_y])
            self.ray_refracted1.setVisible(True)
            
            
            
            angle_incidence2_rad = angle_refraction1_rad
            angle_incidence2_deg = angle_refraction1_deg
            
            
            sin_theta3 = n2 * np.sin(angle_incidence2_rad) / n3
            print(f"sin(theta3) = {sin_theta3:.6f}")
            
            
            tir_at_boundary2 = False
            if critical_angle2 is not None:
                
                if abs(angle_incidence2_rad) > critical_angle2:
                    tir_at_boundary2 = True
                    print(f"TIR at boundary 2! |{angle_incidence2_deg:.2f}°| > {np.degrees(critical_angle2):.2f}°")
            
            
            if abs(sin_theta3) >= 1.0:
                tir_at_boundary2 = True
                print(f"TIR DETECTED at boundary 2! |sin(theta3)| = {abs(sin_theta3):.6f} > 1.0")
            
            if tir_at_boundary2:
                
                self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
                self.angle2_label.setHtml(f"{label_html_style}θ₂: {angle_refraction1_deg:.1f}°</div>")
                self.angle3_label.setHtml(f"{label_html_style}θ₃: TIR</div>")
                
                
                reflection_angle2_rad = -angle_incidence2_rad 
                
                
                reflected2_end_x = self.boundary1  
                reflected2_delta_x = self.boundary2 - reflected2_end_x
                reflected2_end_y = refracted1_end_y + np.tan(reflection_angle2_rad) * reflected2_delta_x * scale_factor
                
                
                self.ray_refracted2.setData([self.boundary2, reflected2_end_x], [refracted1_end_y, reflected2_end_y])
                self.ray_refracted2.setVisible(True)
                
                
                self.reflection_marker2.setData([self.boundary2], [refracted1_end_y])
                self.reflection_marker2.setVisible(True)
                
                return
                
            
            try:
                angle_refraction2_rad = np.arcsin(sin_theta3)
                angle_refraction2_deg = np.degrees(angle_refraction2_rad)
                print(f"Refraction angle at boundary 2: {angle_refraction2_deg:.2f}°")
                
                
                
                refracted2_length = 1000
                refracted2_end_x = self.boundary2 + refracted2_length
                refracted2_end_y = refracted1_end_y + np.tan(angle_refraction2_rad) * refracted2_length * scale_factor
                
                
                self.ray_refracted2.setData([self.boundary2, refracted2_end_x], [refracted1_end_y, refracted2_end_y])
                self.ray_refracted2.setVisible(True)
                
                
                self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
                self.angle2_label.setHtml(f"{label_html_style}θ₂: {angle_refraction1_deg:.1f}°</div>")
                self.angle3_label.setHtml(f"{label_html_style}θ₃: {angle_refraction2_deg:.1f}°</div>")

                if not hasattr(self, '_ray_cache'):
                    self._ray_cache = {}

                
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
                
                
                if len(self._ray_cache) > 100:
                    oldest_keys = list(self._ray_cache.keys())[:20]
                    for key in oldest_keys:
                        self._ray_cache.pop(key)

            except Exception as e:
                
                print(f"Exception in calculating refraction angle: {e}")
                
                
                self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
                self.angle2_label.setHtml(f"{label_html_style}θ₂: {angle_refraction1_deg:.1f}°</div>")
                self.angle3_label.setHtml(f"{label_html_style}θ₃: TIR (error)</div>")
                
                
                reflection_angle2_rad = -angle_incidence2_rad
                
                
                reflected2_end_x = self.boundary1
                reflected2_delta_x = self.boundary2 - reflected2_end_x
                reflected2_end_y = refracted1_end_y + np.tan(reflection_angle2_rad) * reflected2_delta_x * scale_factor
                
                
                self.ray_refracted2.setData([self.boundary2, reflected2_end_x], [refracted1_end_y, reflected2_end_y])
                self.ray_refracted2.setVisible(True)
                
                
                self.reflection_marker2.setData([self.boundary2], [refracted1_end_y])
                self.reflection_marker2.setVisible(True)
                
        except Exception as e:
            
            print(f"Error in ray calculation: {type(e).__name__}: {e}")
            
            self.ray_incident.setData([0, self.boundary1], [0, 0])
            self.ray_refracted1.setData([self.boundary1, self.boundary2], [0, 0])
            self.ray_refracted2.setData([self.boundary2, 3000], [0, 0])
            self.ray_incident.setVisible(True)
            self.ray_refracted1.setVisible(True)
            self.ray_refracted2.setVisible(True)
            
            self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
            self.angle2_label.setHtml(f"{label_html_style}θ₂: Error</div>")
            self.angle3_label.setHtml(f"{label_html_style}θ₃: Error</div>")

    def update_angle(self, value):

        self.angle_of_incidence = value
        
        
        cache_key = (self.angle_of_incidence, self.n1, self.n2, self.n3, self.ray_target_y)
        if hasattr(self, '_ray_cache') and cache_key in self._ray_cache:
            ray_data = self._ray_cache[cache_key]
            self.angle1_label.setHtml(ray_data['angle1_html'])
            self.angle2_label.setHtml(ray_data['angle2_html'])
            self.angle3_label.setHtml(ray_data['angle3_html'])
        else:

            
            label_html_style = """<div style="font-family: Arial; font-size: 16pt; font-weight: bold; color: white;">"""
            self.angle1_label.setHtml(f"{label_html_style}θ₁: {self.angle_of_incidence}°</div>")
            
            
            try:
                
                sin_theta2 = self.n1 * np.sin(np.radians(self.angle_of_incidence)) / self.n2
                
                if abs(sin_theta2) <= 1.0:  
                    angle_refraction1_deg = np.degrees(np.arcsin(sin_theta2))
                    self.angle2_label.setHtml(f"{label_html_style}θ₂: {angle_refraction1_deg:.1f}°</div>")
                    
                    
                    sin_theta3 = self.n2 * np.sin(np.radians(angle_refraction1_deg)) / self.n3
                    
                    if abs(sin_theta3) <= 1.0: 
                        angle_refraction2_deg = np.degrees(np.arcsin(sin_theta3))
                        self.angle3_label.setHtml(f"{label_html_style}θ₃: {angle_refraction2_deg:.1f}°</div>")
                    else:
                        self.angle3_label.setHtml(f"{label_html_style}θ₃: TIR</div>")
                else:
                    self.angle2_label.setHtml(f"{label_html_style}θ₂: TIR</div>")
                    self.angle3_label.setHtml(f"{label_html_style}θ₃: N/A</div>")
            except Exception as e:
                print(f"Error updating angle labels: {e}")
            
            
            if self.show_ray_mode:
                self.update_ray_lines()

    def set_ray_target_y(self, y_value):
        
        self.ray_target_y = y_value
        if self.show_ray_mode:
            self.update_ray_lines()

    def cleanup(self):
        
        if hasattr(self, 'timer') and self.timer.isActive():
            self.timer.stop()
        
        
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
        
        
        self.setWindowTitle("Light Wave Simulation")

        central_widget = QWidget()
        central_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setCentralWidget(central_widget)
        self.main_layout = QVBoxLayout(central_widget)
        self.main_layout.setContentsMargins(0, 0, 0, 0) 
        self.main_layout.setSpacing(0)

        
        self.medium_presets = {
            'Air': {'n': 1.0003, 'color': '#22222259'}, 
            'Water': {'n': 1.33, 'color': '#1E90FF59'},  
            'Glass (Crown)': {'n': 1.52, 'color': '#87CEEB59'}, 
            'Glass (Flint)': {'n': 1.62, 'color': '#9370DB59'}, 
            'Diamond': {'n': 2.42, 'color': '#E0E0E059'},  
            'Acrylic': {'n': 1.49, 'color': '#98FB9859'}, 
            'Glycerine': {'n': 1.47, 'color': '#F0E68C59'}, 
            'Ethanol': {'n': 1.36, 'color': '#D8BFD859'},  
            'Quartz': {'n': 1.54, 'color': '#DEB88759'}, 
            'Sapphire': {'n': 1.77, 'color': '#4682B459'}  
        }

        
        
        self.scenario_materials = {
            'Air → Water → Glass': ('Air', 'Water', 'Glass (Crown)'),
            'Air → Glass → Water': ('Air', 'Glass (Crown)', 'Water'),
            'Water → Air → Glass': ('Water', 'Air', 'Glass (Crown)'),
            'Air → Diamond → Glass': ('Air', 'Diamond', 'Glass (Crown)'),
            'Glass → Air → Water': ('Glass (Crown)', 'Air', 'Water')
        }
        
        
        
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
        
  
 
        
        
        self.wave_widget = WaveSimulationWidget()
        self.wave_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.main_layout.addWidget(self.wave_widget, 1)
        
        
        buttons_container = QWidget()
        buttons_layout = QHBoxLayout(buttons_container)
        buttons_layout.setContentsMargins(10, 5, 10, 5)

        
        self.energy_analysis_label = QLabel()
        self.energy_analysis_label.setStyleSheet("""
            background-color: rgba(0,0,0,0.8);
            padding: 8px;
            border-radius: 8px;
            border: 1px solid #3498db;
            color: white;
        """)
        self.energy_analysis_label.setMinimumWidth(300)
        self.energy_analysis_label.setMinimumHeight(60)
        self.energy_analysis_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)

        
        left_buttons = QWidget()
        left_buttons_layout = QVBoxLayout(left_buttons)
        left_buttons_layout.setContentsMargins(0, 0, 0, 0)
        left_buttons_layout.setSpacing(5)

        
        self.main_wave_check = ModernCheckBox("Show Main Wave")
        self.main_wave_check.setChecked(True)
        self.main_wave_check.toggled.connect(self.toggle_main_wave)
        
        self.reflections_check = ModernCheckBox("Show Reflections")
        self.reflections_check.setChecked(False)
        self.reflections_check.toggled.connect(self.toggle_reflections)
        
        
        left_buttons_layout.addWidget(self.main_wave_check)
        left_buttons_layout.addWidget(self.reflections_check)
        
        self.play_pause_button = PlayPauseButton()
        self.play_pause_button.clicked.connect(self.toggle_animation)


        buttons_layout.addWidget(left_buttons, 0)  
        buttons_layout.addWidget(self.energy_analysis_label, 0, Qt.AlignLeft)
        buttons_layout.addStretch(1)
        buttons_layout.addWidget(self.play_pause_button, 0, Qt.AlignCenter) 
        buttons_layout.addStretch(2)

  
        self.main_layout.addWidget(buttons_container)

 
        self.wave_widget.setContextMenuPolicy(Qt.NoContextMenu)

        
        self.wave_widget.medium_presets = self.medium_presets

        
        controls_container = QWidget()
        controls_container.setMaximumHeight(200)  
        controls_layout = QHBoxLayout()
        controls_layout.setContentsMargins(10, 10, 10, 10) 
        controls_layout.setSpacing(10) 
        controls_container.setLayout(controls_layout)
        self.main_layout.addWidget(controls_container)

        
        self.wave_widget.set_energy_label(self.energy_analysis_label)

        
        self.setup_wave_controls(controls_layout)
        
        self.wave_widget.medium1_color = self.medium_presets['Air']['color']
        self.wave_widget.medium2_color = self.medium_presets['Water']['color']
        self.wave_widget.medium3_color = self.medium_presets['Glass (Crown)']['color']
        self.wave_widget.toggle_main_wave(True)
        self.showMaximized()

    def resizeEvent(self, event):
        
        super().resizeEvent(event)
        
        if hasattr(self, 'pause_container') and hasattr(self, 'wave_widget'):
            self.pause_container.move(20, self.wave_widget.height() - 100)

    def toggle_animation(self):
        paused = self.play_pause_button.isChecked()
        self.wave_widget.toggle_pause(paused)
        self.play_pause_button.update_icon(paused)


    def toggle_main_wave(self, enabled):
        
        if hasattr(self.wave_widget, 'toggle_main_wave'):
            self.wave_widget.toggle_main_wave(enabled)
    
    def toggle_reflections(self, enabled):
        
        if hasattr(self.wave_widget, 'toggle_reflections'):
            self.wave_widget.toggle_reflections(enabled)
     
    def setup_wave_controls(self, layout):
        
        wave_group = QGroupBox("Wave Controls")
        wave_layout = QVBoxLayout()
        wave_group.setLayout(wave_layout)
        layout.addWidget(wave_group)
        
        
        
        
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
        wavelength_layout.addWidget(self.frequency_slider)  
        wavelength_layout.addWidget(self.frequency_value_label)
        wave_layout.addLayout(wavelength_layout)
        
        
        amplitude_layout = QHBoxLayout()
        amplitude_label = QLabel("Brightness/Amplitude:")
        self.amplitude_slider = BlueSlider()
        self.amplitude_slider.setMinimum(1)
        self.amplitude_slider.setMaximum(10)
        self.amplitude_slider.setValue(5)
        self.amplitude_value = QLabel("5")
        
        amplitude_layout.addWidget(amplitude_label)
        amplitude_layout.addWidget(self.amplitude_slider)
        amplitude_layout.addWidget(self.amplitude_value)
        wave_layout.addLayout(amplitude_layout)
        
        
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
        
        
        medium1_group = QGroupBox("Medium 1")
        medium1_layout = QVBoxLayout()
        medium1_group.setLayout(medium1_layout)
        layout.addWidget(medium1_group)
        
        
        self.medium1_combo = QComboBox()
        for medium in sorted(self.medium_presets.keys()):
            self.medium1_combo.addItem(medium)
        medium1_layout.addWidget(self.medium1_combo)
        

        
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

        
        medium2_group = QGroupBox("Medium 2")
        medium2_layout = QVBoxLayout()
        medium2_group.setLayout(medium2_layout)
        layout.addWidget(medium2_group)
        
        
        self.medium2_combo = QComboBox()
        for medium in sorted(self.medium_presets.keys()):
            self.medium2_combo.addItem(medium)
        medium2_layout.addWidget(self.medium2_combo)
        

        
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


        
        medium3_group = QGroupBox("Medium 3")
        medium3_layout = QVBoxLayout()
        medium3_group.setLayout(medium3_layout)
        layout.addWidget(medium3_group)
        

        
        self.medium3_combo = QComboBox()
        for medium in sorted(self.medium_presets.keys()):
            self.medium3_combo.addItem(medium)
        medium3_layout.addWidget(self.medium3_combo)
        
        
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
        
        enabled = state == Qt.Checked
        
        self.wave_widget.toggle_interference(enabled)
        
        
        if enabled:
            self.wave_widget.update_animation()
        
    def update_frequency(self, value):
        
        self.frequency_value_label.setText(f"{value} THz")
        
        wavelength_nm = 299792.458 / value
        self.wavelength_value.setText(f"{wavelength_nm:.1f} nm")
        
        self.wave_widget.frequency = value
        self.wave_widget.wavelength = wavelength_nm

        r, g, b = frequency_to_rgb(value)
        wave_color = QColor(r, g, b)
        

        self.wave_widget.update_wavelength(value)
        
    def update_amplitude(self, value):
        
        self.amplitude_value.setText(str(value))
        self.wave_widget.update_amplitude(value)
        
    def update_speed(self, value):
        
        self.speed_value.setText(str(value))
        self.wave_widget.update_speed(value)
        
    def update_n1(self, value):
        
        n = value / 100
        self.n1_value.setText(f"{n:.4f}")
        self.wave_widget.update_n1(n)
        self.wave_widget.update_plot()
        if self.wave_widget.show_ray_mode:
            self.wave_widget.update_ray_lines()
            
    def update_n2(self, value):
        
        n = value / 100
        self.n2_value.setText(f"{n:.4f}")
        self.wave_widget.update_n2(n)
        self.wave_widget.update_plot()
        if self.wave_widget.show_ray_mode:
            self.wave_widget.update_ray_lines()

    def update_n3(self, value):
        
        n = value / 100
        self.n3_value.setText(f"{n:.4f}")
        self.wave_widget.update_n3(n)
        self.wave_widget.update_plot()
        if self.wave_widget.show_ray_mode:
            self.wave_widget.update_ray_lines()

    def update_medium1(self, medium_name):
        
        if medium_name in self.medium_presets:
            
            n1 = self.medium_presets[medium_name]['n']
            self.n1_slider.setValue(int(n1 * 100))

            self.wave_widget.update_n1(n1)
            
            self.wave_widget.medium1_name = medium_name
            self.wave_widget.medium1_color = self.medium_presets[medium_name]['color']
            
            self.wave_widget.update_plot()
            
    def update_medium2(self, medium_name):
        
        if medium_name in self.medium_presets:
            
            n2 = self.medium_presets[medium_name]['n']
            self.n2_slider.setValue(int(n2 * 100))
            
            self.wave_widget.update_n2(n2)
            
            self.wave_widget.medium2_name = medium_name
            self.wave_widget.medium2_color = self.medium_presets[medium_name]['color']
            
            self.wave_widget.update_plot()
            
    def update_medium3(self, medium_name):
        
        if medium_name in self.medium_presets:
            
            n3 = self.medium_presets[medium_name]['n']
            self.n3_slider.setValue(int(n3 * 100))
            
            self.wave_widget.update_n3(n3)
            
            self.wave_widget.medium3_name = medium_name
            self.wave_widget.medium3_color = self.medium_presets[medium_name]['color']
            
            self.wave_widget.update_plot()
        
    def toggle_white_light(self, state):
        
        enabled = state == Qt.Checked
        self.wave_widget.toggle_white_light(enabled)

    def update_superposition_enabled(self, state):
        
        self.superposition_check.setEnabled(state == Qt.Checked)
        if state != Qt.Checked:
            self.superposition_check.setChecked(False)
            
    def toggle_superposition(self, state):
        
        enabled = state == Qt.Checked
        self.wave_widget.toggle_superposition(enabled)    
        
    def apply_scenario(self):
        
        scenario_name = self.scenario_combo.currentText()
        if scenario_name in self.scenario_materials and hasattr(self, 'medium_presets'):
            medium1, medium2, medium3 = self.scenario_materials[scenario_name]
            
            self.medium1_combo.setCurrentText(medium1)
            self.medium2_combo.setCurrentText(medium2)
            self.medium3_combo.setCurrentText(medium3)
                         
            self.wave_widget.update_plot()

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
        
        self.angle_value.setText(f"{value}°")
        
        self.wave_widget.update_angle(value)
        
        self.wave_widget.update_plot()
        
        if self.wave_widget.show_ray_mode:
            self.wave_widget.update_ray_lines()
            
    def update_ray_target(self, value):
        
        y_value = value / 100.0
        
        self.ray_target_value.setText(f"{y_value:.1f}")
        self.wave_widget.set_ray_target_y(y_value)
            
    def toggle_ray_mode(self, state):
        
        enabled = state == Qt.Checked
        self.wave_widget.toggle_ray_mode(enabled)

    def closeEvent(self, event):
        
        if hasattr(self.wave_widget, 'timer') and self.wave_widget.timer.isActive():
            self.wave_widget.timer.stop()
        
        if hasattr(self.wave_widget, 'clear'):
            self.wave_widget.clear()
        
        if hasattr(self.wave_widget, 'wave_curves'):
            for curve_item in self.wave_widget.wave_curves:
                if isinstance(curve_item, tuple) and len(curve_item) == 2:
                    curve, _ = curve_item
                    if hasattr(self.wave_widget, 'removeItem'):
                        self.wave_widget.removeItem(curve)
        
        event.accept()
        
        import gc
        gc.collect()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = LightSimulationApp()
    window.show()
    app.exec_()
    if hasattr(window, 'wave_widget') and hasattr(window.wave_widget, 'cleanup'):
        window.wave_widget.cleanup()