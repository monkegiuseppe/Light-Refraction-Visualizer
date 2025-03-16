#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QSlider>
#include <QPushButton>
#include <QComboBox>
#include <QGroupBox>
#include <QCheckBox>
#include <QTimer>
#include <cmath>
#include <map>
#include <vector>
#include <tuple>
#include "../Header/qcustomplot.h"
#include "../Header/CustomWidgets.h"
#include "../Header/ColorUtils.h"


// Wave simulation widget
class WaveSimulationWidget : public QCustomPlot {
    Q_OBJECT
public:
    WaveSimulationWidget(QWidget* parent = nullptr);
    ~WaveSimulationWidget();

    // Medium properties
    std::map<QString, std::map<QString, QVariant>> medium_presets;
    QString medium1_name;
    QString medium2_name;
    QString medium3_name;
    QColor medium1_color;
    QColor medium2_color;
    QColor medium3_color;

    // Wave parameters
    double frequency;
    double wavelength;
    double amplitude;
    double speed;
    double angle_of_incidence;
    double n1;
    double n2;
    double n3;
    double ray_target_y;

    // Boundaries
    double boundary1;
    double boundary2;

    // Display flags
    bool show_ray_mode;
    bool white_light;
    bool show_interference;
    bool paused;

public slots:
    void update_plot();
    void update_animation();
    void update_wavelength(double freq);
    void update_amplitude(double amp);
    void update_speed(double spd);
    void update_n1(double n);
    void update_n2(double n);
    void update_n3(double n);
    void toggle_ray_mode(bool enabled);
    void toggle_white_light(bool enabled);
    void toggle_interference(bool enabled);
    void toggle_pause(bool paused);
    void update_angle(double angle);
    void set_ray_target_y(double y);
    void cleanup();

private:
    QTimer* timer;
    double phase;
    QVector<double> x;
    QCPGraph* wave_curve;
    std::vector<std::pair<QCPGraph*, double>> wave_curves;
    std::vector<double> prism_frequencies;

    // Medium rectangles
    QCPItemRect* medium1_rect;
    QCPItemRect* medium2_rect;
    QCPItemRect* medium3_rect;

    // Boundary lines
    QCPItemLine* boundary1_line;
    QCPItemLine* boundary2_line;

    // Ray lines
    QCPItemLine* ray_incident;
    QCPItemLine* ray_refracted1;
    QCPItemLine* ray_refracted2;
    QCPItemText* angle1_label;
    QCPItemText* angle2_label;
    QCPItemText* angle3_label;
    QCPGraph* reflection_marker1;
    QCPGraph* reflection_marker2;


    QVector<double> calculate_wave(double freq);
    void create_white_light_curves();
    void update_ray_lines();
};

WaveSimulationWidget::WaveSimulationWidget(QWidget* parent) : QCustomPlot(parent) {
    // Initialize parameters
    frequency = 545;
    wavelength = 299792.458 / frequency;
    amplitude = 2;
    speed = 0.5;
    angle_of_incidence = 0;
    n1 = 1.0003;
    n2 = 1.33;
    n3 = 1.52;
    show_ray_mode = false;
    white_light = false;
    show_interference = false;
    paused = false;
    ray_target_y = 0.0;
    phase = 0.0;
 

    // Set up plot appearance
    setBackground(QColor(45, 45, 48));
    xAxis->setBasePen(QPen(Qt::white));
    yAxis->setBasePen(QPen(Qt::white));
    xAxis->setTickPen(QPen(Qt::white));
    yAxis->setTickPen(QPen(Qt::white));
    xAxis->setSubTickPen(QPen(Qt::white));
    yAxis->setSubTickPen(QPen(Qt::white));
    xAxis->setTickLabelColor(Qt::white);
    yAxis->setTickLabelColor(Qt::white);
    
    // Set up plot ranges
    xAxis->setRange(0, 3000);
    yAxis->setRange(-5, 5);
    
    setAntialiasedElements(QCP::aeAll);
    setPlottingHint(QCP::phFastPolylines, true);
    setPlottingHint(QCP::phCacheLabels, true);
    setInteraction(QCP::iRangeDrag, false);
    setInteraction(QCP::iRangeZoom, false);
    setOpenGl(true); // Use OpenGL acceleration if available
    setNoAntialiasingOnDrag(false); // Keep antialiasing even during updates
    setBufferDevicePixelRatio(1.0);

    // Set up boundaries
    boundary1 = 1000;
    boundary2 = 2000;
    
    // Set up medium names and colors
    medium1_name = "Air";
    medium2_name = "Water";
    medium3_name = "Glass (Crown)";
    medium1_color = QColor(34, 34, 34, 89);
    medium2_color = QColor(30, 144, 255, 89);
    medium3_color = QColor(135, 206, 235, 89);
    
    // Create medium rectangles
    medium1_rect = new QCPItemRect(this);
    medium1_rect->setPen(Qt::NoPen);
    medium1_rect->setBrush(QBrush(medium1_color));
    medium1_rect->topLeft->setCoords(0, 10);
    medium1_rect->bottomRight->setCoords(boundary1, -10);
    
    medium2_rect = new QCPItemRect(this);
    medium2_rect->setPen(Qt::NoPen);
    medium2_rect->setBrush(QBrush(medium2_color));
    medium2_rect->topLeft->setCoords(boundary1, 10);
    medium2_rect->bottomRight->setCoords(boundary2, -10);
    
    medium3_rect = new QCPItemRect(this);
    medium3_rect->setPen(Qt::NoPen);
    medium3_rect->setBrush(QBrush(medium3_color));
    medium3_rect->topLeft->setCoords(boundary2, 10);
    medium3_rect->bottomRight->setCoords(3000, -10);
    
    // Create boundary lines
    boundary1_line = new QCPItemLine(this);
    boundary1_line->setPen(QPen(Qt::white, 2, Qt::DashLine));
    boundary1_line->start->setCoords(boundary1, 10);
    boundary1_line->end->setCoords(boundary1, -10);
    
    boundary2_line = new QCPItemLine(this);
    boundary2_line->setPen(QPen(Qt::white, 2, Qt::DashLine));
    boundary2_line->start->setCoords(boundary2, 10);
    boundary2_line->end->setCoords(boundary2, -10);
    
    // Create ray lines (initially hidden)
    ray_incident = new QCPItemLine(this);
    ray_incident->setPen(QPen(Qt::white, 2, Qt::DashLine));
    ray_incident->setVisible(false);
    
    ray_refracted1 = new QCPItemLine(this);
    ray_refracted1->setPen(QPen(Qt::white, 2, Qt::DashLine));
    ray_refracted1->setVisible(false);
    
    ray_refracted2 = new QCPItemLine(this);
    ray_refracted2->setPen(QPen(Qt::white, 2, Qt::DashLine));
    ray_refracted2->setVisible(false);
    
    // Create angle labels
    angle1_label = new QCPItemText(this);
    angle1_label->setPositionAlignment(Qt::AlignLeft|Qt::AlignBottom);
    angle1_label->position->setCoords(500, -6);
    angle1_label->setText("θ₁: 0°");
    angle1_label->setColor(Qt::white);
    angle1_label->setFont(QFont("Arial", 12, QFont::Bold));
    
    angle2_label = new QCPItemText(this);
    angle2_label->setPositionAlignment(Qt::AlignLeft|Qt::AlignBottom);
    angle2_label->position->setCoords(1500, -6);
    angle2_label->setText("θ₂: 0°");
    angle2_label->setColor(Qt::white);
    angle2_label->setFont(QFont("Arial", 12, QFont::Bold));
    
    angle3_label = new QCPItemText(this);
    angle3_label->setPositionAlignment(Qt::AlignLeft|Qt::AlignBottom);
    angle3_label->position->setCoords(2500, -6);
    angle3_label->setText("θ₃: 0°");
    angle3_label->setColor(Qt::white);
    angle3_label->setFont(QFont("Arial", 12, QFont::Bold));
    
    // Create reflection markers
    reflection_marker1 = addGraph();
    reflection_marker1->setPen(QPen(Qt::red, 10));
    reflection_marker1->setLineStyle(QCPGraph::lsNone);
    reflection_marker1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 10));
    reflection_marker1->setVisible(false);
    
    reflection_marker2 = addGraph();
    reflection_marker2->setPen(QPen(Qt::red, 10));
    reflection_marker2->setLineStyle(QCPGraph::lsNone);
    reflection_marker2->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 10));
    reflection_marker2->setVisible(false);
    // Create wave curve
    wave_curve = addGraph();
    wave_curve->setPen(QPen(frequencyToRGB(frequency), 3));
    wave_curve->setAntialiased(true);
    wave_curve->setAdaptiveSampling(true);
    wave_curve->setScatterStyle(QCPScatterStyle::ssNone);
    wave_curve->setLineStyle(QCPGraph::lsLine);

 
    // Initialize x values
    x.reserve(12000);
    for (int i = 0; i < 12000; i++) {
    x.append(i * 0.6); // Scale to keep the same visual range
}
    
    // Set up timer for animation
    timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, &WaveSimulationWidget::update_animation);
    timer->start(16);  // ~60 fps
    
    // Initialize white light frequencies
    prism_frequencies = {400, 450, 500, 550, 600, 650, 700, 750};
    
    // Initialize medium presets
    medium_presets["Air"] = {{"n", 1.0003}, {"color", QColor(34, 34, 34, 89)}};
    medium_presets["Water"] = {{"n", 1.33}, {"color", QColor(30, 144, 255, 89)}};
    medium_presets["Glass (Crown)"] = {{"n", 1.52}, {"color", QColor(135, 206, 235, 89)}};
    medium_presets["Glass (Flint)"] = {{"n", 1.62}, {"color", QColor(147, 112, 219, 89)}};
    medium_presets["Diamond"] = {{"n", 2.42}, {"color", QColor(224, 224, 224, 89)}};
    medium_presets["Acrylic"] = {{"n", 1.49}, {"color", QColor(152, 251, 152, 89)}};
    medium_presets["Glycerine"] = {{"n", 1.47}, {"color", QColor(240, 230, 140, 89)}};
    medium_presets["Ethanol"] = {{"n", 1.36}, {"color", QColor(216, 191, 216, 89)}};
    medium_presets["Quartz"] = {{"n", 1.54}, {"color", QColor(222, 184, 135, 89)}};
    medium_presets["Sapphire"] = {{"n", 1.77}, {"color", QColor(70, 130, 180, 89)}};
    
    // Initial wave calculation
    QVector<double> wave_data = calculate_wave(frequency);
    wave_curve->setData(x, wave_data);
}

WaveSimulationWidget::~WaveSimulationWidget() {
    cleanup();
}


void WaveSimulationWidget::cleanup() {
    if (timer && timer->isActive()) {
        timer->stop();
    }
    
    // Clear all plot items
    clearItems();
    clearGraphs();
}


QVector<double> WaveSimulationWidget::calculate_wave(double freq) {
    static QVector<double> wave;
    if (wave.size() != x.size()) {
        wave.resize(x.size());
    }
   
    
    // Calculate wave parameters
    double wl = 299792.458 / freq;
    double k1 = (2 * M_PI * n1) / wl;
    double k2 = (2 * M_PI * n2) / wl;
    double k3 = (2 * M_PI * n3) / wl;
    double phase_factor = speed * this->phase;
    double phase1 = k1 * boundary1 - phase_factor;
    double phase2 = k2 * (boundary2 - boundary1) + phase1;
    
    // Calculate angles once
    double angle_rad = angle_of_incidence * M_PI / 180.0;
    double sin_refraction1 = sin(angle_rad) * n1 / n2;
    bool total_reflection1 = false;
    bool total_reflection2 = false;
    double angle_refraction1 = 0;
    
    if (!total_reflection1) {
        angle_refraction1 = asin(sin_refraction1);
        const double sin_refraction2 = sin(angle_refraction1) * n2 / n3;
        total_reflection2 = fabs(sin_refraction2) > 1.0;
    }
    
    // Process regions in batches for better cache locality and vectorization
    const int size = x.size();
    const int b1_idx = std::lower_bound(x.begin(), x.end(), boundary1) - x.begin();
    const int b2_idx = std::lower_bound(x.begin(), x.end(), boundary2) - x.begin();
    
    // Medium 1 (vectorization-friendly loop)
    for (int i = 0; i < b1_idx; ++i) {
        wave[i] = amplitude * sin(k1 * x[i] - phase_factor);
    }
    
    // Medium 2
    if (total_reflection1) {
        for (int i = b1_idx; i < b2_idx; ++i) {
            const double x2 = x[i] - boundary1;
            wave[i] = amplitude * sin(k1 * x2 + phase1);
        }
    } else {
        for (int i = b1_idx; i < b2_idx; ++i) {
            const double x2 = x[i] - boundary1;
            wave[i] = amplitude * sin(k2 * x2 + phase1);
        }
    }
    
    // Medium 3
    if (total_reflection1) {
        for (int i = b2_idx; i < size; ++i) {
            wave[i] = 0;
        }
    } else if (total_reflection2) {
        for (int i = b2_idx; i < size; ++i) {
            const double x3 = x[i] - boundary2;
            wave[i] = amplitude * sin(k2 * x3 + phase2);
        }
    } else {
        for (int i = b2_idx; i < size; ++i) {
            const double x3 = x[i] - boundary2;
            wave[i] = amplitude * sin(k3 * x3 + phase2);
        }
    }
    
    return wave;
}

// Fix the update_animation method to properly declare last_phase
void WaveSimulationWidget::update_animation() {
    if (!paused) {
        phase += 0.05 * speed;
        
        // Only update the plot if it's visible
        if (isVisible()) {
            // Use a more efficient approach with frame limiting
            static QElapsedTimer frameTimer;
            static bool firstFrame = true;
            static double last_phase = 0.0;
            static double phase_threshold = 0.1; // Properly declare last_phase
            
            bool needsUpdate = std::abs(phase - last_phase) > phase_threshold;

            if (firstFrame || frameTimer.elapsed() > 16) { // ~60 FPS target
                firstFrame = false;
                frameTimer.restart();
                
                // Update wave data
                QVector<double> wave_data = calculate_wave(frequency);
                wave_curve->setData(x, wave_data, true);
                
                // Update white light curves if enabled
                if (white_light) {
                    for (auto& curve_pair : wave_curves) {
                        QCPGraph* curve = curve_pair.first;
                        double freq = curve_pair.second;
                        QVector<double> wave = calculate_wave(freq);
                        curve->setData(x, wave, true);
                    }
                }
                
                last_phase = phase;
                
                // Use a more efficient replot method
                replot(QCustomPlot::rpQueuedReplot);
            }
        }
    }
}

void WaveSimulationWidget::update_plot() {
    // Only update what's necessary
    bool needsReplot = false;
    
    // Check if medium rectangles need updating
    if (medium1_rect->bottomRight->coords().x() != boundary1 ||
        medium2_rect->topLeft->coords().x() != boundary1 ||
        medium2_rect->bottomRight->coords().x() != boundary2 ||
        medium3_rect->topLeft->coords().x() != boundary2) {
        
        // Update medium rectangles
        medium1_rect->setBrush(QBrush(medium1_color));
        medium1_rect->topLeft->setCoords(0, 10);
        medium1_rect->bottomRight->setCoords(boundary1, -10);
        
        medium2_rect->setBrush(QBrush(medium2_color));
        medium2_rect->topLeft->setCoords(boundary1, 10);
        medium2_rect->bottomRight->setCoords(boundary2, -10);
        
        medium3_rect->setBrush(QBrush(medium3_color));
        medium3_rect->topLeft->setCoords(boundary2, 10);
        medium3_rect->bottomRight->setCoords(3000, -10);
        
        // Update boundary lines
        boundary1_line->start->setCoords(boundary1, 10);
        boundary1_line->end->setCoords(boundary1, -10);
        
        boundary2_line->start->setCoords(boundary2, 10);
        boundary2_line->end->setCoords(boundary2, -10);
        
        needsReplot = true;
    }
    
    // Update wave data
    QVector<double> wave_data = calculate_wave(frequency);
    wave_curve->setData(x, wave_data, true); 
    needsReplot = true;
    
    // Update ray lines if enabled
    if (show_ray_mode) {
        update_ray_lines();
    }
    
    // Only replot if necessary
    if (needsReplot) {
        replot(QCustomPlot::rpQueuedReplot);
    }
}

void WaveSimulationWidget::update_wavelength(double freq) {
    frequency = freq;
    wavelength = 299792.458 / freq;
    
    // Update wave color
    QColor wave_color = frequencyToRGB(freq);
    wave_curve->setPen(QPen(wave_color, 4));
    
    // Update wave data
    QVector<double> wave_data = calculate_wave(freq);
    wave_curve->setData(x, wave_data);
    
    replot(QCustomPlot::rpQueuedReplot);
}

void WaveSimulationWidget::update_amplitude(double amp) {
    amplitude = amp;
    
    // Update wave data
    QVector<double> wave_data = calculate_wave(frequency);
    wave_curve->setData(x, wave_data);
    
    replot(QCustomPlot::rpQueuedReplot);
}

void WaveSimulationWidget::update_speed(double spd) {
    speed = spd;
}

void WaveSimulationWidget::update_n1(double n) {
    n1 = n;
    
    // Update wave data
    QVector<double> wave_data = calculate_wave(frequency);
    wave_curve->setData(x, wave_data);
    
    // Update ray lines if enabled
    if (show_ray_mode) {
        update_ray_lines();
    }
    
    replot(QCustomPlot::rpQueuedReplot);
}

void WaveSimulationWidget::update_n2(double n) {
    n2 = n;
    
    // Update wave data
    QVector<double> wave_data = calculate_wave(frequency);
    wave_curve->setData(x, wave_data);
    
    // Update ray lines if enabled
    if (show_ray_mode) {
        update_ray_lines();
    }
    
    replot(QCustomPlot::rpQueuedReplot);
}

void WaveSimulationWidget::update_n3(double n) {
    n3 = n;
    
    // Update wave data
    QVector<double> wave_data = calculate_wave(frequency);
    wave_curve->setData(x, wave_data);
    
    // Update ray lines if enabled
    if (show_ray_mode) {
        update_ray_lines();
    }
    
    replot(QCustomPlot::rpQueuedReplot);
}

void WaveSimulationWidget::toggle_ray_mode(bool enabled) {
    show_ray_mode = enabled;
    
    ray_incident->setVisible(enabled);
    ray_refracted1->setVisible(enabled);
    ray_refracted2->setVisible(enabled);
    angle1_label->setVisible(enabled);
    angle2_label->setVisible(enabled);
    angle3_label->setVisible(enabled);
    
    if (enabled) {
        update_ray_lines();
    }
    
    replot(QCustomPlot::rpQueuedReplot);
}

void WaveSimulationWidget::toggle_white_light(bool enabled) {
    white_light = enabled;
    
    if (enabled) {
        // Create white light curves
        create_white_light_curves();
        wave_curve->setVisible(false);
    } else {
        // Hide white light curves
        for (auto& curve_pair : wave_curves) {
            curve_pair.first->setVisible(false);
        }
        wave_curve->setVisible(true);
    }
    
    replot(QCustomPlot::rpQueuedReplot);
}

void WaveSimulationWidget::create_white_light_curves() {
    // Clear existing curves
    for (auto& curve_pair : wave_curves) {
        removeGraph(curve_pair.first);
    }
    wave_curves.clear();
    
    // Create new curves for each frequency
    for (double freq : prism_frequencies) {
        QCPGraph* curve = addGraph();
        QColor color = frequencyToRGB(freq);
        curve->setPen(QPen(color, 4));
        
        QVector<double> wave_data = calculate_wave(freq);
        curve->setData(x, wave_data);
        
        wave_curves.push_back(std::make_pair(curve, freq));
    }
}

void WaveSimulationWidget::toggle_interference(bool enabled) {
    show_interference = enabled;
    // Implementation would go here
}

void WaveSimulationWidget::toggle_pause(bool pause) {
    paused = pause;
}

void WaveSimulationWidget::update_angle(double angle) {
    angle_of_incidence = angle;
    
    // Update wave data
    QVector<double> wave_data = calculate_wave(frequency);
    wave_curve->setData(x, wave_data);
    
    // Update ray lines if enabled
    if (show_ray_mode) {
        update_ray_lines();
    }
    
    // Update angle labels
    double angle_rad = angle_of_incidence * M_PI / 180.0;
    double sin_refraction1 = sin(angle_rad) * n1 / n2;
    double angle_refraction1 = 0;
    double sin_refraction2 = 0;
    double angle_refraction2 = 0;
    
    if (fabs(sin_refraction1) <= 1.0) {
        angle_refraction1 = asin(sin_refraction1);
        sin_refraction2 = sin(angle_refraction1) * n2 / n3;
        
        if (fabs(sin_refraction2) <= 1.0) {
            angle_refraction2 = asin(sin_refraction2);
        }
    }
    
    angle1_label->setText(QString("θ₁: %1°").arg(angle_of_incidence, 0, 'f', 1));
    angle2_label->setText(QString("θ₂: %1°").arg(angle_refraction1 * 180.0 / M_PI, 0, 'f', 1));
    angle3_label->setText(QString("θ₃: %1°").arg(angle_refraction2 * 180.0 / M_PI, 0, 'f', 1));
    
    replot(QCustomPlot::rpQueuedReplot);
}

void WaveSimulationWidget::set_ray_target_y(double y) {
    ray_target_y = y;
    
    if (show_ray_mode) {
        update_ray_lines();
    }
}

void WaveSimulationWidget::update_ray_lines() {
    // Calculate angles
    double angle_rad = angle_of_incidence * M_PI / 180.0;
    double sin_refraction1 = sin(angle_rad) * n1 / n2;
    double angle_refraction1 = 0;
    double sin_refraction2 = 0;
    double angle_refraction2 = 0;
    bool total_reflection1 = false;
    bool total_reflection2 = false;
    
    if (fabs(sin_refraction1) > 1.0) {
        total_reflection1 = true;
        angle_refraction1 = M_PI / 2;
    } else {
        angle_refraction1 = asin(sin_refraction1);
        sin_refraction2 = sin(angle_refraction1) * n2 / n3;
        
        if (fabs(sin_refraction2) > 1.0) {
            total_reflection2 = true;
            angle_refraction2 = M_PI / 2;
        } else {
            angle_refraction2 = asin(sin_refraction2);
        }
    }
    
    // Calculate ray paths
    double start_x = 0;
    double start_y = ray_target_y;
    
    // Calculate incident ray
    double dx = boundary1 - start_x;
    double dy = dx * tan(angle_rad);
    double end_y = start_y + dy;
    
    ray_incident->start->setCoords(start_x, start_y);
    ray_incident->end->setCoords(boundary1, end_y);
    
    // Calculate first refracted ray
    if (total_reflection1) {
        // Total internal reflection at first boundary
        double reflection_angle = M_PI - angle_rad;
        double reflect_dx = boundary1;
        double reflect_dy = reflect_dx * tan(reflection_angle);
        ray_refracted1->start->setCoords(boundary1, end_y);
        ray_refracted1->end->setCoords(0, end_y + reflect_dy);
        
        // Show reflection marker
        QVector<double> rx(1), ry(1);
        rx[0] = boundary1;
        ry[0] = end_y;
        reflection_marker1->setData(rx, ry);
        reflection_marker1->setVisible(true);
        
        // Hide second refracted ray and reflection marker
        ray_refracted2->setVisible(false);
        reflection_marker2->setVisible(false);
    } else {
        // Regular refraction at first boundary
        double dx2 = boundary2 - boundary1;
        double dy2 = dx2 * tan(angle_refraction1);
        double end_y2 = end_y + dy2;
        
        ray_refracted1->start->setCoords(boundary1, end_y);
        ray_refracted1->end->setCoords(boundary2, end_y2);
        
        reflection_marker1->setVisible(false);
        
        // Calculate second refracted ray
        if (total_reflection2) {
            // Total internal reflection at second boundary
            double reflection_angle2 = M_PI - angle_refraction1;
            double reflect_dx2 = boundary2 - boundary1;
            double reflect_dy2 = reflect_dx2 * tan(reflection_angle2);
            ray_refracted2->start->setCoords(boundary2, end_y2);
            ray_refracted2->end->setCoords(boundary1, end_y2 + reflect_dy2);
            
            // Show reflection marker
            QVector<double> rx(1), ry(1);
            rx[0] = boundary2;
            ry[0] = end_y2;
            reflection_marker2->setData(rx, ry);
            reflection_marker2->setVisible(true);
        } else {
            // Regular refraction at second boundary
            double dx3 = 3000 - boundary2;
            double dy3 = dx3 * tan(angle_refraction2);
            double end_y3 = end_y2 + dy3;
            
            ray_refracted2->start->setCoords(boundary2, end_y2);
            ray_refracted2->end->setCoords(3000, end_y3);
            
            reflection_marker2->setVisible(false);
        }
        
        ray_refracted2->setVisible(true);
    }
}

// Main application class
class LightSimulationApp : public QMainWindow {
    Q_OBJECT
public:
    LightSimulationApp(QWidget* parent = nullptr);
    ~LightSimulationApp();

private slots:
    void toggle_animation();
    void update_frequency(int value);
    void update_amplitude(int value);
    void update_speed(int value);
    void update_angle(int value);
    void update_ray_target(int value);
    void update_n1(int value);
    void update_n2(int value);
    void update_n3(int value);
    void update_medium1(const QString& medium_name);
    void update_medium2(const QString& medium_name);
    void update_medium3(const QString& medium_name);
    void toggle_white_light(int state);
    void toggle_interference(int state);
    void toggle_ray_mode(int state);
    void apply_scenario();

private:
    void setup_wave_controls(QHBoxLayout* layout);
    
    QVBoxLayout* main_layout;
    WaveSimulationWidget* wave_widget;
    PlayPauseButton* play_pause_button;
    
    // Wave controls
    ColoredSlider* frequency_slider;
    QLabel* frequency_value_label;
    QLabel* wavelength_value;
    BlueSlider* amplitude_slider;
    QLabel* amplitude_value;
    BlueSlider* speed_slider;
    QLabel* speed_value;
    BlueSlider* angle_slider;
    QLabel* angle_value;
    BlueSlider* ray_target_slider;
    QLabel* ray_target_value;
    
    // Medium controls
    QComboBox* medium1_combo;
    BlueSlider* n1_slider;
    QLabel* n1_value;
    QComboBox* medium2_combo;
    BlueSlider* n2_slider;
    QLabel* n2_value;
    QComboBox* medium3_combo;
    BlueSlider* n3_slider;
    QLabel* n3_value;
    
    // Mode controls
    QCheckBox* white_light_check;
    QCheckBox* interference_check;
    QCheckBox* ray_mode_check;
    
    // Scenario controls
    QComboBox* scenario_combo;
    
    // Medium presets
    std::map<QString, std::map<QString, QVariant>> medium_presets;
    
    // Scenario presets
    std::map<QString, std::tuple<QString, QString, QString>> scenario_materials;
};

LightSimulationApp::LightSimulationApp(QWidget* parent) : QMainWindow(parent) {
    setWindowTitle("Light Simulation");
    
    // Initialize medium presets
    medium_presets["Air"] = {{"n", 1.0003}, {"color", QColor(34, 34, 34, 89)}};
    medium_presets["Water"] = {{"n", 1.33}, {"color", QColor(30, 144, 255, 89)}};
    medium_presets["Glass (Crown)"] = {{"n", 1.52}, {"color", QColor(135, 206, 235, 89)}};
    medium_presets["Glass (Flint)"] = {{"n", 1.62}, {"color", QColor(147, 112, 219, 89)}};
    medium_presets["Diamond"] = {{"n", 2.42}, {"color", QColor(224, 224, 224, 89)}};
    medium_presets["Acrylic"] = {{"n", 1.49}, {"color", QColor(152, 251, 152, 89)}};
    medium_presets["Glycerine"] = {{"n", 1.47}, {"color", QColor(240, 230, 140, 89)}};
    medium_presets["Ethanol"] = {{"n", 1.36}, {"color", QColor(216, 191, 216, 89)}};
    medium_presets["Quartz"] = {{"n", 1.54}, {"color", QColor(222, 184, 135, 89)}};
    medium_presets["Sapphire"] = {{"n", 1.77}, {"color", QColor(70, 130, 180, 89)}};
    
    // Initialize scenario presets
    scenario_materials["Air → Water → Glass"] = std::make_tuple("Air", "Water", "Glass (Crown)");
    scenario_materials["Air → Glass → Water"] = std::make_tuple("Air", "Glass (Crown)", "Water");
    scenario_materials["Water → Air → Glass"] = std::make_tuple("Water", "Air", "Glass (Crown)");
    scenario_materials["Air → Diamond → Glass"] = std::make_tuple("Air", "Diamond", "Glass (Crown)");
    scenario_materials["Glass → Air → Water"] = std::make_tuple("Glass (Crown)", "Air", "Water");
    
    // Set dark mode stylesheet
    setStyleSheet(
        "QMainWindow, QWidget {"
        "    background-color: #2D2D30;"
        "    color: #FFFFFF;"
        "}"
        "QGroupBox {"
        "    border: 1px solid #3F3F46;"
        "    border-radius: 5px;"
        "    margin-top: 10px;"
        "    font-weight: bold;"
        "}"
        "QGroupBox::title {"
        "    subcontrol-origin: margin;"
        "    left: 10px;"
        "    padding: 0 5px 0 5px;"
        "}"
        "QPushButton {"
        "    background-color: #007ACC;"
        "    border: none;"
        "    border-radius: 3px;"
        "    padding: 5px;"
        "    color: white;"
        "}"
        "QPushButton:hover {"
        "    background-color: #1C97EA;"
        "}"
        "QComboBox, QLineEdit, QSpinBox, QDoubleSpinBox {"
        "    background-color: #333337;"
        "    border: 1px solid #3F3F46;"
        "    border-radius: 3px;"
        "    padding: 2px;"
        "}"
        "QSlider::groove:horizontal {"
        "    border: 1px solid #999999;"
        "    height: 8px;"
        "    background: #333337;"
        "    margin: 2px 0;"
        "    border-radius: 2px;"
        "}"
        "QSlider::handle:horizontal {"
        "    background: #007ACC;"
        "    border: 1px solid #5c5c5c;"
        "    width: 18px;"
        "    margin: -2px 0;"
        "    border-radius: 3px;"
        "}"
    );
    
    // Create central widget and main layout
    QWidget* central_widget = new QWidget(this);
    setCentralWidget(central_widget);
    main_layout = new QVBoxLayout(central_widget);
    
    // Create wave widget
    wave_widget = new WaveSimulationWidget();
    wave_widget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    main_layout->addWidget(wave_widget, 1);
    
    // Create play/pause button
    play_pause_button = new PlayPauseButton();
    connect(play_pause_button, &QPushButton::clicked, this, &LightSimulationApp::toggle_animation);
    
    // Add button to centered container
    QWidget* button_container = new QWidget();
    QHBoxLayout* button_layout = new QHBoxLayout(button_container);
    button_layout->addStretch(1);
    button_layout->addWidget(play_pause_button);
    button_layout->addStretch(1);
    button_layout->setContentsMargins(0, 5, 0, 5);
    
    main_layout->addWidget(button_container);
    
    // Share medium presets with wave widget
    wave_widget->medium_presets = medium_presets;
    
    // Create controls container
    QWidget* controls_container = new QWidget();
    controls_container->setMaximumHeight(200);
    QHBoxLayout* controls_layout = new QHBoxLayout();
    controls_layout->setContentsMargins(10, 10, 10, 10);
    controls_layout->setSpacing(10);
    controls_container->setLayout(controls_layout);
    main_layout->addWidget(controls_container);
    
    // Add wave controls
    setup_wave_controls(controls_layout);
    
    // Set initial medium colors
    wave_widget->medium1_color = medium_presets["Air"]["color"].value<QColor>();
    wave_widget->medium2_color = medium_presets["Water"]["color"].value<QColor>();
    wave_widget->medium3_color = medium_presets["Glass (Crown)"]["color"].value<QColor>();
    

    showMaximized();
}

LightSimulationApp::~LightSimulationApp() {
    // Cleanup is handled by Qt's parent-child relationship
}

void LightSimulationApp::toggle_animation() {
    bool paused = play_pause_button->isChecked();
    wave_widget->toggle_pause(paused);
    play_pause_button->updateIcon(paused);
}

void LightSimulationApp::setup_wave_controls(QHBoxLayout* layout) {
    // Left controls group (wave properties)
    QGroupBox* wave_group = new QGroupBox("Wave Controls");
    QVBoxLayout* wave_layout = new QVBoxLayout();
    wave_group->setLayout(wave_layout);
    
    // Frequency slider with color gradient
    QHBoxLayout* wavelength_layout = new QHBoxLayout();
    QLabel* frequency_label = new QLabel("Frequency (THz):");
    frequency_slider = new ColoredSlider();
    frequency_slider->setMinimum(400);
    frequency_slider->setMaximum(790);
    frequency_slider->setValue(545);
    frequency_slider->setTickPosition(QSlider::TicksBelow);
    frequency_slider->setTickInterval(50);
    
    wavelength_value = new QLabel("550 nm");
    frequency_value_label = new QLabel("545 THz");
    
    wavelength_layout->addWidget(frequency_label);
    wavelength_layout->addWidget(frequency_slider);
    wavelength_layout->addWidget(frequency_value_label);
    wave_layout->addLayout(wavelength_layout);
    
    // Amplitude slider
    QHBoxLayout* amplitude_layout = new QHBoxLayout();
    QLabel* amplitude_label = new QLabel("Amplitude:");
    amplitude_slider = new BlueSlider();
    amplitude_slider->setMinimum(1);
    amplitude_slider->setMaximum(5);
    amplitude_slider->setValue(2);
    amplitude_value = new QLabel("2");
    
    amplitude_layout->addWidget(amplitude_label);
    amplitude_layout->addWidget(amplitude_slider);
    amplitude_layout->addWidget(amplitude_value);
    wave_layout->addLayout(amplitude_layout);
    
    // Speed slider
    QHBoxLayout* speed_layout = new QHBoxLayout();
    QLabel* speed_label = new QLabel("Wave Speed:");
    speed_slider = new BlueSlider();
    speed_slider->setMinimum(1);
    speed_slider->setMaximum(10);
    speed_slider->setValue(5);
    speed_value = new QLabel("0.5");
    
    speed_layout->addWidget(speed_label);
    speed_layout->addWidget(speed_slider);
    speed_layout->addWidget(speed_value);
    wave_layout->addLayout(speed_layout);
    
    // Angle of incidence slider
    QHBoxLayout* angle_layout = new QHBoxLayout();
    QLabel* angle_label = new QLabel("Angle of Incidence:");
    angle_slider = new BlueSlider();
    angle_slider->setMinimum(-90);
    angle_slider->setMaximum(90);
    angle_slider->setValue(0);
    angle_value = new QLabel("0°");
    
    angle_layout->addWidget(angle_label);
    angle_layout->addWidget(angle_slider);
    angle_layout->addWidget(angle_value);
    wave_layout->addLayout(angle_layout);
    
    // Ray target Y position slider
    QHBoxLayout* ray_target_layout = new QHBoxLayout();
    QLabel* ray_target_label = new QLabel("Ray Target Y:");
    ray_target_slider = new BlueSlider();
    ray_target_slider->setMinimum(-100);
    ray_target_slider->setMaximum(100);    // Calculate incident ray
        ray_target_slider->setMinimum(-100);
        ray_target_slider->setMaximum(100);
        ray_target_slider->setValue(0);
        ray_target_value = new QLabel("0.0");
        
        ray_target_layout->addWidget(ray_target_label);
        ray_target_layout->addWidget(ray_target_slider);
        ray_target_layout->addWidget(ray_target_value);
        wave_layout->addLayout(ray_target_layout);
        
        // Mode checkboxes
        QHBoxLayout* mode_layout = new QHBoxLayout();
        white_light_check = new QCheckBox("White Light");
        interference_check = new QCheckBox("Interference");
        ray_mode_check = new QCheckBox("Ray Mode");
        
        mode_layout->addWidget(white_light_check);
        mode_layout->addWidget(interference_check);
        mode_layout->addWidget(ray_mode_check);
        wave_layout->addLayout(mode_layout);
        
        // Add wave group to main layout
        layout->addWidget(wave_group);
        
        // Medium 1 group
        QGroupBox* medium1_group = new QGroupBox("Medium 1");
        QVBoxLayout* medium1_layout = new QVBoxLayout();
        medium1_group->setLayout(medium1_layout);
        
        // Medium 1 selection
        medium1_combo = new QComboBox();
        for (const auto& medium : medium_presets) {
            medium1_combo->addItem(medium.first);
        }
        medium1_layout->addWidget(medium1_combo);
        
        // Medium 1 n slider
        QHBoxLayout* n1_layout = new QHBoxLayout();
        QLabel* n1_label = new QLabel("n₁:");
        n1_slider = new BlueSlider();
        n1_slider->setMinimum(100);
        n1_slider->setMaximum(300);
        n1_slider->setValue(100);
        n1_value = new QLabel("1.0003");
        
        n1_layout->addWidget(n1_label);
        n1_layout->addWidget(n1_slider);
        n1_layout->addWidget(n1_value);
        medium1_layout->addLayout(n1_layout);
        
        // Set initial medium
        medium1_combo->setCurrentText("Air");
        
        // Add medium 1 group to main layout
        layout->addWidget(medium1_group);
        
        // Medium 2 group
        QGroupBox* medium2_group = new QGroupBox("Medium 2");
        QVBoxLayout* medium2_layout = new QVBoxLayout();
        medium2_group->setLayout(medium2_layout);
        
        // Medium 2 selection
        medium2_combo = new QComboBox();
        for (const auto& medium : medium_presets) {
            medium2_combo->addItem(medium.first);
        }
        medium2_layout->addWidget(medium2_combo);
        
        // Medium 2 n slider
        QHBoxLayout* n2_layout = new QHBoxLayout();
        QLabel* n2_label = new QLabel("n₂:");
        n2_slider = new BlueSlider();
        n2_slider->setMinimum(100);
        n2_slider->setMaximum(300);
        n2_slider->setValue(133);
        n2_value = new QLabel("1.33");
        
        n2_layout->addWidget(n2_label);
        n2_layout->addWidget(n2_slider);
        n2_layout->addWidget(n2_value);
        medium2_layout->addLayout(n2_layout);
        
        // Set initial medium
        medium2_combo->setCurrentText("Water");
        
        // Add medium 2 group to main layout
        layout->addWidget(medium2_group);
        
        // Medium 3 group
        QGroupBox* medium3_group = new QGroupBox("Medium 3");
        QVBoxLayout* medium3_layout = new QVBoxLayout();
        medium3_group->setLayout(medium3_layout);
        
        // Medium 3 selection
        medium3_combo = new QComboBox();
        for (const auto& medium : medium_presets) {
            medium3_combo->addItem(medium.first);
        }
        medium3_layout->addWidget(medium3_combo);
        
        // Medium 3 n slider
        QHBoxLayout* n3_layout = new QHBoxLayout();
        QLabel* n3_label = new QLabel("n₃:");
        n3_slider = new BlueSlider();
        n3_slider->setMinimum(100);
        n3_slider->setMaximum(300);
        n3_slider->setValue(152);
        n3_value = new QLabel("1.52");
        
        n3_layout->addWidget(n3_label);
        n3_layout->addWidget(n3_slider);
        n3_layout->addWidget(n3_value);
        medium3_layout->addLayout(n3_layout);
        
        // Set initial medium
        medium3_combo->setCurrentText("Glass (Crown)");
        
        // Add medium 3 group to main layout
        layout->addWidget(medium3_group);
        
        // Scenario presets group
        QGroupBox* scenario_group = new QGroupBox("Presets");
        QVBoxLayout* scenario_layout = new QVBoxLayout();
        scenario_group->setLayout(scenario_layout);
        
        // Scenario selection
        scenario_combo = new QComboBox();
        for (const auto& scenario : scenario_materials) {
            scenario_combo->addItem(scenario.first);
        }
        scenario_layout->addWidget(scenario_combo);
        
        // Apply button
        QPushButton* apply_button = new QPushButton("Apply Preset");
        scenario_layout->addWidget(apply_button);
        
        // Add scenario group to main layout
        layout->addWidget(scenario_group);
        
        // Connect signals
        connect(frequency_slider, &QSlider::valueChanged, this, &LightSimulationApp::update_frequency);
        connect(amplitude_slider, &QSlider::valueChanged, this, &LightSimulationApp::update_amplitude);
        connect(speed_slider, &QSlider::valueChanged, this, &LightSimulationApp::update_speed);
        connect(angle_slider, &QSlider::valueChanged, this, &LightSimulationApp::update_angle);
        connect(ray_target_slider, &QSlider::valueChanged, this, &LightSimulationApp::update_ray_target);
        
        connect(n1_slider, &QSlider::valueChanged, this, &LightSimulationApp::update_n1);
        connect(n2_slider, &QSlider::valueChanged, this, &LightSimulationApp::update_n2);
        connect(n3_slider, &QSlider::valueChanged, this, &LightSimulationApp::update_n3);
        
        connect(medium1_combo, &QComboBox::currentTextChanged, this, &LightSimulationApp::update_medium1);
        connect(medium2_combo, &QComboBox::currentTextChanged, this, &LightSimulationApp::update_medium2);
        connect(medium3_combo, &QComboBox::currentTextChanged, this, &LightSimulationApp::update_medium3);
        
        connect(white_light_check, &QCheckBox::stateChanged, this, &LightSimulationApp::toggle_white_light);
        connect(interference_check, &QCheckBox::stateChanged, this, &LightSimulationApp::toggle_interference);
        connect(ray_mode_check, &QCheckBox::stateChanged, this, &LightSimulationApp::toggle_ray_mode);
        
        connect(apply_button, &QPushButton::clicked, this, &LightSimulationApp::apply_scenario);
    }
    
    void LightSimulationApp::update_frequency(int value) {
        // Update displayed values
        frequency_value_label->setText(QString("%1 THz").arg(value));
        double wavelength_nm = 299792.458 / value;
        wavelength_value->setText(QString("%1 nm").arg(wavelength_nm, 0, 'f', 1));
        
        // Update wave widget
        wave_widget->update_wavelength(value);
    }
    
    void LightSimulationApp::update_amplitude(int value) {
        amplitude_value->setText(QString::number(value));
        wave_widget->update_amplitude(value);
    }
    
    void LightSimulationApp::update_speed(int value) {
        // Convert slider value (1-10) to actual speed (0.1-1.0)
        double actual_speed = value * 0.1;
        speed_value->setText(QString::number(actual_speed, 'f', 1));
        wave_widget->update_speed(actual_speed);
    }
    
    void LightSimulationApp::update_angle(int value) {
        angle_value->setText(QString("%1°").arg(value));
        wave_widget->update_angle(value);
    }
    
    void LightSimulationApp::update_ray_target(int value) {
        // Scale the slider value to a smaller range for better control
        double y_value = value / 100.0;
        ray_target_value->setText(QString("%1").arg(y_value, 0, 'f', 1));
        wave_widget->set_ray_target_y(y_value);
    }
    
    void LightSimulationApp::update_n1(int value) {
        double n = value / 100.0;
        n1_value->setText(QString("%1").arg(n, 0, 'f', 4));
        wave_widget->update_n1(n);
    }
    
    void LightSimulationApp::update_n2(int value) {
        double n = value / 100.0;
        n2_value->setText(QString("%1").arg(n, 0, 'f', 4));
        wave_widget->update_n2(n);
    }
    
    void LightSimulationApp::update_n3(int value) {
        double n = value / 100.0;
        n3_value->setText(QString("%1").arg(n, 0, 'f', 4));
        wave_widget->update_n3(n);
    }
    
    void LightSimulationApp::update_medium1(const QString& medium_name) {
        if (medium_presets.find(medium_name) != medium_presets.end()) {
            // Update the refractive index
            double n1 = medium_presets[medium_name]["n"].toDouble();
            n1_slider->setValue(static_cast<int>(n1 * 100));
            
            // Update the wave widget's medium name and color
            wave_widget->medium1_name = medium_name;
            wave_widget->medium1_color = medium_presets[medium_name]["color"].value<QColor>();
            
            // Update the wave widget's plot
            wave_widget->update_plot();
        }
    }
    
    void LightSimulationApp::update_medium2(const QString& medium_name) {
        if (medium_presets.find(medium_name) != medium_presets.end()) {
            // Update the refractive index
            double n2 = medium_presets[medium_name]["n"].toDouble();
            n2_slider->setValue(static_cast<int>(n2 * 100));
            
            // Update the wave widget's medium name and color
            wave_widget->medium2_name = medium_name;
            wave_widget->medium2_color = medium_presets[medium_name]["color"].value<QColor>();
            
            // Update the wave widget's plot
            wave_widget->update_plot();
        }
    }
    
    void LightSimulationApp::update_medium3(const QString& medium_name) {
        if (medium_presets.find(medium_name) != medium_presets.end()) {
            // Update the refractive index
            double n3 = medium_presets[medium_name]["n"].toDouble();
            n3_slider->setValue(static_cast<int>(n3 * 100));
            
            // Update the wave widget's medium name and color
            wave_widget->medium3_name = medium_name;
            wave_widget->medium3_color = medium_presets[medium_name]["color"].value<QColor>();
            
            // Update the wave widget's plot
            wave_widget->update_plot();
        }
    }
    
    void LightSimulationApp::toggle_white_light(int state) {
        bool enabled = (state == Qt::Checked);
        wave_widget->toggle_white_light(enabled);
    }
    
    void LightSimulationApp::toggle_interference(int state) {
        bool enabled = (state == Qt::Checked);
        wave_widget->toggle_interference(enabled);
    }
    
    void LightSimulationApp::toggle_ray_mode(int state) {
        bool enabled = (state == Qt::Checked);
        wave_widget->toggle_ray_mode(enabled);
    }
    
    void LightSimulationApp::apply_scenario() {
        QString scenario_name = scenario_combo->currentText();
        if (scenario_materials.find(scenario_name) != scenario_materials.end()) {
            auto [medium1, medium2, medium3] = scenario_materials[scenario_name];
            
            // Update medium selections
            medium1_combo->setCurrentText(medium1);
            medium2_combo->setCurrentText(medium2);
            medium3_combo->setCurrentText(medium3);
        }
    }
    
    // Main function
    int main(int argc, char *argv[]) {
        QApplication app(argc, argv);
        LightSimulationApp window;
        window.show();
        return app.exec();
    }
    
    #include "LightVisC.moc"
