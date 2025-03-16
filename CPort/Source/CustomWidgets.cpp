#include "../Header/CustomWidgets.h"

ColoredSlider::ColoredSlider(QWidget* parent) : QSlider(Qt::Horizontal, parent) {
    // Set up basic properties
}

void ColoredSlider::paintEvent(QPaintEvent* event) {
    // Draw the base slider
    QSlider::paintEvent(event);
    
    // Add color gradient
    QPainter painter(this);
    QStyleOptionSlider opt;
    initStyleOption(&opt);
    
    QRect groove = style()->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderGroove, this);
    
    // Create gradient from violet to red (visible spectrum)
    QLinearGradient gradient(groove.left(), 0, groove.right(), 0);
    gradient.setColorAt(0.0, QColor(148, 0, 211)); // Violet
    gradient.setColorAt(0.2, QColor(75, 0, 130));  // Indigo
    gradient.setColorAt(0.4, QColor(0, 0, 255));   // Blue
    gradient.setColorAt(0.6, QColor(0, 255, 0));   // Green
    gradient.setColorAt(0.8, QColor(255, 255, 0)); // Yellow
    gradient.setColorAt(1.0, QColor(255, 0, 0));   // Red
    
    painter.fillRect(groove, gradient);
    
    // Redraw the handle to make sure it's on top
    QRect handle = style()->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderHandle, this);
    painter.setPen(Qt::black);
    painter.setBrush(QColor(0, 122, 204)); // Blue handle
    painter.drawRoundedRect(handle, 3, 3);
}

BlueSlider::BlueSlider(QWidget* parent) : QSlider(Qt::Horizontal, parent) {
    setStyleSheet(
        "QSlider::groove:horizontal {"
        "    border: none;"
        "    height: 10px;"
        "    background: #333337;"
        "    margin: 2px 0;"
        "    border-radius: 5px;"
        "}"
        "QSlider::sub-page:horizontal {"
        "    background: #0088ff;"
        "    border: none;"
        "    height: 10px;"
        "    margin: 2px 0;"
        "    border-radius: 5px;"
        "}"
        "QSlider::handle:horizontal {"
        "    background: white;"
        "    border: none;"
        "    width: 18px;"
        "    margin: -4px 0;"
        "    border-radius: 9px;"
        "}"
    );
}

PlayPauseButton::PlayPauseButton(QWidget* parent) : QPushButton(parent), m_checked(false) {
    setFixedSize(40, 40);
    setCheckable(true);
    updateIcon(false);
}

bool PlayPauseButton::isChecked() const {
    return m_checked;
}

void PlayPauseButton::updateIcon(bool paused) {
    m_checked = paused;
    if (paused) {
        setStyleSheet(
            "QPushButton {"
            "    background-color: #007ACC;"
            "    border-radius: 20px;"
            "    border: none;"
            "    color: white;"
            "    font-size: 16px;"
            "    font-weight: bold;"
            "    text-align: center;"
            "}"
            "QPushButton:hover {"
            "    background-color: #1C97EA;"
            "}"
        );
        setText("▶");
    } else {
        setStyleSheet(
            "QPushButton {"
            "    background-color: #007ACC;"
            "    border-radius: 20px;"
            "    border: none;"
            "    color: white;"
            "    font-size: 16px;"
            "    font-weight: bold;"
            "    text-align: center;"
            "}"
            "QPushButton:hover {"
            "    background-color: #1C97EA;"
            "}"
        );
        setText("❚❚");
    }
}