#ifndef CUSTOMWIDGETS_H
#define CUSTOMWIDGETS_H

#include <QSlider>
#include <QPushButton>
#include <QPainter>
#include <QStyle>
#include <QStyleOptionSlider>

// A slider with blue styling
class BlueSlider : public QSlider {
    Q_OBJECT
public:
    BlueSlider(QWidget* parent = nullptr);
};

// A slider with color gradient for wavelength visualization
class ColoredSlider : public QSlider {
    Q_OBJECT
public:
    ColoredSlider(QWidget* parent = nullptr);
    
protected:
    void paintEvent(QPaintEvent* event) override;
};

// A play/pause button that toggles state
class PlayPauseButton : public QPushButton {
    Q_OBJECT
public:
    PlayPauseButton(QWidget* parent = nullptr);
    bool isChecked() const;
    
public slots:
    void updateIcon(bool paused);
    
private:
    bool m_checked;
};

#endif // CUSTOMWIDGETS_H