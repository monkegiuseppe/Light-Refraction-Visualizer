QT       += core gui printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17
CONFIG += optimize_full


# Enable OpenMP for parallel processing
QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

# Add these optimization flags
QMAKE_CXXFLAGS += -O3 -march=native -mtune=native -ffast-math

SOURCES += \
    Source/LightVisC.cpp \
    Source/CustomWidgets.cpp \
    Source/ColorUtils.cpp \
    Source/qcustomplot.cpp

HEADERS += \
    Header/CustomWidgets.h \
    Header/ColorUtils.h \
    Header/qcustomplot.h

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target