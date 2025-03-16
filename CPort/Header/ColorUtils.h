#ifndef COLORUTILS_H
#define COLORUTILS_H

#include <QColor>

// Convert wavelength (in nm) to RGB
QColor wavelengthToRGB(double wavelength);

// Convert frequency (in THz) to RGB
QColor frequencyToRGB(double frequency);

#endif // COLORUTILS_H