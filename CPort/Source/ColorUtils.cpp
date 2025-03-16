#include "../Header/ColorUtils.h"

QColor wavelengthToRGB(double wavelength) {
    // Convert wavelength (in nm) to RGB
    int r = 0, g = 0, b = 0;
    
    if (wavelength >= 380 && wavelength < 440) {
        r = static_cast<int>((440 - wavelength) / (440 - 380) * 255);
        g = 0;
        b = 255;
    } else if (wavelength >= 440 && wavelength < 490) {
        r = 0;
        g = static_cast<int>((wavelength - 440) / (490 - 440) * 255);
        b = 255;
    } else if (wavelength >= 490 && wavelength < 510) {
        r = 0;
        g = 255;
        b = static_cast<int>((510 - wavelength) / (510 - 490) * 255);
    } else if (wavelength >= 510 && wavelength < 580) {
        r = static_cast<int>((wavelength - 510) / (580 - 510) * 255);
        g = 255;
        b = 0;
    } else if (wavelength >= 580 && wavelength < 645) {
        r = 255;
        g = static_cast<int>((645 - wavelength) / (645 - 580) * 255);
        b = 0;
    } else if (wavelength >= 645 && wavelength <= 780) {
        r = 255;
        g = 0;
        b = 0;
    }
    
    return QColor(r, g, b);
}

QColor frequencyToRGB(double frequency) {
    // Convert frequency (in THz) to wavelength (in nm) and then to RGB
    double wavelength = 299792.458 / frequency;
    return wavelengthToRGB(wavelength);
}