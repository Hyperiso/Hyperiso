#ifndef __ISPECTRUMCALCULATOR_H__
#define __ISPECTRUMCALCULATOR_H__

#include <filesystem>
#include "General.h"

namespace fs = std::filesystem;

class ISpectrumCalculator {
public:
    virtual void calculate_spectrum(fs::path in_lha_path, fs::path out_spectrum_path, Model model) = 0;

    virtual ~ISpectrumCalculator() = default;
};

#endif // __ISPECTRUMCALCULATOR_H__
