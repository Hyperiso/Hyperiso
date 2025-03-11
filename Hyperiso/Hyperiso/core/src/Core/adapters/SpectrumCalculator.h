#ifndef __SPECTRUMCALCULATOR_H__
#define __SPECTRUMCALCULATOR_H__

#include "ISpectrumCalculator.h"
#include "Interface.h"

class SpectrumCalculator : public ISpectrumCalculator {
public:
    void calculate_spectrum(fs::path in_lha_path, fs::path out_spectrum_path, Model model) override;
};


#endif // __SPECTRUMCALCULATOR_H__
