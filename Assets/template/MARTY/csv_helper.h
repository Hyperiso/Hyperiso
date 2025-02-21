#include <string>
#include <complex>
#include <fstream>
#include "common.h"
#include "libcomplexop.h"
#include "csl/initSanitizer.h"

void writeWilsonCoefficients(const std::string& coefficientName, 
                             std::complex<double> value, 
                             double Q_match, 
                             const std::string& fileName = "SM_wilson.csv");


namespace c7_sm {

void readParams(std::ifstream& inputFile, std::map<std::string, csl::InitSanitizer<real_t>*>& real, std::map<std::string, csl::InitSanitizer<complex_t>*>& complex);

}