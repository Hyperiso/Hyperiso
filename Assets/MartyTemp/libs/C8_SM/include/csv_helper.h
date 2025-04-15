#include <string>
#include <complex>
#include <fstream>
#include <map>
#include "common.h"
#include "libcomplexop.h"
#include "csl/initSanitizer.h"

using real_t = double;
using complex_t = std::complex<double>;
void writeWilsonCoefficients(const std::string& coefficientName, 
                             std::complex<double> value, 
                             double Q_match, 
                             const std::string& fileName = "SM_wilson.csv");



void readParams(std::ifstream& inputFile, std::map<std::string, csl::InitSanitizer<real_t>*>& real, std::map<std::string, csl::InitSanitizer<complex_t>*>& complex);
