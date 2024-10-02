#include <string>
#include <complex>

void writeWilsonCoefficients(const std::string& coefficientName, 
                             std::complex<double> value, 
                             double Q_match, 
                             const std::string& fileName = "SM_wilson.csv");