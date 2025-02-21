#ifndef COMPARECSV_H
#define COMPARECSV_H
#include "Logger.h"
#include <string>
#include <vector>

bool compareCSV(const std::string& file1, const std::string& file2, double tolerance);
std::vector<std::vector<std::string>> readCSV(const std::string& fileName);
bool compareValues(double val1, double val2, double tolerance);

#endif // COMPARE_CSV_H