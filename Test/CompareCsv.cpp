#include "compare_csv.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>

bool compareValues(double val1, double val2, double tolerance) {
    return std::fabs(val1 - val2) <= tolerance;
}

std::vector<std::vector<double>> readCSV(const std::string& fileName) {
    std::vector<std::vector<double>> data;
    std::ifstream file(fileName);
    std::string line, cell;
    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::vector<double> row;
        while (std::getline(lineStream, cell, ',')) {
            row.push_back(std::stod(cell));
        }
        data.push_back(row);
    }
    return data;
}

bool compareCSV(const std::string& file1, const std::string& file2, double tolerance) {
    auto data1 = readCSV(file1);
    auto data2 = readCSV(file2);

    if (data1.size() != data2.size() || data1[0].size() != data2[0].size()) {
        return false;
    }

    for (size_t i = 1; i < data1.size(); ++i) { // start from 1 to skip headers
        for (size_t j = 0; j < data1[i].size(); ++j) {
            if (!compareValues(data1[i][j], data2[i][j], tolerance)) {
                std::cerr << "Mismatch at row " << i << ", column " << j
                          << ": " << data1[i][j] << " != " << data2[i][j] << std::endl;
                return false;
            }
        }
    }
    return true;
}