#include "CompareCsv.h"
#include <fstream>
#include <sstream>

#include <cmath>
#include <iostream>


bool compareValues(double val1, double val2, double tolerance) {
    if (val2 == 0) {
        return std::fabs(val1) <= tolerance;
    }
    return std::fabs((val1 - val2) / val2) <= tolerance;
}

std::vector<std::vector<std::string>> readCSV(const std::string& fileName) {
    std::vector<std::vector<std::string>> data;
    std::ifstream file(fileName);
    if (!file.is_open()) {
        LOG_ERROR("Could not open file: ", fileName);
    }

    std::string line, cell;
    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::vector<std::string> row;
        while (std::getline(lineStream, cell, ',')) {
            row.push_back(cell);
        }
        data.push_back(row);
    }
    return data;
}

bool compareCSV(const std::string& file1, const std::string& file2, double tolerance) {
    auto data1 = readCSV(file1);
    auto data2 = readCSV(file2);
    LOG_INFO(file1);
    if (data1.size() != data2.size() || data1[0].size() != data2[0].size()) {
        std::cerr << "File sizes are different." << std::endl;
        return false;
    }
    LOG_INFO("TOUT VA BIEEEEN");
    for (size_t i = 1; i < data1.size(); ++i) { // start from 1 to skip headers
        for (size_t j = 0; j < data1[i].size(); ++j) {
            double val1 = std::stod(data1[i][j]);
            double val2 = std::stod(data2[i][j]);
            if (!compareValues(val1, val1, tolerance)) {
                std::cerr << "Mismatch at row " << i << ", column " << j
                          << ": " << data1[i][j] << " != " << data2[i][j] << " with tolerance " << tolerance << std::endl;
                return false;
            }
        }
    }
    return true;
}
