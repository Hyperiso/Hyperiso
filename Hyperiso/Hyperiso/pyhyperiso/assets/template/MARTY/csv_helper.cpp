#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <complex>
#include <algorithm>
#include "csv_helper.h"

std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

void writeWilsonCoefficients(const std::string& coefficientName, 
                             std::complex<double> value, 
                             double Q_match, 
                             const std::string& fileName) {
    
    std::ifstream file(fileName);
    std::vector<std::vector<std::string>> data;
    std::vector<std::string> headers;
    std::map<double, int> qMatchRowMap;
    bool coefficientExists = false;
    bool fileExists = file.good();
    
    if (fileExists && file.is_open()) {
        std::string line;
        bool isFirstLine = true;

        while (std::getline(file, line)) {
            auto tokens = split(line, ',');
            if (isFirstLine) {
                headers = tokens;
                isFirstLine = false;
            } else {
                data.push_back(tokens);
                qMatchRowMap[std::stod(tokens[0])] = data.size() - 1;
            }
        }

        for (const auto& header : headers) {
            if (header == coefficientName + "_real" || header == coefficientName + "_img") {
                coefficientExists = true;
                break;
            }
        }

        file.close();
    } else {
        headers.push_back("Q_match");
    }

    if (!coefficientExists) {
        headers.push_back(coefficientName + "_real");
        headers.push_back(coefficientName + "_img");

        for (auto& row : data) {
            row.push_back("NaN");
            row.push_back("NaN");
        }
    }

    if (qMatchRowMap.find(Q_match) != qMatchRowMap.end()) {
        int rowIndex = qMatchRowMap[Q_match];
        int realIndex = std::distance(headers.begin(), std::find(headers.begin(), headers.end(), coefficientName + "_real"));
        int imgIndex = std::distance(headers.begin(), std::find(headers.begin(), headers.end(), coefficientName + "_img"));

        std::stringstream ss_r;
        ss_r << std::scientific << value.real();
        data[rowIndex][realIndex] = ss_r.str();
        
        std::stringstream ss_i;
        ss_i << std::scientific << value.imag();
        data[rowIndex][imgIndex] = ss_i.str();
    } else {
        std::vector<std::string> newRow(headers.size(), "NaN");
        
        newRow[0] = std::to_string(Q_match); 
        
        int realIndex = std::distance(headers.begin(), std::find(headers.begin(), headers.end(), coefficientName + "_real"));
        int imgIndex = std::distance(headers.begin(), std::find(headers.begin(), headers.end(), coefficientName + "_img"));

        std::stringstream ss_r;
        ss_r << std::scientific << value.real();
        newRow[realIndex] = ss_r.str();
        
        std::stringstream ss_i;
        ss_i << std::scientific << value.imag();
        newRow[imgIndex] = ss_i.str();

        data.push_back(newRow);
    }

    std::ofstream outFile(fileName);
    
    for (size_t i = 0; i < headers.size(); ++i) {
        outFile << headers[i];
        if (i < headers.size() - 1) {
            outFile << ",";
        }
    }
    outFile << "\n";

    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            outFile << row[i];
            if (i < row.size() - 1) {
                outFile << ",";
            }
        }
        outFile << "\n";
    }

    outFile.close();
}


void readParams(std::ifstream& inputFile, 
                std::map<std::string, csl::InitSanitizer<real_t>*>& real, 
                std::map<std::string, csl::InitSanitizer<complex_t>*>& complex) {

    std::string line;
    std::map<std::string, std::complex<double>> complex_;

    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        std::string key;
        double value;

        if (std::getline(iss, key, ',') && iss >> value) {
            if (key.size() > 4 && (key.substr(key.size() - 4) == "_rel" || key.substr(key.size() - 4) == "_img")) {
                std::string baseKey = key.substr(0, key.size() - 4);
                if (complex.find(baseKey) != complex.end()) {
                    if (key.substr(key.size() - 4) == "_rel") {
                        if (complex_.find(baseKey) != complex_.end()) {
                            complex_[baseKey] += std::complex<double>(value, 0);
                            *complex[baseKey] = complex_[baseKey];
                        } else {
                            complex_[baseKey] = std::complex<double>(value, 0);
                        }
                    } else if (key.substr(key.size() - 4) == "_img") {
                        if (complex_.find(baseKey) != complex_.end()) {
                            complex_[baseKey] += std::complex<double>(0, value);
                            *complex[baseKey] = complex_[baseKey];
                        } else {
                            complex_[baseKey] = std::complex<double>(0, value);
                        }
                    }
                } else {
                    std::cerr << "Warning: Complex parameter " << baseKey << " not found in the map." << std::endl;
                }

            } else {
                std::cout << "key : " << key << std::endl;
                std::cout << "value : " << value << std::endl;
                if (real.find(key) != real.end()) {
                    *real[key] = value;
                } else {
                    std::cerr << "Warning: Real parameter " << key << " not found in the map." << std::endl;
                }
            }
        }
    }
}
