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

        data[rowIndex][realIndex] = std::to_string(value.real());
        data[rowIndex][imgIndex] = std::to_string(value.imag());
    } else {
        std::vector<std::string> newRow(headers.size(), "NaN");
        
        newRow[0] = std::to_string(Q_match); 
        
        int realIndex = std::distance(headers.begin(), std::find(headers.begin(), headers.end(), coefficientName + "_real"));
        int imgIndex = std::distance(headers.begin(), std::find(headers.begin(), headers.end(), coefficientName + "_img"));

        newRow[realIndex] = std::to_string(value.real());
        newRow[imgIndex] = std::to_string(value.imag());

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


// int main() {
//     std::string C1 = "C1";
//     std::string C7 = "C7";
//     std::string C9 = "C9";

//     std::complex<double> valueC1_1(1.1, 0.1);
//     std::complex<double> valueC1_2(1.2, 0.2);

//     std::complex<double> valueC7_1(7.1, -0.7);
//     std::complex<double> valueC7_2(7.2, -0.8);

//     std::complex<double> valueC9_1(9.1, 0.9); 

//     std::cout << "Test 1: Ajout du coefficient C1 à Q_match = 100\n";
//     writeWilsonCoefficients(C1, valueC1_1, 100);

//     std::cout << "Test 2: Ajout du coefficient C1 à Q_match = 200\n";
//     writeWilsonCoefficients(C1, valueC1_2, 200);

//     std::cout << "Test 3: Ajout du coefficient C7 à Q_match = 100\n";
//     writeWilsonCoefficients(C7, valueC7_1, 100);

//     std::cout << "Test 4: Ajout du coefficient C7 à Q_match = 300\n";
//     writeWilsonCoefficients(C7, valueC7_2, 300);

//     std::cout << "Test 5: Ajout du coefficient C9 à Q_match = 150\n";
//     writeWilsonCoefficients(C9, valueC9_1, 150);

//     writeWilsonCoefficients("C10", {42,42}, 500);
//     std::cout << "Test terminé. Consultez le fichier SM_wilson.csv pour vérifier les résultats.\n";

//     return 0;
// }