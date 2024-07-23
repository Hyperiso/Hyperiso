#ifndef JSON_PARSER_H
#define JSON_PARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// Structure pour les valeurs expérimentales
struct Value {
    std::string name;
    double central_value;
    double stat_error;
    double syst_error;
};

// Structure pour les corrélations
struct Correlation {
    std::string name1;
    std::string name2;
    double value;
};

// Fonctions de parsing JSON
std::string trim(const std::string& str);
std::string remove_quotes(const std::string& str);
std::vector<std::string> split(const std::string& str, char delimiter);
void read_json(const std::string& filename, std::vector<Value>& values, std::vector<Correlation>& correlations);

#endif // JSON_PARSER_H
