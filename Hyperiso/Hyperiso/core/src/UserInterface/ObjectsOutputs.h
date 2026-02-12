#ifndef OBJECTS_OUTPUTS_H
#define OBJECTS_OUTPUTS_H

#include <string>
#include <vector>
#include <unordered_map>
#include <variant>

using Value = std::variant<std::monostate, double, int64_t, bool, std::string>;

struct ScanVar {
    std::string name;  
    std::string block;
    int pdg = 0;
    double min = 0;
    double max = 0;
    double step = 0;
};

struct DataPoint {
    std::vector<double> x;

    std::unordered_map<std::string, Value> y;
};

struct DataSet {
    std::unordered_map<std::string, Value> meta;
    std::vector<ScanVar> vars;

    std::vector<std::string> outputs_schema;

    std::vector<DataPoint> points;
};

enum class OutputFormat { TERMINAL, CSV, JSON };

struct OutputSpec {
    std::vector<std::string> keys;

    bool pretty_json = true;
    bool csv_write_header = true;
    char csv_sep = ',';
};

#endif