#include "json_parser.h"
#include "config.hpp"

int main() {

    std::string root_path = project_root.data();
    std::string filename = root_path +"Test/InputFiles/data_test.json";
    std::vector<Value> values;
    std::vector<Correlation> correlations;


    read_json(filename, values, correlations);


    return 0;
}
