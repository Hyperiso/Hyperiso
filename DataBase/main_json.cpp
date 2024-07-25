#include "json_parser.h"

int main() {
    std::string filename = "data_exp.json";
    std::vector<Value> values;
    std::vector<Correlation> correlations;
    // Lire les données du fichier JSON
    read_json(filename, values, correlations);


    return 0;
}
