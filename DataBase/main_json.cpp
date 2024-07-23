#include "json_parser.h"

int main() {
    std::string filename = "data.json";
    std::vector<Value> values;
    std::vector<Correlation> correlations;
    // Lire les données du fichier JSON
    read_json(filename, values, correlations);

    // // Afficher les valeurs lues
    // std::cout << "Values:" << std::endl;
    // for (const auto& val : values) {
    //     std::cout << "Name: " << val.name
    //               << ", Central Value: " << val.central_value
    //               << ", Stat Error: " << val.stat_error
    //               << ", Syst Error: " << val.syst_error
    //               << std::endl;
    // }

    // // Afficher les corrélations lues
    // std::cout << "Correlations:" << std::endl;
    // for (const auto& corr : correlations) {
    //     std::cout << "Name1: " << corr.name1
    //               << ", Name2: " << corr.name2
    //               << ", Value: " << corr.value
    //               << std::endl;
    // }

    return 0;
}
