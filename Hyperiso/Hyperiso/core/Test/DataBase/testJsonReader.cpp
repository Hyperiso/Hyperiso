#include "json_parser.h"
#include "config.hpp"
#include <iostream>
#include <cassert>

int main() {
    std::string root_path = project_root.data();
    std::string filename = root_path + "Test/InputFiles/data_test.json";
    
    std::vector<Value> values;
    std::vector<Correlation> correlations;

    read_json(filename, values, correlations);

    assert(values.size() == 2);
    assert(values[0].name == "BRuntag_Bsmumu");
    assert(values[0].central_value == 2.65e-09);
    assert(values[0].stat_error == 4.3e-10);
    assert(values[0].syst_error == 0.0);

    assert(values[1].name == "BR_Bdmumu");
    assert(values[1].central_value == 1.09e-10);
    assert(values[1].stat_error == 7.4e-11);
    assert(values[1].syst_error == 0.0);

    assert(correlations.size() == 1);
    assert(correlations[0].name1 == "BRuntag_Bsmumu");
    assert(correlations[0].name2 == "BR_Bdmumu");
    assert(correlations[0].value == -0.2);

    std::cout << "All tests passed" << std::endl;

    return 0;
}
