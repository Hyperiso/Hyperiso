// #include "json_parser.h"
#include <cassert>

int main() {
    // assert(trim("   test   ") == "test");
    // assert(trim("test") == "test");
    // assert(trim("   ") == "");

    // assert(remove_quotes("\"test\"") == "test");
    // assert(remove_quotes("\"test,") == "test");
    // assert(remove_quotes("test") == "test");

    // auto tokens = split("a,b,c", ',');
    // assert(tokens.size() == 3);
    // assert(tokens[0] == "a");
    // assert(tokens[1] == "b");
    // assert(tokens[2] == "c");

    // std::ofstream tempFile("test_data.json");
    // tempFile << R"({
    //     "values": [
    //         {
    //             "name": "value1",
    //             "central_value": 10.5,
    //             "stat_error": 0.5,
    //             "syst_error": 1.0
    //         },
    //         {
    //             "name": "value2",
    //             "central_value": 20.0,
    //             "stat_error": 1.0,
    //             "syst_error": 2.0
    //         }
    //     ],
    //     "correlations": [
    //         {
    //             "name1": "value1",
    //             "name2": "value2",
    //             "value": 0.8
    //         }
    //     ]
    // })";
    // tempFile.close();

    // std::vector<Value> values;
    // std::vector<Correlation> correlations;

    // read_json("test_data.json", values, correlations);

    // assert(values.size() == 2);
    // assert(values[0].name == "value1");
    // assert(values[0].central_value == 10.5);
    // assert(values[0].stat_error == 0.5);
    // assert(values[0].syst_error == 1.0);

    // assert(values[1].name == "value2");
    // assert(values[1].central_value == 20.0);
    // assert(values[1].stat_error == 1.0);
    // assert(values[1].syst_error == 2.0);

    // assert(correlations.size() == 1);
    // assert(correlations[0].name1 == "value1");
    // assert(correlations[0].name2 == "value2");
    // assert(correlations[0].value == 0.8);

    // read_json("invalid_file.json", values, correlations);
    // assert(values.size() == 2);
    // assert(correlations.size() == 1);

    // std::cout << "All test succeeded." << std::endl;

    // std::remove("test_data.json");

    return 0;
}
