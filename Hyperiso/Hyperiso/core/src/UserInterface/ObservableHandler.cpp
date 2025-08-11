#include "ObservableHandler.h"


//TODO : refactor
int handleObservableOptions(int argc, char* argv[]) {
    return 0;
    // ArgParser parser;

    // parser.addArgument(ArgumentBuilder()
    //                        .setLongName("observables")
    //                        .setShortName("os")
    //                        .setHelpText("Comma-separated list of observable names")
    //                        .setType(ArgType::STRING)
    //                        .setAllowsMultiple(true)
    //                        .build());

    // parser.addArgument(ArgumentBuilder()
    //                        .setLongName("input_file")
    //                        .setShortName("if")
    //                        .setHelpText("Input file path")
    //                        .setType(ArgType::STRING)
    //                        .setDefaultValue("Test/InputFiles/testinput_thdm.lha")
    //                        .setRequired(false)
    //                        .build());

    // parser.addArgument(ArgumentBuilder()
    //                        .setLongName("help")
    //                        .setShortName("h")
    //                        .setHelpText("Show usage information")
    //                        .setType(ArgType::STRING)
    //                        .setRequired(false)
    //                        .build());

    // try {
    //     parser.parse(argc, argv);

    //     // if (parser.getValue("help").empty() == false) {
    //     //     parser.displayHelp();
    //     //     return 0;
    //     // }

    //     std::string input_file = parser.getValue("input_file");
    //     std::vector<std::string> observable_names;

    //     if (parser.getValues("observables").size() > 0) {
    //         std::string obs_list = parser.getValue("observables");
    //         size_t pos = 0;
    //         while ((pos = obs_list.find(',')) != std::string::npos) {
    //             observable_names.push_back(obs_list.substr(0, pos));
    //             obs_list.erase(0, pos + 1);
    //         }
    //         observable_names.push_back(obs_list);
    //     }

    //     MemoryManager::GetInstance()->init(input_file, Model::SM);

    //     ObservableInterface obs_interface;

    //     std::vector<Observables> selected_obs;
    //     for (const auto& name : observable_names) {
    //         selected_obs.push_back(ObservableMapper::enum_elt(name));
    //     }

    //     for (const auto& obs : selected_obs) {
    //         std::cout << ObservableMapper::str(obs) << ": "
    //                   << obs_interface.compute_observable(obs) << "\n";
    //     }

    //     if (!selected_obs.empty()) {
    //         double chi2 = obs_interface.compute_chi2();
    //         std::cout << "Chi-squared: " << chi2 << "\n";
    //     }

    //     return 0;
    // } catch (const std::exception& ex) {
    //     std::cerr << "An error occurred: " << ex.what() << std::endl;
    //     return 1;
    // }
}
