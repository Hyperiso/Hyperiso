#include "ObservableHandler.h"

static void print_observable_usage() {
    std::cout << "Usage: ./main observable [options]\n"
              << "\nOptions:\n"
              << "  --model/-m <model_name>                    : Specify the model (SM, THDM, MSSM, ...)\n"
              << "  --wilson/-w <coefficient_name>             : Specify the Wilson coefficient (e.g., C1, C2, ..., CQ1)\n"
              << "  --group/-g <group_name>                    : Specify the group if multiple coefficients are provided\n"
              << "  --Q_match/-q <value>                       : Set Q_match scale\n"
              << "  --Q/-Q <value>                             : Set Q scale\n"
              << "  --martypath/-M <marty_model_path>          : Marty Model path (default: None)\n"
              << "  --input_file/-if <slha_name>               : input file for parameters spectrum\n"
              << "  --order/-o <LO|NLO|NNLO>                   : Specify the calculation order (default: LO)\n"
              << "  --write_to_flha/-F <path|None>             : Write output to file (default: None)\n"
              << "  --help/-h                                  : Display this help message\n";
}

//TODO : refactor
int handleObservableOptions(int argc, char* argv[]) {
    try {
        ArgParser parser;

        // parser.addArgument(ArgumentBuilder()
        //                        .setLongName("observables")
        //                        .setShortName("os")
        //                        .setHelpText("Comma-separated list of observable names")
        //                        .setType(ArgType::STRING)
        //                        .setAllowsMultiple(true)
        //                        .build());

        parser.addArgument(ArgumentBuilder()
                            .setLongName("input_file")
                            .setShortName("if")
                            .setHelpText("Input file path")
                            .setType(ArgType::STRING)
                            .setDefaultValue("Test/InputFiles/testinput_thdm.lha")
                            .setRequired(false)
                            .build());

        parser.addArgument(
                ArgumentBuilder()
                    .setLongName("model")
                    .setShortName("m")
                    .setHelpText("Specify the model (SM, THDM, MSSM, ...)")
                    .setType(ArgType::STRING)
                    .setRequired(false)
                    .setDefaultValue("SM")
                    .build()
            );

        parser.addArgument(
                ArgumentBuilder()
                    .setLongName("martypath")
                    .setShortName("M")
                    .setHelpText("Marty Model path")
                    .setType(ArgType::STRING)
                    .setRequired(false)
                    .setDefaultValue("None")
                    .build()
            );

        parser.addArgument(
                ArgumentBuilder()
                    .setLongName("model")
                    .setShortName("m")
                    .setHelpText("Specify the model (SM, THDM, MSSM, ...)")
                    .setType(ArgType::STRING)
                    .setRequired(false)
                    .setDefaultValue("SM")
                    .build()
            );

        parser.addArgument(
                ArgumentBuilder()
                    .setLongName("model")
                    .setShortName("m")
                    .setHelpText("Specify the model (SM, THDM, MSSM, ...)")
                    .setType(ArgType::STRING)
                    .setRequired(false)
                    .setDefaultValue("SM")
                    .build()
            );

        parser.addArgument(ArgumentBuilder()
                            .setLongName("help")
                            .setShortName("h")
                            .setHelpText("Show usage information")
                            .setType(ArgType::STRING)
                            .setRequired(false)
                            .build());

        parser.parse(argc, argv);

        // Vérifier la sous-commande
        const auto pos = parser.getPositionalValues();
        const std::string cmd = pos.empty() ? "" : pos[0];

        // Gestion help simple
        // (Avec ton parser actuel, il faut passer une valeur: --help true)
        bool wantHelp = false;
        try {
            wantHelp = (parser.getValue("help") == "true" || parser.getValue("help") == "1");
        } catch (...) {
            // ignore si absent
        }

        if (wantHelp) {
            print_observable_usage();
            if (!cmd.empty()) {
                std::cerr << "\nUnknown subcommand: " << cmd << "\n";
                return 1;
            }
            return 0;
        }

        // ---- Ici tu mettras le fonctionnement réel plus tard ----
        // Pour l'instant, on montre juste ce qui a été parsé.

        std::cout << "[observable] parsed options:\n";

        try { std::cout << "  model         = " << parser.getValue("model") << "\n"; } catch (...) {}
        try {
            auto ws = parser.getValues("wilson");
            std::cout << "  observable        = ";
            for (size_t i = 0; i < ws.size(); ++i) std::cout << ws[i] << (i + 1 < ws.size() ? ", " : "");
            std::cout << "\n";
        } catch (...) {}

        try { std::cout << "  observable         = " << parser.getValue("observable") << "\n"; } catch (...) {std::cout << std::endl;}
        try { std::cout << "  group         = " << parser.getValue("group") << "\n"; } catch (...) {std::cout << std::endl;}
        try { std::cout << "  Q_match       = " << parser.getValue("Q_match") << "\n"; } catch (...) {std::cout << std::endl;}
        try { std::cout << "  Q             = " << parser.getValue("Q") << "\n"; } catch (...) {std::cout << std::endl;}
        try { std::cout << "  martypath     = " << parser.getValue("martypath") << "\n"; } catch (...) {std::cout << std::endl;}
        try { std::cout << "  input_file    = " << parser.getValue("input_file") << "\n"; } catch (...) {std::cout << std::endl;}
        try { std::cout << "  order         = " << parser.getValue("order") << "\n"; } catch (...) {std::cout << std::endl;}
        try { std::cout << "  write_to_flha = " << parser.getValue("write_to_flha") << "\n"; } catch (...) {std::cout << std::endl;}


        HyperisoMaster hyp = HyperisoMaster();

        HyperisoConfig config;

        if (parser.getValue("martypath") == "None") {
            std::cout << parser.getValue("model") << std::endl;
            if (parser.getValue("model") == "SM") {
                config.model = Model::SM;
            } else if (parser.getValue("model") == "THDM") {
                config.model = Model::THDM;
            } else if (parser.getValue("model") == "MSSM") {
                config.model = Model::SUSY;
            } else if (parser.getValue("model") == "NMSSM") {
                config.model = Model::SUSY;
            } else {
                std::cout << "error" << std::endl;
            }
        } else { 
            config.model = Model::MARTY;
        }
        std::string lha_path = parser.getValue("input_file");
        
        hyp.init(lha_path, config);
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n\n";
        print_observable_usage();
        return 1;
    }

    return 0;
}
