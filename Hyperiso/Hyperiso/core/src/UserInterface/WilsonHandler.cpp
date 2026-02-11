#include "WilsonHandler.h"
#include <iostream>
#include <string>
#include <complex>
#include <map>
#include <memory>
#include <vector>
#include "MartyWilson.h"
#include "ArgsParser.h"
#include "WilsonInterface.h"

static void print_wilson_usage() {
    std::cout << "Usage: ./main wilson [options]\n"
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

int handleWilsonOptions(int argc, char* argv[]) {
    try {
        ArgParser parser;

        // Sous-commande positionnelle: "wilson"
        // parser.addArgument(
        //     ArgumentBuilder()
        //         .setLongName("command")
        //         .setHelpText("Subcommand (wilson)")
        //         .setPositional(true)
        //         .setRequired(true)
        //         .build()
        // );

        // --model/-m
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

        // --wilson/-w (multi)
        parser.addArgument(
            ArgumentBuilder()
                .setLongName("wilson")
                .setShortName("w")
                .setHelpText("Specify Wilson coefficient(s) (e.g., C1, C2, ..., CQ1)")
                .setType(ArgType::STRING)
                .setAllowsMultiple(true)
                .setRequired(false)
                .build()
        );

        // --group/-g
        parser.addArgument(
            ArgumentBuilder()
                .setLongName("group")
                .setShortName("g")
                .setHelpText("Group name if multiple coefficients are provided")
                .setType(ArgType::STRING)
                .setRequired(false)
                .build()
        );

        // --Q_match/-q
        parser.addArgument(
            ArgumentBuilder()
                .setLongName("Q_match")
                .setShortName("q")
                .setHelpText("Set Q_match scale")
                .setType(ArgType::DOUBLE)
                .setRequired(false)
                .setDefaultValue("81")
                .build()
        );

        // --Q/-Q (oui shortName="Q" est OK)
        parser.addArgument(
            ArgumentBuilder()
                .setLongName("Q")
                .setShortName("Q")
                .setHelpText("Set Q scale")
                .setType(ArgType::DOUBLE)
                .setRequired(false)
                .setDefaultValue("4.18")
                .build()
        );

        // --martypath/-M (default None)
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

        // --input_file/-if
        parser.addArgument(
            ArgumentBuilder()
                .setLongName("input_file")
                .setShortName("if") // ton parser supporte les shortName multi-caractères
                .setHelpText("Input file for parameters spectrum")
                .setType(ArgType::STRING)
                .setRequired(false)
                .build()
        );

        // --order/-o (default LO) avec AllowedValuesValidator
        parser.addArgument(
            ArgumentBuilder()
                .setLongName("order")
                .setShortName("o")
                .setHelpText("Specify the calculation order (LO|NLO|NNLO)")
                .setType(ArgType::STRING)
                .setRequired(false)
                .setDefaultValue("LO")
                .addValidator(std::make_shared<AllowedValuesValidator>(
                    std::vector<std::string>{"LO","NLO","NNLO"}, ArgType::STRING
                ))
                .build()
        );

        // --write_to_flha/-F (default None)
        parser.addArgument(
            ArgumentBuilder()
                .setLongName("write_to_flha")
                .setShortName("F")
                .setHelpText("Write output to given file path or None")
                .setType(ArgType::STRING)
                .setRequired(false)
                .setDefaultValue("None")
                .build()
        );

        // --help/-h
        parser.addArgument(
            ArgumentBuilder()
                .setLongName("help")
                .setShortName("h")
                .setHelpText("Display help")
                .setType(ArgType::STRING)     // ton parser ne gère pas encore les flags sans valeur
                .setRequired(false)
                .setDefaultValue("false")     // astuce: si l'utilisateur met --help true
                .build()
        );

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
            print_wilson_usage();
            if (!cmd.empty()) {
                std::cerr << "\nUnknown subcommand: " << cmd << "\n";
                return 1;
            }
            return 0;
        }

        // ---- Ici tu mettras le fonctionnement réel plus tard ----
        // Pour l'instant, on montre juste ce qui a été parsé.

        std::cout << "[wilson] parsed options:\n";

        try { std::cout << "  model         = " << parser.getValue("model") << "\n"; } catch (...) {}
        try {
            auto ws = parser.getValues("wilson");
            std::cout << "  wilson        = ";
            for (size_t i = 0; i < ws.size(); ++i) std::cout << ws[i] << (i + 1 < ws.size() ? ", " : "");
            std::cout << "\n";
        } catch (...) {}

        try { std::cout << "  wilson         = " << parser.getValue("wilson") << "\n"; } catch (...) {std::cout << std::endl;}
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
        
        WilsonInterface wi;
        WilsonBuildConfig wbc;

        wbc.matching_scale = parser.get<double>("Q_match");
        wbc.hadronic_scale = parser.get<double>("Q");
        wbc.order = OrderMapper::enum_elt(parser.getValue("order"));

        std::unordered_set<WCoefId> coefs;

        if (parser.exists("wilson")) {
            for (auto w : parser.getValues("wilson")) {
                coefs.emplace(WCoefMapper::id_of(w));
                auto w_enum = WCoefMapper::id_of(w);
                auto group = WCoefMapper::group_of(w_enum);
                wbc.groups.emplace(GroupMapper::to_id(group));
            }
        } else if (parser.exists("group")) {
            WGroupId grp = GroupMapper::id_of(parser.getValue("group"));
            wbc.groups.emplace(grp);
            auto coefs_temp = WCoefMapper::get_group(GroupMapper::enum_of(grp).value());
            for (auto coef : coefs_temp) {
                coefs.emplace(WCoefMapper::to_id(coef));
            }
        } else {
            std::cout << "error" << std::endl;
        }

        wi.build(wbc);
        std::cout << std::endl;
        
        for (auto coef : coefs) {
            std::cout << "SM Wilson Coefficient " << WCoefMapper::str(coef) << " at LO : "<< wi.getM(WCoefMapper::group_of(coef), WCoefMapper::enum_of(coef).value(), QCDOrder::LO, ContributionType::SM) << std::endl;
            std::cout << "BSM Wilson Coefficient " << WCoefMapper::str(coef) << " at LO : "<< wi.getM(WCoefMapper::group_of(coef), WCoefMapper::enum_of(coef).value(), QCDOrder::LO, ContributionType::BSM) << std::endl;

            if (wbc.order > QCDOrder::LO) {
                std::cout << "SM Wilson Coefficient " << WCoefMapper::str(coef) << " at NLO : "<< wi.getM(WCoefMapper::group_of(coef), WCoefMapper::enum_of(coef).value(), QCDOrder::NLO, ContributionType::SM) << std::endl;
                std::cout << "BSM Wilson Coefficient " << WCoefMapper::str(coef) << " at NLO : "<< wi.getM(WCoefMapper::group_of(coef), WCoefMapper::enum_of(coef).value(), QCDOrder::NLO, ContributionType::BSM) << std::endl;
                if (wbc.order > QCDOrder::NLO) {
                    std::cout << "SM Wilson Coefficient " << WCoefMapper::str(coef) << " at NNLO : "<< wi.getM(WCoefMapper::group_of(coef), WCoefMapper::enum_of(coef).value(), QCDOrder::NNLO, ContributionType::SM) << std::endl;
                    std::cout << "BSM Wilson Coefficient " << WCoefMapper::str(coef) << " at NNLO : "<< wi.getM(WCoefMapper::group_of(coef), WCoefMapper::enum_of(coef).value(), QCDOrder::NNLO, ContributionType::BSM) << std::endl;
                }
            }
            std::cout << std::endl;
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n\n";
        print_wilson_usage();
        return 1;
    }


}