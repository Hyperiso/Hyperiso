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
#include "YamlInputReader.h"
#include "UserParameterProxy.h"
#include "ObjectsOutputs.h"
#include "WilsonExtractor.h"
#include "CsvWriter.h"
#include "JsonWriter.h"
#include "ScanRunner.h"
#include "MakeWriter.h"

static void print_wilson_usage() {
    std::cout << "Usage: ./main wilson [options]\n"
              << "\nOptions:\n"
              << "  --model/-m <model_name>                    : Specify the model (SM, THDM, MSSM, ...)\n"
              << "  --wilson/-w <coefficient_name>             : Specify the Wilson coefficient (e.g., C1, C2, ..., CQ1)\n"
              << "  --group/-g <group_name>                    : Specify the group if multiple coefficients are provided\n"
              << "  --Q_match/-QM <value>                       : Set Q_match scale\n"
              << "  --Q/-Q <value>                             : Set Q scale\n"
              << "  --martypath/-M <marty_model_path>          : Marty Model path (default: None)\n"
              << "  --input_file/-if <slha_name>               : input file for parameters spectrum\n"
              << "  --order/-q <LO|NLO|NNLO>                   : Specify the calculation order (default: LO)\n"
              << "  --write_to_flha/-F <path|None>             : Write output to file (default: None)\n"
              << "  --help/-h                                  : Display this help message\n";
}

int handleWilsonOptions(int argc, char* argv[]) {
    try {
        ArgParser parser;

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
                .setShortName("QM")
                .setHelpText("Set Q_match scale")
                .setType(ArgType::DOUBLE)
                .setRequired(false)
                .setDefaultValue("81")
                .build()
        );

        // --Q/-Q
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
                .setShortName("if")
                .setHelpText("Input file for parameters spectrum")
                .setType(ArgType::STRING)
                .setRequired(false)
                .build()
        );

        // --order/-o (default LO)
        parser.addArgument(
            ArgumentBuilder()
                .setLongName("qcd_order")
                .setShortName("q")
                .setHelpText("Specify the calculation order in QCD (LO|NLO|NNLO)")
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

        parser.addArgument(
            ArgumentBuilder()
                .setLongName("scan")
                .setShortName("s")
                .setHelpText("Is the calculation a scan, if set, a path to an input file (.yaml/.yml) need to be provided, with the format min/max/step for each parameter.")
                .setType(ArgType::STRING)
                .setRequired(false)
                .setDefaultValue("None")
                .build()
        );


        parser.addArgument(
            ArgumentBuilder()
                .setLongName("output")
                .setShortName("o")
                .setHelpText("In which format you want the output (single point or scan). Options are terminal (default), csv, json, lha.")
                .setType(ArgType::STRING)
                .setRequired(false)
                .setDefaultValue("terminal")
                .build()
        );

        // --help/-h
        parser.addArgument(
            ArgumentBuilder()
                .setLongName("help")
                .setShortName("h")
                .setHelpText("Display help")
                .setType(ArgType::STRING)   
                .setRequired(false)
                .setDefaultValue("false")
                .build()
        );

        parser.parse(argc, argv);

        const auto pos = parser.getPositionalValues();
        const std::string cmd = pos.empty() ? "" : pos[0];


        bool wantHelp = false;
        try {
            wantHelp = (parser.getValue("help") == "true" || parser.getValue("help") == "1");
        } catch (...) {
        }

        if (wantHelp) {
            print_wilson_usage();
            if (!cmd.empty()) {
                std::cerr << "\nUnknown subcommand: " << cmd << "\n";
                return 1;
            }
            return 0;
        }

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
        try { std::cout << "  qcd_order         = " << parser.getValue("qcd_order") << "\n"; } catch (...) {std::cout << std::endl;}
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
        
        std::string scan_yaml_path = parser.getValue("scan");

        std::vector<YamlScanParam> scan_params;
        if (scan_yaml_path != "None") {
            scan_params = YamlInputReader(scan_yaml_path).get_scan_params();
        }
        std::cout << scan_yaml_path << " hey : " << std::endl;

        WilsonInterface wi;
        WilsonBuildConfig wbc;

        wbc.matching_scale = parser.get<double>("Q_match");
        wbc.hadronic_scale = parser.get<double>("Q");
        wbc.order = OrderMapper::enum_elt(parser.getValue("qcd_order"));

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
        if (scan_params.size() == 0) {
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
        } else {
            std::cout << "doing scan now.." << std::endl;
            std::shared_ptr<UserParameterProxy> upp;
            if (config.model == Model::SM) {
                upp = std::make_shared<UserParameterProxy>(std::vector<ParameterType>{ParameterType::SM, ParameterType::WILSON});
            } else {
                upp = std::make_shared<UserParameterProxy>(std::vector<ParameterType>{ParameterType::SM, ParameterType::BSM, ParameterType::WILSON});
            }

            auto wilsonExtractor = std::make_shared<WilsonCoeffExtractor>(wi, wbc, coefs);

            OutputSpec spec2;
            spec2.csv_write_header = true;

            std::vector<YamlScanParam> scan_params = YamlInputReader("scan.yml").get_scan_params();


            ScanRunner runner(*upp, scan_params, wilsonExtractor);

            DataSet ds2 = runner.run(spec2);

            ds2.meta["model"] = std::string("SM"); 
            ds2.meta["Q_match"] = wbc.matching_scale;
            ds2.meta["Q"] = wbc.hadronic_scale;
            ds2.meta["qcd_order"] = std::string("LO");

            std::string output_format = parser.getValue("output");

            if (output_format.ends_with("csv")) {
                auto w = make_writer(OutputFormat::CSV, output_format);
                w->write(ds2, spec2);
                std::cout << "Wrote output_format\n";
            } else if (output_format.ends_with("json")) {
                auto w = make_writer(OutputFormat::JSON, output_format);
                w->write(ds2, spec2);
                std::cout << "Wrote output_format\n";
            } else {
                std::cout << "bonjour" << std::endl;
            }

            
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n\n";
        print_wilson_usage();
        return 1;
    }


}