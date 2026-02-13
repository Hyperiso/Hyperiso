#include "ObservableHandler.h"
#include "YamlInputReader.h"
#include "UserParameterProxy.h"
#include "ObjectsOutputs.h"
#include "ScanRunner.h"
#include "MakeWriter.h"
#include "ObservableExtractor.h"

static void print_observable_usage() {
    std::cout << "Usage: ./main observable [options]\n"
              << "\nOptions:\n"
              << "  --model/-m <model_name>                    : Specify the model (SM, THDM, MSSM, ...)\n"
              << "  --observables/-os <observable_name>             : Specify the Wilson coefficient (e.g., C1, C2, ..., CQ1)\n"
              << "  --decay/-d <group_name>                    : Specify the group if multiple coefficients are provided\n"
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

        parser.addArgument(ArgumentBuilder()
                               .setLongName("observables")
                               .setShortName("os")
                               .setHelpText("List of observable names")
                               .setType(ArgType::STRING)
                               .setAllowsMultiple(true)
                               .build());

        parser.addArgument(ArgumentBuilder()
                               .setLongName("decay")
                               .setShortName("d")
                               .setHelpText("Decay name to get all the observables concerned (use if no observables have been provided)")
                               .setType(ArgType::STRING)
                               .setAllowsMultiple(false)
                               .setRequired(false)
                               .build());

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
                .setLongName("qcd_order")
                .setShortName("q")
                .setHelpText("Specify the calculation order (LO|NLO|NNLO)")
                .setType(ArgType::STRING)
                .setRequired(false)
                .setDefaultValue("LO")
                .addValidator(std::make_shared<AllowedValuesValidator>(
                    std::vector<std::string>{"LO","NLO","NNLO"}, ArgType::STRING
                ))
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

        parser.addArgument(ArgumentBuilder()
                            .setLongName("help")
                            .setShortName("h")
                            .setHelpText("Show usage information")
                            .setType(ArgType::STRING)
                            .setRequired(false)
                            .build());

        parser.parse(argc, argv);

        const auto pos = parser.getPositionalValues();
        const std::string cmd = pos.empty() ? "" : pos[0];


        bool wantHelp = false;
        try {
            wantHelp = (parser.getValue("help") == "true" || parser.getValue("help") == "1");
        } catch (...) {
        }

        if (wantHelp) {
            print_observable_usage();
            if (!cmd.empty()) {
                std::cerr << "\nUnknown subcommand: " << cmd << "\n";
                return 1;
            }
            return 0;
        }

        std::cout << "[observable] parsed options:\n";

        try { std::cout << "  model         = " << parser.getValue("model") << "\n"; } catch (...) {}
        try {
            auto ws = parser.getValues("observable");
            std::cout << "  observable        = ";
            for (size_t i = 0; i < ws.size(); ++i) std::cout << ws[i] << (i + 1 < ws.size() ? ", " : "");
            std::cout << "\n";
        } catch (...) {std::cout << "error" << std::endl;}

        try { std::cout << "  decay         = " << parser.getValue("decay") << "\n"; } catch (...) {std::cout << std::endl;}
        try { std::cout << "  observable         = " << parser.getValue("observable") << "\n"; } catch (...) {std::cout << std::endl;}
        try { std::cout << "  martypath     = " << parser.getValue("martypath") << "\n"; } catch (...) {std::cout << std::endl;}
        try { std::cout << "  input_file    = " << parser.getValue("input_file") << "\n"; } catch (...) {std::cout << std::endl;}
        try { std::cout << "  order         = " << parser.getValue("qcd_order") << "\n"; } catch (...) {std::cout << std::endl;}
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
        
        ObservableInterface oi;

        auto obs_str = parser.getValues("observables");
        auto order = OrderMapper::enum_elt(parser.getValue("qcd_order"));
        std::map<ObservableId, QCDOrder> obs;

        for (auto elem : obs_str) {
            obs[ObservableMapper::id_of(elem)] = order;
            
        }
        oi.add_observables(obs);

        std::string scan_yaml_path = parser.getValue("scan");

        std::vector<YamlScanParam> scan_params;
        if (scan_yaml_path != "None") {
            scan_params = YamlInputReader(scan_yaml_path).get_scan_params();
        }
        std::cout << scan_yaml_path << " hey : " << std::endl;

        if (scan_params.size() == 0) {
            for (auto elem : oi.get_current_observables()) {
                auto val = oi.compute_observable(elem);

                if (val.size() < 2) {
                    std::cout << "Observable " << elem.str() << " = " << val[0].value << std::endl;
                } else {
                    for (auto v : val) {
                        std::cout << "Observable " << elem.str() << " bin["<< v.bin.value().first << "," << v.bin.value().second << "] = " << v.value << std::endl;
                    }
                }
            }
        } else {
            std::cout << "doing scan now.." << std::endl;
            std::shared_ptr<UserParameterProxy> upp;
            if (config.model == Model::SM) {
                upp = std::make_shared<UserParameterProxy>(std::vector<ParameterType>{ParameterType::SM, ParameterType::WILSON});
            } else {
                upp = std::make_shared<UserParameterProxy>(std::vector<ParameterType>{ParameterType::SM, ParameterType::BSM, ParameterType::WILSON});
            }

            auto order = OrderMapper::enum_elt(parser.getValue("qcd_order"));

            auto wilsonExtractor = std::make_shared<ObservableExtractor>(oi, order);

            OutputSpec spec2;
            spec2.csv_write_header = true;

            std::vector<YamlScanParam> scan_params = YamlInputReader("scan.yml").get_scan_params();


            ScanRunner runner(*upp, scan_params, wilsonExtractor);

            DataSet ds2 = runner.run(spec2);

            ds2.meta["model"] = std::string("SM"); 
            // ds2.meta["Q_match"] = wbc.matching_scale;
            // ds2.meta["Q"] = wbc.hadronic_scale;
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
        print_observable_usage();
        return 1;
    }

    return 0;
}
