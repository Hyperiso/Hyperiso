#include <iostream>
#include <unordered_set>
#include <memory>

// inclue ton système d'output (celui qu'on a fait)
#include "CsvWriter.h"
#include "TerminalWriter.h"
#include "JsonWriter.h"
#include "MakeWriter.h"
#include "Include.h"
#include "WilsonInterface.h"
#include "WilsonExtractor.h"
#include "UserParameterProxy.h"
#include "ObjectsOutputs.h"
#include "ScanRunner.h"

int main() {
    // 1) dataset fake
    DataSet ds;
    ds.meta["model"] = std::string("SM");
    ds.meta["Q_match"] = 160.0;
    ds.meta["Q"] = 4.8;
    ds.meta["qcd_order"] = std::string("NLO");

    ds.vars = {
        ScanVar{ "MINPAR:3", "MINPAR", 3, 1.0, 4.0, 1.0 },
        ScanVar{ "EXTPAR:1", "EXTPAR", 1, 10.0, 13.0, 1.0 }
    };

    // colonnes outputs (ordre stable CSV)
    ds.outputs_schema = {
        "WC:C7:LO:SM",
        "WC:C7:LO:BSM",
        "OBS:BR_BtoXsGamma:SM",
        "OBS:BR_BtoXsGamma:SM+BSM"
    };

    // 2) quelques points
    {
        DataPoint p;
        p.x = {1.0, 10.0};
        p.y["WC:C7:LO:SM"] = 0.12;
        p.y["WC:C7:LO:BSM"] = 0.01;
        p.y["OBS:BR_BtoXsGamma:SM"] = 3.15e-4;
        p.y["OBS:BR_BtoXsGamma:SM+BSM"] = 3.22e-4;
        ds.points.push_back(std::move(p));
    }
    {
        DataPoint p;
        p.x = {2.0, 11.0};
        p.y["WC:C7:LO:SM"] = 0.11;
        p.y["WC:C7:LO:BSM"] = 0.02;
        p.y["OBS:BR_BtoXsGamma:SM"] = 3.14e-4;
        p.y["OBS:BR_BtoXsGamma:SM+BSM"] = 3.30e-4;
        ds.points.push_back(std::move(p));
    }

    // 3) écrire
    OutputSpec spec;
    spec.csv_write_header = true;
    spec.csv_sep = ',';
    // spec.keys = {...}; // si tu veux filtrer certains outputs uniquement

    // CSV
    {
        auto w = make_writer(OutputFormat::CSV, "out.csv");
        w->write(ds, spec);
        std::cout << "Wrote out.csv\n";
    }

    // JSON (via DBNode)
    {
        auto w = make_writer(OutputFormat::JSON, "out.json");
        w->write(ds, spec);
        std::cout << "Wrote out.json\n";
    }

    // Terminal (JSON-like)
    {
        auto w = make_writer(OutputFormat::TERMINAL, "");
        w->write(ds, spec);
    }

    std::unordered_set<WCoefId> coefs;

    for (int i = 7; i <= 10; ++i) {
        coefs.emplace(WCoefMapper::id_of("C" + std::to_string(i)));
    }

    HyperisoMaster hyp = HyperisoMaster();

    HyperisoConfig config;

    
    std::string lha_path = "lha/si_input.flha";

    hyp.init(lha_path, config);

    WilsonInterface wi;
    WilsonBuildConfig wbc;

    wbc.matching_scale = 81;
    wbc.hadronic_scale = 4.18;
    wbc.order = QCDOrder::LO;
    wbc.groups = {GroupMapper::to_id(WGroup::B)};
    wi.build(wbc);

    auto wilsonExtractor =
        std::make_shared<WilsonCoeffExtractor>(wi, wbc, coefs);

    std::shared_ptr<UserParameterProxy> upp;

    if (config.model == Model::SM) {
        upp = std::make_shared<UserParameterProxy>(
            std::vector<ParameterType>{ParameterType::SM, ParameterType::WILSON});
    }
    else {
        upp = std::make_shared<UserParameterProxy>(
            std::vector<ParameterType>{ParameterType::SM,
                                    ParameterType::BSM,
                                    ParameterType::WILSON});
    }

    OutputSpec spec2;
    spec2.csv_write_header = true;

    std::vector<YamlScanParam> scan_params = YamlInputReader("scan.yml").get_scan_params();


    ScanRunner runner(*upp, scan_params, wilsonExtractor);

    DataSet ds2 = runner.run(spec2);

    ds2.meta["model"] = std::string("SM"); 
    ds2.meta["Q_match"] = wbc.matching_scale;
    ds2.meta["Q"] = wbc.hadronic_scale;
    ds2.meta["qcd_order"] = std::string("LO");

    {
        auto w = make_writer(OutputFormat::CSV, "wilson_scan.csv");
        w->write(ds2, spec2);
        std::cout << "Wrote wilson_scan.csv\n";
    }

    {
        auto w = make_writer(OutputFormat::JSON, "wilson_scan.json");
        w->write(ds2, spec2);
        std::cout << "Wrote wilson_scan.json\n";
    }

    // {
    //     auto w = make_writer(OutputFormat::TERMINAL, "");
    //     w->write(ds2, spec2);
    // }
    return 0;
}
