#include "ObservablesUtils.h"

#include <fstream>
#include <iomanip>
#include <filesystem>
#include <unordered_map>
#include <sstream>

#include "Logger.h"
#include "WilsonInterface.h"  
#include "CompareCsv.h" 
#include "Configs.h"
#include "ObsParameterProxy.h"

namespace fs = std::filesystem;

static std::string bin_label(const ObservableValue& ov) {
    if (!ov.bin.has_value()) return "";
    std::ostringstream ss;
    ss << " [" << ov.bin.value().first << ", " << ov.bin.value().second << "]";
    return ss.str();
}

void writeDecayObservablesCsv(ObservableInterface& oi,
                              Decays dec,
                              QCDOrder order,
                              const std::string& out_dir)
{
    fs::create_directories(out_dir);

    oi.add_observables(dec, order, /*add_dependencies=*/false);

    std::vector<std::string> headers;
    std::vector<double> values; values.reserve(64);

    const double Q = ObsParameterProxy(ParameterType::WILSON).get_parameter({ParameterType::WILSON, "B_SCALE", 1})->get_val().real();
    const double alpha_s = QCDHelper::alpha_s(Q);

    headers.push_back("Q");
    headers.push_back("alpha_s");
    values.push_back(Q);
    values.push_back(alpha_s);

    for (auto& o : DecayMapper::get_observables(dec)) {
        const auto obs_vals = oi.compute_observable(o);
        if (obs_vals.size() == 1) {
            std::string h = ObservableMapper::str(o);
            headers.push_back(h);
            values.push_back(obs_vals[0].value);
        } else {
            for (const auto& ov : obs_vals) {
                std::string h = ObservableMapper::str(o) + bin_label(ov);
                headers.push_back(h);
                values.push_back(ov.value);
            }
        }
    }

    const std::string file_name = DecayMapper::str(dec) + ".csv";
    fs::path out_path = fs::path(out_dir) / file_name;
    std::ofstream file(out_path.string());
    if (!file.is_open()) {
        LOG_ERROR("FileError", "Could not open for write: ", out_path.string());
    }

    for (size_t i = 0; i < headers.size(); ++i) {
        if (i) file << ',';
        file << headers[i];
    }
    file << '\n';

    file << std::scientific << std::setprecision(4);
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) file << ',';
        file << values[i];
    }
    file << '\n';
    file.close();

    LOG_INFO("Wrote observables CSV for decay", DecayMapper::str(dec), "->", out_path.string());
}

void writeAllDecaysObservablesCsv(ObservableInterface& oi,
                                  QCDOrder order,
                                  const std::string& out_dir)
{
    for (auto& dec : DecayMapper::get_enum()) {
        if (dec == Decays::M0_Mix) continue;
        LOG_INFO("Exporting observables for decay", DecayMapper::str(dec));
        writeDecayObservablesCsv(oi, dec, order, out_dir);
    }
}

void runObservablesTest(const std::string& lha_path,
                        Model model,
                        QCDOrder order,
                        const std::string& out_dir,
                        const std::string& ref_dir,
                        double tolerance)
{
    HyperisoMaster hyp;
    HyperisoConfig cfg;
    cfg.model = model;
    hyp.init(lha_path, cfg);

    ObservableInterface oi;
    writeAllDecaysObservablesCsv(oi, order, out_dir);

    for (auto& dec : DecayMapper::get_enum()) {
        if (dec == Decays::M0_Mix) continue;

        const std::string f = DecayMapper::str(dec) + ".csv";
        fs::path gen = fs::path(out_dir) / f;
        fs::path ref = fs::path(ref_dir) / f;

        LOG_INFO("Comparing", gen.string(), "vs", ref.string());
        if (!compareCSV(gen.string(), ref.string(), tolerance)) {
            std::cerr << "Test failed for " << DecayMapper::str(dec) << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Test passed for " << DecayMapper::str(dec) << std::endl;
    }
}
