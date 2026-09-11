#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>
#include "MartyFileWriter.h"
#include "FileNameManager.h"

namespace fs = std::filesystem;

static std::string read_all(const fs::path& p){
    std::ifstream f(p); std::stringstream ss; ss << f.rdbuf(); return ss.str();
}

int main() {
    std::cout << "== MartyFileWriter UNIT ==\n";

    const fs::path root = fs::temp_directory_path() / "fnm_fw_unit";
    const std::string templ = (root / "templ").string() + "/";
    const std::string base  = (root / "base").string() + "/";
    const std::string assets= (root / "assets").string() + "/";
    fs::create_directories(root / "templ");

    FileNameManager::setTestingRoots(templ, base, assets);

    auto mgr = FileNameManager::getInstance("C7","SM");

    MartyFileWriter fw("C7","SM");
    const fs::path out = root / "snippet.cpp";
    {
        std::ofstream ofs(out);
        fw.add_argpars(ofs);
        fw.add_input_reader(ofs);
        fw.add_output_writer(ofs);
    }
    const std::string s = read_all(out);

    assert(s.find(mgr->getParamFileName()) != std::string::npos);
    assert(s.find(mgr->getCsvWilsonFileName()) != std::string::npos);
    assert(s.find("--Q_match") != std::string::npos);
    assert(s.find("--param-file") != std::string::npos);
    assert(s.find("--output-file") != std::string::npos);
    assert(s.find("std::ifstream ParamFile(param_file_path)") != std::string::npos);
    assert(s.find("const std::string& path = output_file_path") != std::string::npos);
    assert(s.find("std::string param_file_path") < s.find("std::ifstream ParamFile(param_file_path)"));

    // C9/CP9 keep the physical non-photon branch at the tiny regulator and
    // evaluate only the raw photon diagnostic with reg_prop=1.  The raw linker
    // must never be added directly to the physical coefficient.
    const fs::path c9_out = root / "c9_photon_veto_snippet.cpp";
    {
        MartyFileWriter c9_fw("C9", "THDM", true, false);
        std::ofstream ofs(c9_out);
        c9_fw.add_output_writer(ofs);
    }
    const std::string c9 = read_all(c9_out);
    const auto small_reg = c9.find("*hyperiso_regprop_it->second = 1e-6");
    const auto non_photon = c9.find("auto hyperiso_bsm_non_photon = C9(param)");
    const auto large_reg = c9.find("*hyperiso_regprop_it->second = 1.0");
    const auto photon_raw_zero = c9.find("auto hyperiso_bsm_photon_raw = 0.0 * hyperiso_bsm_non_photon");
    const auto photon_opt_in = c9.find("if (raw_photon_diagnostic)");
    const auto photon_raw = c9.find("hyperiso_bsm_photon_raw = C9_A(param)");
    const auto physical = c9.find("auto hyperiso_bsm_physical = hyperiso_bsm_non_photon");
    const auto raw_sum = c9.find("auto hyperiso_raw_photon_sum = hyperiso_bsm_non_photon + hyperiso_bsm_photon_raw");
    const auto physical_write = c9.find("writeWilsonCoefficients(\"C9\", hyperiso_bsm_physical");
    assert(small_reg != std::string::npos);
    assert(large_reg != std::string::npos);
    assert(non_photon != std::string::npos);
    assert(photon_raw_zero != std::string::npos);
    assert(photon_opt_in != std::string::npos);
    assert(photon_raw != std::string::npos);
    assert(physical != std::string::npos);
    assert(raw_sum != std::string::npos);
    assert(physical_write != std::string::npos);
    assert(small_reg < non_photon);
    assert(non_photon < photon_raw_zero);
    assert(photon_raw_zero < photon_opt_in);
    assert(photon_opt_in < large_reg);
    assert(large_reg < photon_raw);
    assert(photon_raw < physical);
    assert(physical < raw_sum);
    assert(c9.find("hyperiso_bsm_physical = hyperiso_bsm_non_photon + hyperiso_bsm_photon") == std::string::npos);
    assert(c9.find("writeWilsonCoefficients(\"C9\", hyperiso_raw_photon_sum") == std::string::npos);

    FileNameManager::clearTestingRoots();
    std::cout << "UNIT OK\n";
    return 0;
}
