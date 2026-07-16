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
        fw.add_input_reader(ofs);
        fw.add_argpars(ofs);
        fw.add_output_writer(ofs);
    }
    const std::string s = read_all(out);

    assert(s.find(mgr->getParamFileName()) != std::string::npos);
    assert(s.find(mgr->getCsvWilsonFileName()) != std::string::npos);
    assert(s.find("--Q_match") != std::string::npos);

    FileNameManager::clearTestingRoots();
    std::cout << "UNIT OK\n";
    return 0;
}
