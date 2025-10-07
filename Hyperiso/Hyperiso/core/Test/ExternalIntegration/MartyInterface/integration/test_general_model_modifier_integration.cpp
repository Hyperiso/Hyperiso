#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "GeneralModelModifier.h"

namespace fs = std::filesystem;

static void write_file(const fs::path& p, const std::string& content) {
    fs::create_directories(p.parent_path());
    std::ofstream f(p);
    f << content;
    f.flush();
}

static std::string slurp(const fs::path& p) {
    std::ifstream f(p);
    std::stringstream ss; ss << f.rdbuf();
    return ss.str();
}

int main() {
    std::cout << "== GeneralModelModifier INTEGRATION ==\n";

    const fs::path root = fs::temp_directory_path() / "gmm_integ";
    const fs::path zprime_hdr = root / "models" / "zprime.hpp";
    const fs::path thdm_hdr   = root / "models" / "thdm.hpp";

    write_file(zprime_hdr, "class ZPrime_Model : public mty::Model {};\n");
    write_file(thdm_hdr,   "template <int N>\nclass THDM_Model : public mty::Model {};\n");

    const std::vector<std::string> input = {
        "#include <iostream>",
        "int main(){",
        "  SM_Model sm;",
        "  int flag_SM = 0;",
        "  return 0;",
        "}"
    };

    {
        fs::path out = root / "gen_zprime.cpp";
        GeneralModelModifier mod("none", "ZPrime", zprime_hdr.string());
        {
            std::ofstream ofs(out);
            for (auto line : input) {
                mod.modifyLine(line);
                mod.addLine(ofs, line, false);
            }
        }
        const std::string s = slurp(out);

        assert(s.find(zprime_hdr.string()) != std::string::npos);
        assert(s.find("marty.h") != std::string::npos);

        assert(s.find("ZPrime_Model sm;") != std::string::npos);
        assert(s.find("flag_ZPrime = 0;") != std::string::npos);
    }

    {
        fs::path out = root / "gen_thdm.cpp";
        GeneralModelModifier mod("none", "THDM", thdm_hdr.string());
        {
            std::ofstream ofs(out);
            for (auto line : input) {
                mod.modifyLine(line);
                mod.addLine(ofs, line, false);
            }
        }
        const std::string s = slurp(out);
        assert(s.find(thdm_hdr.string()) != std::string::npos);
        assert(s.find("THDM_Model<2> sm;") != std::string::npos);
        assert(s.find("flag_THDM = 0;") != std::string::npos);
    }

    std::cout << "INTEGRATION OK\n";
    return 0;
}
