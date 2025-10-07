#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>
#include "GeneralModelModifier.h"

namespace fs = std::filesystem;

static void write_file(const fs::path& p, const std::string& content) {
    fs::create_directories(p.parent_path());
    std::ofstream f(p);
    f << content;
    f.flush();
}

int main() {
    std::cout << "== GeneralModelModifier UNIT ==\n";

    const fs::path root = fs::temp_directory_path() / "gmm_unit";
    const fs::path thdm_hdr = root / "models" / "thdm.hpp";
    const fs::path zprime_hdr = root / "models" / "zprime.hpp";

    write_file(thdm_hdr,
        "template <int N>\n"
        "class THDM_Model : public mty::Model {};\n"
    );

    write_file(zprime_hdr,
        "class ZPrime_Model : public mty::Model {};\n"
    );

    {
        GeneralModelModifier mod(/*wilson*/"none", /*model*/"THDM", thdm_hdr.string());
        std::string line = "SM_Model sm;";
        mod.modifyLine(line);
        assert(line.find("THDM_Model<2> sm;") != std::string::npos);
    }

    {
        GeneralModelModifier mod("none", "ZPrime", zprime_hdr.string());
        std::string line = "SM_Model sm;";
        mod.modifyLine(line);
        assert(line.find("ZPrime_Model sm;") != std::string::npos);
    }

    {
        GeneralModelModifier mod("none", "ZPrime", zprime_hdr.string());
        std::string line = "int flag_SM = 0;";
        mod.modifyLine(line);
        // remplace le premier “SM” rencontré → flag_ZPrime
        assert(line.find("flag_ZPrime = 0;") != std::string::npos);
    }

    {
        GeneralModelModifier mod("none", "ZPrime", zprime_hdr.string());
        std::ostringstream oss;
        fs::path out = root / "out.hpp";
        {
            std::ofstream f(out);
            std::string l = "#include <iostream>";
            mod.addLine(f, l, /*addBefore*/false);
            mod.addLine(f, "// nothing", false);
        }
        std::ifstream ifs(out);
        std::stringstream buf; buf << ifs.rdbuf();
        std::string s = buf.str();

        assert(s.find("#include <iostream>") != std::string::npos);
        assert(s.find(zprime_hdr.string()) != std::string::npos);
        assert(s.find("marty.h") != std::string::npos);
        assert(s.find("// nothing") != std::string::npos);
    }

    std::cout << "UNIT OK\n";
    return 0;
}
