#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>
#include "ParamWriter.h"

namespace fs = std::filesystem;

static std::string slurp(const fs::path& p){
    std::ifstream f(p);
    std::stringstream ss; ss << f.rdbuf();
    return ss.str();
}

int main(){
    std::cout << "== ParamWriter UNIT ==\n";
    ParamWriter pw;

    std::unordered_map<std::string,double> ps{
        {"alpha", 1.0}, {"beta", 2.5}, {"gamma", 0.0},
        {"precision_probe", 0.123456789012345}
    };
    const fs::path out = fs::temp_directory_path() / "pw_unit.csv";
    {
        std::ofstream ofs(out);
        pw.writeParams(ofs, ps);
    }

    const auto s = slurp(out);

    assert(s.find("alpha,1")  != std::string::npos);
    assert(s.find("beta,2.5") != std::string::npos);
    assert(s.find("gamma,0")  != std::string::npos);
    assert(s.find("precision_probe,0.123457") != std::string::npos);
    assert(s.find("precision_probe,0.123456789") == std::string::npos);

    std::cout << "✅ UNIT OK\n";
    return 0;
}
