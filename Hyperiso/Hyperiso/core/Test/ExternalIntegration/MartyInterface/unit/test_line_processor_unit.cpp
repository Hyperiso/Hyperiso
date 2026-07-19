#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

#include "LineProcessor.h"
#include "ModelWriter.h"
#include "MartyFileWriter.h"
#include "IncludeManager.hpp"
#include "FileNameManager.h"

namespace fs = std::filesystem;

static std::string slurp(const fs::path& p){
    std::ifstream f(p);
    std::stringstream ss; ss << f.rdbuf();
    return ss.str();
}

int main(){
    std::cout << "== LineProcessor/ModelWriter UNIT ==\n";

    const fs::path root = fs::temp_directory_path() / "lp_mw_unit";
    const std::string templ = (root / "templ").string() + "/";
    const std::string base  = (root / "base").string() + "/";
    const std::string assets= (root / "assets").string() + "/";
    std::filesystem::create_directories(root / "templ");

    FileNameManager::setTestingRoots(templ, base, assets);
    (void)FileNameManager::getInstance("W1","SM");

    const fs::path in  = root / "in.cpp";
    const fs::path out = root / "out.cpp";
    {
        std::ofstream ofs(in);
        ofs << "#include <iostream>\n";
        ofs << "using namespace std;\n";
        ofs << "int main(){\n";
        ofs << "  // body\n";
        ofs << "  return 0;\n";
        ofs << "}\n";
    }

    IncludeManager im;
    MartyFileWriter fw("W1","SM");
    LineProcessor lp(im, fw, /*forceMode*/false);
    ParamWriter pw;
    ModelWriter mw(lp, pw);

    {
        std::ifstream ifs(in);
        std::ofstream ofs(out);
        mw.writeModel(ifs, ofs);
    }

    const auto s = slurp(out);

    assert(s.find("#include <fstream>") != std::string::npos);
    assert(s.find("#include \"csv_helper.h\"") != std::string::npos);
    assert(s.find("int main(int argc, char** argv)") != std::string::npos);
    assert(s.find("param_t param;") != std::string::npos);
    assert(s.find("--Q_match") != std::string::npos);
    assert(s.find("--param-file") != std::string::npos);
    assert(s.find("--output-file") != std::string::npos);
    assert(s.find("std::string param_file_path") < s.find("std::ifstream ParamFile(param_file_path)"));

    const fs::path csv = root / "params.csv";
    {
        std::ofstream ofs(csv);
        std::unordered_map<std::string,double> ps{{"p1",1.2},{"p2",3.4}};
        mw.writeParam(ofs, ps);
    }
    const auto sc = slurp(csv);
    assert(sc.find("p1,1.2") != std::string::npos);
    assert(sc.find("p2,3.4") != std::string::npos);

    FileNameManager::clearTestingRoots();
    std::cout << "UNIT OK\n";
    return 0;
}
