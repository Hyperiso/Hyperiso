#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include "ModelFileChecker.h"

namespace fs = std::filesystem;

static void write_file(const fs::path& p, const std::string& content) {
    fs::create_directories(p.parent_path());
    std::ofstream f(p);
    f << content;
    f.flush();
}

int main() {
    std::cout << "== ModelFileChecker UNIT ==\n";

    const fs::path root = fs::temp_directory_path() / "mfc_unit";
    const fs::path templ = root / "THDM_model.hpp";
    const fs::path ntempl= root / "ZPrime_model.hpp";

    write_file(templ,
        "template <int N>\n"
        "class THDM_Model : public mty::Model {};\n"
    );

    write_file(ntempl,
        "class ZPrime_Model : public mty::Model {};\n"
    );

    {
        ModelFileChecker c(templ.string());
        assert(c.isAnyModelTemplate() == true);
    }

    {
        ModelFileChecker c(ntempl.string());
        assert(c.isAnyModelTemplate() == false);
    }

    bool threw = false;
    try {
        ModelFileChecker c((root/"missing.hpp").string());
        (void)c.isAnyModelTemplate();
    } catch (const std::runtime_error&) { threw = true; }
    assert(threw);

    std::cout << "UNIT OK\n";
    return 0;
}
