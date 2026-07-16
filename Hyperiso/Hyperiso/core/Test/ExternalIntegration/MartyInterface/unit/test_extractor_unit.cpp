#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include "Extractor.h"

namespace fs = std::filesystem;

static void write_file(const fs::path& p, const std::string& s) {
    fs::create_directories(p.parent_path());
    std::ofstream f(p); f << s; f.flush();
}

int main() {
    std::cout << "== Extractor UNIT ==\n";
    const fs::path f = fs::temp_directory_path() / "extractor_unit" / "source.cpp";

    write_file(f,
        " // noise\n"
        "csl::InitSanitizer<real_t> g_L { \"g_L\" };\n"
        "csl::InitSanitizer<real_t> alphaS { \"alphaS\" };\n"
        "csl::InitSanitizer<complex_t> Yb { \"Yb\" };\n"
        " // another line\n"
    );

    auto v = Extractor::extract(f.string());
    assert(v.size() == 3);

    auto has = [&](const std::string& type, const std::string& name, bool cpx){
        for (auto& p : v) if (p.type==type && p.name==name && p.complex==cpx) return true;
        return false;
    };
    assert(has("g_L","g_L",false));
    assert(has("alphaS","alphaS",false));
    assert(has("Yb","Yb",true));

    std::cout << "✅ UNIT OK\n";
    return 0;
}
