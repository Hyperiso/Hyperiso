#include <cassert>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <string>

#include "MakeCompilerStrategy.h"

namespace fs = std::filesystem;

static void write_file(const fs::path& p, const std::string& content) {
    fs::create_directories(p.parent_path());
    std::ofstream f(p);
    f << content;
}

static std::string slurp(const fs::path& p) {
    std::ifstream f(p);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

int main() {
    std::cout << "== CompilerStrategy INTEGRATION ==\n";

    const auto stamp = std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    fs::path root = fs::temp_directory_path() / ("compiler_integ_" + stamp);
    fs::path src  = root / "src";
    fs::path bin  = root / "bin";
    fs::create_directories(src);
    fs::create_directories(bin);

    fs::path out = bin / "runner.sh";
    fs::path runlog = root / "run.log";
    fs::path compiled_marker = root / "compiled.marker";

    std::string mk;
    mk += "all:\n";
    mk += "\t@mkdir -p " + bin.string() + "\n";
    mk += "\t@echo '#!/bin/sh' > " + out.string() + "\n";
    mk += "\t@echo 'printf \"%s \" \"$$@\" > " + runlog.string() + "' >> " + out.string() + "\n";
    mk += "\t@echo 'printf \"\\n\" >> " + runlog.string() + "' >> " + out.string() + "\n";
    mk += "\t@chmod +x " + out.string() + "\n";
    mk += "\t@echo compiled > " + compiled_marker.string() + "\n";
    write_file(src / "Makefile", mk);

    MakeCompilerStrategy strat("SM", "none");

    strat.set_Q_match(80.379);
    strat.compile_run(src.string(), out.string());

    assert(fs::exists(out) && fs::file_size(out) > 0);
    assert(fs::exists(compiled_marker));

    std::string log1 = slurp(runlog);

    assert(log1.find("-Q") != std::string::npos);
    assert(log1.find("80.379") != std::string::npos);

    fs::remove(compiled_marker);
    strat.compile_run(src.string(), out.string());

    assert(!fs::exists(compiled_marker));

    std::string log2 = slurp(runlog);
    assert(log2.find("-Q") != std::string::npos);

    strat.set_Q_match(91.1876);
    strat.compile_run(src.string(), out.string());
    std::string log3 = slurp(runlog);
    assert(log3.find("91.1876") != std::string::npos || log3.find("91.187") != std::string::npos);

    std::cout << "INTEGRATION OK\n";
    return 0;
}
