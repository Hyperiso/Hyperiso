#include <cassert>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <stdexcept>

#include "CompilerStrategy.h"
#include "MakeCompilerStrategy.h"

namespace fs = std::filesystem;

static void write_file(const fs::path& p, const std::string& content) {
    fs::create_directories(p.parent_path());
    std::ofstream f(p);
    f << content;
}

static void make_executable(const fs::path& p) {
    std::error_code ec;
    fs::permissions(p,
        fs::perms::owner_exec | fs::perms::owner_write | fs::perms::owner_read,
        fs::perm_options::add, ec);
}

int main() {
    std::cout << "== CompilerStrategy UNIT ==\n";

    fs::path tmp = fs::temp_directory_path() / ("compiler_unit_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(tmp);

    fs::path ok = tmp / "ok.sh";
    write_file(ok, "#!/bin/sh\nexit 0\n");
    make_executable(ok);
    assert(executeCommand(ok.string()) == true);

    fs::path ko = tmp / "ko.sh";
    write_file(ko, "#!/bin/sh\nexit 2\n");
    make_executable(ko);

    bool threw_on_failure = false;
    try {
        (void)executeCommand(ko.string());
    } catch (const std::runtime_error&) {
        threw_on_failure = true;
    }
    assert(threw_on_failure);

    MakeCompilerStrategy strat("SM", "none");

    fs::path missing = tmp / "missing.bin";
    assert(strat.check_if_compile(missing.string()) == false);

    fs::path empty = tmp / "empty.bin";
    write_file(empty, "");
    assert(strat.check_if_compile(empty.string()) == false);

    fs::path nonempty = tmp / "nonempty.bin";
    write_file(nonempty, "x");
    assert(strat.check_if_compile(nonempty.string()) == true);

    std::cout << "UNIT OK\n";
    return 0;
}
