#include "GppCompilerStrategy.h"
#include "MartyRuntimeConfig.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cctype>
#include <filesystem>
#include <sstream>
#include <stdexcept>

namespace {

std::string trim(std::string value) {
    const auto not_space = [](unsigned char c) { return !std::isspace(c); };
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
    value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(), value.end());
    return value;
}

std::string command_output(const std::string& command) {
    std::array<char, 256> buffer{};
    std::string result;

    FILE* pipe = popen((command + " 2>/dev/null").c_str(), "r");
    if (!pipe) {
        return {};
    }

    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        result += buffer.data();
    }

    const int status = pclose(pipe);
    if (status != 0) {
        return {};
    }

    return trim(result);
}

std::filesystem::path find_libgfortran() {
    const std::string from_gpp = command_output("g++ -print-file-name=libgfortran.so");
    if (!from_gpp.empty() && from_gpp != "libgfortran.so" && std::filesystem::exists(from_gpp)) {
        return from_gpp;
    }

    const std::string from_gfortran = command_output("gfortran -print-file-name=libgfortran.so");
    if (!from_gfortran.empty() && from_gfortran != "libgfortran.so" && std::filesystem::exists(from_gfortran)) {
        return from_gfortran;
    }

    return {};
}

std::filesystem::path require_libgfortran() {
    const auto lib = find_libgfortran();
    if (!lib.empty()) {
        return lib;
    }

    throw std::runtime_error(
        "MARTY generated-code compilation requires libgfortran, but the linker "
        "cannot find libgfortran.so.\n"
        "Install the Fortran development package, for example on Ubuntu/Debian:\n"
        "  sudo apt update && sudo apt install gfortran\n"
        "Then verify with:\n"
        "  g++ -print-file-name=libgfortran.so\n"
        "The command must print an absolute path, not just 'libgfortran.so'."
    );
}

} // namespace

void GppCompilerStrategy::compile_run(const std::string& sourceFile, const std::string& outputBinary) {
    if (!this->check_if_compile(outputBinary)) {
        this->compile(sourceFile, outputBinary);
    }

    std::string command_run = "cd " + FileNameManager::getInstance(wilson, model)->getOutputDir() + " && " + outputBinary;
    executeCommand(command_run);
}

void GppCompilerStrategy::compile(const std::string& sourceFile, const std::string& outputBinary) {
    const auto marty = MartyRuntimeConfig::require_available("GppCompilerStrategy::compile");
    if (!marty.valid) {
        return;
    }

    const auto libgfortran = require_libgfortran();

    const std::string include_dir = MartyRuntimeConfig::shell_quote(marty.include_dir);
    const std::string lib_dir = MartyRuntimeConfig::shell_quote(marty.lib_dir);
    const std::string libgfortran_path = MartyRuntimeConfig::shell_quote(libgfortran);

    std::string command_compile = "g++ -o " + outputBinary + " " + sourceFile
        + " -I" + include_dir
        + " -L" + lib_dir
        + " -Wl,-rpath," + lib_dir
        + " -lmarty " + libgfortran_path;

    executeCommand(command_compile);
}
