#include "CompilerStrategy.h"

#include <stdexcept>
#include <sys/wait.h>

bool executeCommand(const std::string& command) {
    std::array<char, 128> buffer;
    std::string result;

    FILE* pipe = popen((command + " 2>&1").c_str(), "r");
    if (!pipe) {
        throw std::runtime_error("Error while opening command pipe for: " + command);
    }

    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        result += buffer.data();
    }

    const int status = pclose(pipe);
    const bool ok = WIFEXITED(status) && WEXITSTATUS(status) == 0;

    if (!ok) {
        std::string message = "Command failed:\n" + command + "\n";
        if (!result.empty()) {
            message += "Output:\n" + result;
        }
        throw std::runtime_error(message);
    }

    return true;
}

bool CompilerStrategy::check_if_compile(const std::string& outputBinary) {
    struct stat buffer;
    if (stat(outputBinary.c_str(), &buffer) != 0) {
        return false;
    }
    if (buffer.st_size == 0) {
        return false;
    }
    LOG_DEBUG("Already compiled !");
    return true;
}
