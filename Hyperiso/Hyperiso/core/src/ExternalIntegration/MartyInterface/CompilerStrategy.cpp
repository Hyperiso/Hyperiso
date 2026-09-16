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

bool executeCommandStreaming(const std::string& command) {
    std::array<char, 512> buffer;
    std::string result;

    FILE* pipe = popen((command + " 2>&1").c_str(), "r");
    if (!pipe) {
        throw std::runtime_error("Error while opening command pipe for: " + command);
    }

    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        const std::string chunk(buffer.data());
        result += chunk;

        // MARTY prints one line for every vertex-set candidate/amplitude.  For
        // large models this can mean thousands of low-value lines and drowns
        // out HyperIso's higher-level analytical progress messages.  Keep the
        // complete output in `result` so failures still report it verbatim,
        // but suppress these repetitive lines from the successful live stream.
        const bool noisy_marty_enumeration =
            chunk.find("possible sets of vertices found") != std::string::npos ||
            chunk.find("new particle amplitude found (set of vertices") != std::string::npos ||
            chunk.find("total particle amplitudes found") != std::string::npos;

        if (!noisy_marty_enumeration) {
            std::cout << chunk << std::flush;
        }
    }

    const int status = pclose(pipe);
    const bool ok = WIFEXITED(status) && WEXITSTATUS(status) == 0;

    if (!ok) {
        std::string message = "Command failed:\n" + command + "\n";
        if (!result.empty()) {
            message += "Output (already streamed above):\n" + result;
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
