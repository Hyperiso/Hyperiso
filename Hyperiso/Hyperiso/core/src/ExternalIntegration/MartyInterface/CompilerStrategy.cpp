#include "CompilerStrategy.h"

bool executeCommand(const std::string& command) {
    std::array<char, 128> buffer;
    std::string result;
    FILE* pipe = popen((command + " 2>&1").c_str(), "r");
    if (!pipe) {
        std::cerr << "Error while pipe oppening !" << std::endl;
        return false;
    }
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        result += buffer.data();
    }
    int returnCode = pclose(pipe);
    
    if (returnCode != 0) {
        std::cerr << "Error while command execution :\n" << result;
        return false;
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
    std::cout << "Already compiled !" << std::endl;
    return true;
}