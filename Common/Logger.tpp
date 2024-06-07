#include "Logger.h"

// Logger* Logger::instance = nullptr;

// Logger* Logger::getInstance() {
//     if (!Logger::instance) {
//         Logger::instance = new Logger();
//     }
//     return Logger::instance;
// }

// void Logger::setLevel(LogLevel level) {
//     this->level = level;
// }

// void Logger::setLogFile(const std::string& filename) {
//     logFile.open(filename);
// }

// std::string Logger::toString(LogLevel level) {
//     switch (level) {
//         case LogLevel::INFO:  return "INFO";
//         case LogLevel::WARN:  return "WARN";
//         case LogLevel::ERROR: return "ERROR";
//         case LogLevel::DEBUG: return "DEBUG";
//         case LogLevel::TRACE: return "TRACE";
//         case LogLevel::VERBOSE: return "VERBOSE";
//         default: return "UNKNOWN";
//     }
// }

template<typename T>
void Logger::logMessage(std::ostream& os, T value) {
    os << value;
}

template<typename T, typename... Args>
void Logger::logMessage(std::ostream& os, T value, Args... args) {
    os << value << ", ";
    logMessage(os, args...);
}

template<typename... Args>
void Logger::log(LogLevel messageLevel, Args... args) {
    if (messageLevel >= level) {
        std::ostringstream messageStream;
        logMessage(messageStream, args...);
        std::string message = messageStream.str();

        if (logFile.is_open()) {
            logFile << "[" << toString(messageLevel) << "] " << message << std::endl;
        }
        std::cout << "[" << toString(messageLevel) << "] " << message << std::endl;
    }
}
