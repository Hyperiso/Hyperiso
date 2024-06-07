#include "Logger.h"

Logger* Logger::instance = nullptr;

Logger* Logger::getInstance() {
    if (!Logger::instance) {
        Logger::instance = new Logger();
    }
    return Logger::instance;
}

void Logger::setLevel(LogLevel level) {
    this->level = level;
}

void Logger::setLogFile(const std::string& filename) {
    logFile.open(filename);
}

std::string Logger::toString(LogLevel level) {
    switch (level) {
        case LogLevel::INFO:  return "INFO";
        case LogLevel::WARN:  return "WARN";
        case LogLevel::ERROR: return "ERROR";
        case LogLevel::DEBUG: return "DEBUG";
        case LogLevel::TRACE: return "TRACE";
        case LogLevel::VERBOSE: return "VERBOSE";
        default: return "UNKNOWN";
    }
}
