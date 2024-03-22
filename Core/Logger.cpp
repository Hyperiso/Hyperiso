#include "Logger.h"
#include <iostream>

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

void Logger::log(LogLevel messageLevel, const std::string& message) {
    if (messageLevel >= level) {
        // Écriture dans le fichier de log si le fichier est ouvert
        if (logFile.is_open()) {
            logFile << "[" << toString(messageLevel) << "] " << message << std::endl;
        }
        // Écriture dans la sortie standard
        std::cout << "[" << toString(messageLevel) << "] " << message << std::endl;
    }
}

void Logger::info(const std::string& message) {
    log(LogLevel::INFO, message);
}

void Logger::warn(const std::string& message) {
    log(LogLevel::WARN, message);
}

void Logger::error(const std::string& message) {
    log(LogLevel::ERROR, message);
}

void Logger::debug(const std::string& message) {
    log(LogLevel::DEBUG, message);
}

void Logger::trace(const std::string& message) {
    log(LogLevel::TRACE, message);
}

void Logger::verbose(const std::string& message) {
    log(LogLevel::VERBOSE, message);
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

Logger* Logger::instance = nullptr;