#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <sstream>

class Logger {
public:
    enum class LogLevel {
        TRACE,
        VERBOSE,
        DEBUG,
        INFO,
        WARN,
        ERROR
    };

    static Logger* getInstance();

    void setLevel(LogLevel level);
    void setLogFile(const std::string& filename);
    template<typename... Args>
    void log(LogLevel messageLevel, Args... args);

    template<typename... Args>
    void info(Args... args) { log(LogLevel::INFO, args...); }
    template<typename... Args>
    void warn(Args... args) { log(LogLevel::WARN, args...); }
    template<typename... Args>
    void error(Args... args) { log(LogLevel::ERROR, args...); }
    template<typename... Args>
    void debug(Args... args) { log(LogLevel::DEBUG, args...); }
    template<typename... Args>
    void trace(Args... args) { log(LogLevel::TRACE, args...); }
    template<typename... Args>
    void verbose(Args... args) { log(LogLevel::VERBOSE, args...); }

private:
    static Logger* instance;
    Logger() {}
    LogLevel level = LogLevel::INFO;
    std::ofstream logFile;
    std::string toString(LogLevel level);
    
    template<typename T>
    void logMessage(std::ostream& os, T value);
    template<typename T, typename... Args>
    void logMessage(std::ostream& os, T value, Args... args);
};

#include "Logger.tpp"

#endif // LOGGER_H

