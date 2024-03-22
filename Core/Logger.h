#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <iostream>
#include <fstream>

class Logger {
public:
    enum class LogLevel {
        INFO,
        WARN,
        ERROR,
        DEBUG,
        TRACE,
        VERBOSE
    };

    static Logger* getInstance();

    void setLevel(LogLevel level);
    void setLogFile(const std::string& filename);
    void log(LogLevel messageLevel, const std::string& message);

    void info(const std::string& message);
    void warn(const std::string& message);
    void error(const std::string& message);
    void debug(const std::string& message);
    void trace(const std::string& message);
    void verbose(const std::string& message);
private:
    static Logger* instance;
    Logger() {}
    LogLevel level = LogLevel::INFO;
    std::ofstream logFile;
    std::string toString(LogLevel level);
};

#endif // LOGGER_H
