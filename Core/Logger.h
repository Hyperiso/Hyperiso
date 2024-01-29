#ifndef LOGGER_H
#define LOGGER_H

#include <string>

class Logger {
public:
    enum class LogLevel {
        INFO,
        WARN,
        ERROR,
        DEBUG
    };

    static Logger& getInstance();

    void setLevel(LogLevel level);
    void log(LogLevel messageLevel, const std::string& message);

    void info(const std::string& message);
    void warn(const std::string& message);
    void error(const std::string& message);
    void debug(const std::string& message);

private:
    Logger();
    LogLevel level = LogLevel::INFO;
    std::string toString(LogLevel level);
};

#endif // LOGGER_H
