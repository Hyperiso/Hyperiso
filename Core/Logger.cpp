#include <string>
#include <iostream>

class Logger {
public:
    enum class LogLevel {
        INFO,
        WARN,
        ERROR,
        DEBUG
    };

    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }

    void setLevel(LogLevel level) {
        this->level = level;
    }

    void log(LogLevel messageLevel, const std::string& message) {
        if (messageLevel >= level) {
            // Implémentation de la logique de log, exemple:
            std::cout << "[" << toString(messageLevel) << "] " << message << std::endl;
        }
    }

    // Méthodes spécifiques pour chaque niveau
    void info(const std::string& message) {
        log(LogLevel::INFO, message);
    }

    void warn(const std::string& message) {
        log(LogLevel::WARN, message);
    }

    void error(const std::string& message) {
        log(LogLevel::ERROR, message);
    }

    void debug(const std::string& message) {
        log(LogLevel::DEBUG, message);
    }

private:
    Logger() {}
    LogLevel level = LogLevel::INFO;

    std::string toString(LogLevel level) {
        switch (level) {
            case LogLevel::INFO:  return "INFO";
            case LogLevel::WARN:  return "WARN";
            case LogLevel::ERROR: return "ERROR";
            case LogLevel::DEBUG: return "DEBUG";
            default: return "UNKNOWN";
        }
    }
};
