#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <sstream>
#include <ctime>
#include <filesystem>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <unordered_map>

#define LOG_INFO(...) Logger::getInstance()->log(Logger::LogLevel::INFO, __FILE__, __LINE__, __func__, __VA_ARGS__)
#define LOG_WARN(...) Logger::getInstance()->log(Logger::LogLevel::WARN, __FILE__, __LINE__, __func__, __VA_ARGS__)
#define LOG_ERROR(...) Logger::getInstance()->log(Logger::LogLevel::ERROR, __FILE__, __LINE__, __func__, __VA_ARGS__)
#define LOG_DEBUG(...) Logger::getInstance()->log(Logger::LogLevel::DEBUG, __FILE__, __LINE__, __func__, __VA_ARGS__)
#define LOG_TRACE(...) Logger::getInstance()->log(Logger::LogLevel::TRACE, __FILE__, __LINE__, __func__, __VA_ARGS__)
#define LOG_VERBOSE(...) Logger::getInstance()->log(Logger::LogLevel::VERBOSE, __FILE__, __LINE__, __func__, __VA_ARGS__)


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
    void setLogDirectory(const std::string& directory, std::size_t maxSize = 1048576);
    
    template<typename... Args>
    void log(LogLevel messageLevel, const char* file, int line, const char* func, Args... args);

    ~Logger();

private:
    static Logger* instance;
    Logger();
    LogLevel level = LogLevel::INFO;
    std::string logDirectory{"logs"};
    std::size_t maxFileSize;
    std::unordered_map<std::thread::id, std::ofstream> logFiles;
    std::mutex logFileMutex;

    void rotateLogFile(std::thread::id threadId);

    std::string toString(LogLevel level);
    std::string currentDateTime();
    std::ofstream& getLogFile();

    template<typename T>
    void logMessage(std::ostream& os, T value);
    template<typename T, typename... Args>
    void logMessage(std::ostream& os, T value, Args... args);

    template<typename... Args>
    void logToFile(std::ostream& os, LogLevel messageLevel, const char* file, int line, const char* func, Args... args);

    template<typename... Args> 
    void logToTerminal(std::ostream& os, LogLevel messageLevel, Args... args); 

    std::queue<std::string> logQueue;
    std::thread loggingThread;
    std::mutex queueMutex;
    std::condition_variable conditionVar;
    std::atomic<bool> exitFlag;

    void processQueue();
};

#include "Logger.tpp"

#endif // LOGGER_H

