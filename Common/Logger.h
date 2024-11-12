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
#include "LoggingException.h"

/**
 * @def LOG_INFO(...)
 * @brief Macro for logging informational messages.
 *
 * @def LOG_WARN(...)
 * @brief Macro for logging warning messages.
 *
 * @def LOG_ERROR(type, ...)
 * @brief Macro for logging error messages and terminating the application.
 *
 * @def LOG_DEBUG(...)
 * @brief Macro for logging debug messages.
 *
 * @def LOG_TRACE(...)
 * @brief Macro for logging trace messages.
 *
 * @def LOG_VERBOSE(...)
 * @brief Macro for logging verbose messages.
 */
#define LOG_INFO(...) Logger::getInstance()->log(Logger::LogLevel::INFO, __FILE__, __LINE__, __func__, "", __VA_ARGS__)
#define LOG_WARN(...) Logger::getInstance()->log(Logger::LogLevel::WARN, __FILE__, __LINE__, __func__, "",__VA_ARGS__)
#define LOG_ERROR(type, ...) do { \
    Logger::getInstance()->log(Logger::LogLevel::ERROR, __FILE__, __LINE__, __func__, type, __VA_ARGS__); \
    std::terminate(); \
} while(0)
#define LOG_DEBUG(...) Logger::getInstance()->log(Logger::LogLevel::DEBUG, __FILE__, __LINE__, __func__,"", __VA_ARGS__)
#define LOG_TRACE(...) Logger::getInstance()->log(Logger::LogLevel::TRACE, __FILE__, __LINE__, __func__,"", __VA_ARGS__)
#define LOG_VERBOSE(...) Logger::getInstance()->log(Logger::LogLevel::VERBOSE, __FILE__, __LINE__, __func__,"", __VA_ARGS__)

/**
 * @class Logger
 * @brief A singleton logging utility supporting multiple logging levels and concurrent logging to file and terminal.
 */
class Logger {
public:
    /**
     * @enum LogLevel
     * @brief Logging levels supported by the Logger.
     */
    enum class LogLevel {
        TRACE,      /**< Trace level for detailed memory information. */
        VERBOSE,    /**< Verbose level for general debug messages. */
        DEBUG,      /**< Debug level for debugging messages. */
        INFO,       /**< Info level for informational messages. */
        WARN,       /**< Warn level for warning messages. */
        ERROR       /**< Error level for critical error messages. */
    };

    /**
     * @brief Retrieves the singleton instance of the Logger.
     * @return Pointer to the singleton Logger instance.
     */
    static Logger* getInstance();

    /**
     * @brief Sets the logging level.
     * @param level The minimum logging level to be logged.
     */
    void setLevel(LogLevel level);

    /**
     * @brief Sets the directory for storing log files.
     * @param directory Directory path where log files are stored.
     * @param maxSize Maximum allowed size for each log file in bytes, default is 10048576.
     */
    void setLogDirectory(const std::string& directory, std::size_t maxSize = 10048576);

    /**
     * @brief Enables or disables the Logger.
     * @param enabled Boolean flag to enable or disable logging.
     */
    void setEnabled(bool enabled);

    /**
     * @brief Retrieves the log filename for the specified thread.
     * @param threadId The thread ID.
     * @return Log filename associated with the specified thread.
     */
    std::string getLogFilenameForThread(std::thread::id threadId);

    /**
     * @brief Retrieves a log filename with a suffix.
     * @param baseName Base name for the log file.
     * @param suffix Integer suffix for the filename.
     * @return Full filename with suffix.
     */
    std::string getLogFilenameWithSuffix(const std::string& baseName, int suffix);
    
    std::thread::id threadId; /**< The thread ID for logging purposes. */
    
    /**
     * @brief Logs a message with the specified level and details.
     * @param messageLevel The logging level of the message.
     * @param file Source file where the log is made.
     * @param line Line number in the source file.
     * @param func Function name.
     * @param errorType Optional error type string for error-level logs.
     * @param args The message to log.
     * @throws std::runtime_error If logging to file fails.
     */
    template<typename... Args>
    void log(LogLevel messageLevel, const char* file, int line, const char* func, const char* errorType, Args... args) noexcept(false);
    
    /**
     * @brief Destructor for Logger, ensures graceful shutdown of logging.
     */
    ~Logger();

private:
    static Logger* instance; /**< Singleton instance of Logger. */

    /**
     * @brief Private constructor for singleton pattern.
     */
    Logger();

    LogLevel level = LogLevel::INFO; /**< Current logging level. */
    std::string logDirectory{"logs"}; /**< Directory path for storing log files. */
    std::size_t maxFileSize{10048576}; /**< Maximum size for each log file. */
    std::unordered_map<std::thread::id, std::ofstream> logFiles; /**< Map of log files by thread ID. */
    std::mutex logFileMutex; /**< Mutex to protect log file access. */
    bool enabled{true}; /**< Flag to enable or disable logging. */

    /**
     * @brief Rotates the log file if the maximum file size is exceeded.
     * @param threadId The thread ID associated with the log file.
     */
    void rotateLogFile(std::thread::id threadId);

    /**
     * @brief Converts a LogLevel to a string representation.
     * @param level The LogLevel to convert.
     * @return String representation of the logging level.
     */
    std::string toString(LogLevel level);

    /**
     * @brief Gets the current date and time as a string.
     * @return Current date and time in formatted string.
     */
    std::string currentDateTime();

    /**
     * @brief Retrieves or opens a log file for a specific thread.
     * @param threadId The thread ID.
     * @return Output stream for the log file.
     * @throws std::runtime_error If unable to open log file.
     */
    std::ofstream& getLogFile(std::thread::id threadId);

    

    template<typename T>
    void logMessage(std::ostream& os, T value);
    template<typename T, typename... Args>
    void logMessage(std::ostream& os, T value, Args... args);

    template<typename... Args>
    void logToFile(std::ostream& os, LogLevel messageLevel, const char* file, int line, const char* func, const char* errorType, Args... args);

    template<typename... Args> 
    void logToTerminal(std::ostream& os, LogLevel messageLevel, const char* errorType, Args... args); 

    std::queue<std::string> logQueue; /**< Queue for holding log messages before writing to file. */
    std::thread loggingThread; /**< Thread for processing log queue asynchronously. */
    std::mutex queueMutex; /**< Mutex to protect access to log queue. */
    std::condition_variable conditionVar; /**< Condition variable to synchronize queue processing. */
    std::atomic<bool> exitFlag; /**< Flag to signal termination of logging thread. */

    /**
     * @brief Processes the log queue, writing each entry to file.
     */
    void processQueue();

};

#include "Logger.tpp"

#endif // LOGGER_H

