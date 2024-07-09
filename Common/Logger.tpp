#include "Logger.h"


template<typename T>
void Logger::logMessage(std::ostream& os, T value) {
    os << value;
}

template<typename T, typename... Args>
void Logger::logMessage(std::ostream& os, T value, Args... args) {
    os << value << " ";
    logMessage(os, args...);
}

template<typename... Args>
void Logger::logToFile(std::ostream& os, LogLevel messageLevel, const char* file, int line, const char* func, Args... args) {
    os << "[" << currentDateTime() << "] [" << toString(messageLevel) << "] [" 
       << file << ":" << line << " " << func << "] [Thread " << std::this_thread::get_id() << "] ";
    logMessage(os, args...);
}

template<typename... Args>
void Logger::logToTerminal(std::ostream& os, LogLevel messageLevel, Args... args) {
    os << "[" << currentDateTime() << "] [" << toString(messageLevel) << "] ";
    logMessage(os, args...);
}

template<typename... Args>
void Logger::log(LogLevel messageLevel, const char* file, int line, const char* func, Args... args) {
    if (!enabled || messageLevel < level) {
        return;
    }
    if (messageLevel >= level) {
        std::ostringstream fileMessageStream;
        std::ostringstream terminalMessageStream;

        logToFile(fileMessageStream, messageLevel, file, line, func, args...);
        logToTerminal(terminalMessageStream, messageLevel, args...);

        std::string fileMessage = fileMessageStream.str();
        std::string terminalMessage = terminalMessageStream.str();

        {
            std::lock_guard<std::mutex> lock(queueMutex);
            logQueue.push(fileMessage);
        }
        conditionVar.notify_one();

        std::cout << terminalMessage << std::endl;
    }
}
