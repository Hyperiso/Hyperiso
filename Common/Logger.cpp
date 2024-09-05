#include "Logger.h"
#include <filesystem>
#include <iostream>

Logger* Logger::instance = nullptr;

Logger* Logger::getInstance() {
    if (!Logger::instance) {
        Logger::instance = new Logger();
    }
    return Logger::instance;
}

Logger::Logger() : exitFlag(false) {
    if (!std::filesystem::exists(this->logDirectory)) {
        if(std::filesystem::create_directories(this->logDirectory)) {

        }
        else {
            std::cerr << "Error during logDirectory creation. " << std::endl;
        }
    }
    loggingThread = std::thread(&Logger::processQueue, this);
    this->threadId = loggingThread.get_id();
}

Logger::~Logger() {
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        exitFlag = true;
    }
    conditionVar.notify_all();
    loggingThread.join();

    for (auto& logFilePair : logFiles) {
        logFilePair.second.close();
    }
}

void Logger::setLevel(LogLevel level) {
    this->level = level;
}

void Logger::setLogDirectory(const std::string& directory, std::size_t maxSize) {
    logDirectory = directory;
    maxFileSize = maxSize;
    if (!std::filesystem::exists(logDirectory)) {
        std::filesystem::create_directories(logDirectory);
    }
}

void Logger::setEnabled(bool enabled) {
    this->enabled = enabled;
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


std::string Logger::currentDateTime() {
    std::time_t now = std::time(nullptr);
    char buf[80];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d_%H", std::localtime(&now));
    return std::string(buf);
}

std::string Logger::getLogFilenameWithSuffix(const std::string& baseName, int suffix) {
    std::ostringstream filename;
    filename << baseName << "_" << suffix << ".log";
    return filename.str();
}

std::string Logger::getLogFilenameForThread(std::thread::id threadId) {
    std::ostringstream filename;
    filename << logDirectory << "/log_" << currentDateTime() << "_thread_" << threadId;
    return filename.str();
}

std::ofstream& Logger::getLogFile(std::thread::id threadId) {
    std::lock_guard<std::mutex> lock(logFileMutex);

    auto it = logFiles.find(threadId);
    if (it == logFiles.end()) {
        std::string filepath = getLogFilenameForThread(threadId) + ".log";
        logFiles[threadId].open(filepath, std::ios_base::app);
        if (!logFiles[threadId].is_open()) {
            throw std::runtime_error("Cannot open log file for thread");
        }
    }
    return logFiles[threadId];
}

void Logger::rotateLogFile(std::thread::id threadId) {
    std::lock_guard<std::mutex> lock(logFileMutex);

    auto& logFile = logFiles[threadId];
    logFile.close();

    std::string baseFilepath = getLogFilenameForThread(threadId);
    std::string oldFilepath = baseFilepath + ".log";

    int suffix = 1;
    std::string newFilename = getLogFilenameWithSuffix(baseFilepath, suffix);
    while (std::filesystem::exists(newFilename)) {
        ++suffix;
        newFilename = getLogFilenameWithSuffix(baseFilepath, suffix);
    }
    std::filesystem::rename(oldFilepath, newFilename);

    logFiles[threadId].open(oldFilepath, std::ios_base::app);
    if (!logFiles[threadId].is_open()) {
        throw std::runtime_error("Cannot open new log file for thread");
    }
}

void Logger::processQueue() {
    while (!exitFlag) {
        std::unique_lock<std::mutex> lock(queueMutex);
        conditionVar.wait(lock, [this] { return !logQueue.empty() || exitFlag; });

        while (!logQueue.empty()) {
            std::string logEntry = logQueue.front();
            logQueue.pop();
            lock.unlock();

            // this->threadId = std::this_thread::get_id();
            std::string logFilePath = getLogFilenameForThread(this->threadId) + ".log";

            if (std::filesystem::exists(logFilePath) && std::filesystem::file_size(logFilePath) >= maxFileSize) {
                rotateLogFile(threadId);
                std::cout << logFilePath << " " <<  std::filesystem::file_size(logFilePath) << std::endl;
                logFilePath = getLogFilenameForThread(threadId) + ".log"; // Update the file path after rotation
            }

            std::ofstream& logFile = getLogFile(this->threadId);
            logFile << logEntry << std::endl;

            lock.lock();
        }
    }

    // Drain the queue when exiting
    while (!logQueue.empty()) {
        std::string logEntry = logQueue.front();
        logQueue.pop();

        // std::thread::id threadId = std::this_thread::get_id();
        std::string logFilePath = getLogFilenameForThread(threadId) + ".log";

        if (std::filesystem::exists(logFilePath) && std::filesystem::file_size(logFilePath) >= maxFileSize) {
            rotateLogFile(threadId);
            logFilePath = getLogFilenameForThread(threadId) + ".log"; // Update the file path after rotation
        }

        std::ofstream& logFile = getLogFile(threadId);
        logFile << logEntry << std::endl;
    }
}
