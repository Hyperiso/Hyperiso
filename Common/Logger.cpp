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
    loggingThread = std::thread(&Logger::processQueue, this);
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
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %X", std::localtime(&now));
    return std::string(buf);
}

std::ofstream& Logger::getLogFile() {
    std::thread::id threadId = std::this_thread::get_id();
    std::lock_guard<std::mutex> lock(logFileMutex);

    auto it = logFiles.find(threadId);
    if (it == logFiles.end()) {
        std::ostringstream filename;
        filename << logDirectory << "/log_" << currentDateTime() << "_thread_" << threadId << ".log";
        std::string filepath = filename.str();
        logFiles[threadId].open(filepath, std::ios_base::app);
        if (!logFiles[threadId].is_open()) {
            throw std::runtime_error("Cannot open log file for thread");
        }
    }
    return logFiles[threadId];
}

void Logger::rotateLogFile(std::thread::id threadId) {
    std::lock_guard<std::mutex> lock(logFileMutex);
    std::ostringstream filename;
    filename << logDirectory << "/log_" << currentDateTime() << "_thread_" << threadId << ".log";
    std::string filepath = filename.str();

    if (std::filesystem::exists(filepath) && std::filesystem::file_size(filepath) >= maxFileSize) {
        std::ofstream& logFile = logFiles[threadId];
        logFile.close();
        std::ostringstream newFilename;
        newFilename << logDirectory << "/log_" << currentDateTime() << "_thread_" << threadId << ".log";
        std::filesystem::rename(filepath, newFilename.str());
        logFile.open(filepath, std::ios_base::app);
        if (!logFile.is_open()) {
            throw std::runtime_error("Cannot open log file after rotation for thread");
        }
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

            std::thread::id threadId = std::this_thread::get_id();
            rotateLogFile(threadId);

            std::ofstream& logFile = getLogFile();
            logFile << logEntry << std::endl;

            lock.lock();
        }
    }

    // Drain the queue when exiting
    while (!logQueue.empty()) {
        std::string logEntry = logQueue.front();
        logQueue.pop();

        std::thread::id threadId = std::this_thread::get_id();
        rotateLogFile(threadId);

        std::ofstream& logFile = getLogFile();
        logFile << logEntry << std::endl;
    }
}
