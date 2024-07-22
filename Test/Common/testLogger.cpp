#include "Logger.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <exception>
#include <csignal>
#include <cstdlib>
#include <sys/wait.h>

void testTerminate() {
    std::cout << "std::terminate called as expected." << std::endl;
    // std::_Exit(0);
    throw std::runtime_error("yeah");
    abort();
}

std::string readLogFile(const std::string& filepath) {
    std::ifstream file(filepath);
    std::ostringstream ss;
    ss << file.rdbuf();
    return ss.str();
}

void testBasicLogging() {
    Logger* logger = Logger::getInstance();
    logger->setLogDirectory("test_logs");
    logger->setLevel(Logger::LogLevel::TRACE);
    
    LOG_INFO("This is an info log message");
    LOG_WARN("This is a warning log message");
    LOG_DEBUG("This is a debug log message");

    std::this_thread::sleep_for(std::chrono::seconds(1));

    std::string logFile = logger->getLogFilenameForThread(logger->threadId)  + ".log";
    std::string logContent = readLogFile(logFile);
    LOG_INFO(logFile);
    assert(logContent.find("This is an info log message") != std::string::npos);
    assert(logContent.find("This is a warning log message") != std::string::npos);
    assert(logContent.find("This is a debug log message") != std::string::npos);

    std::cout << "Basic logging test passed!" << std::endl;
}

void testLogError() {
    Logger* logger = Logger::getInstance();
    logger->setLogDirectory("test_logs");
    logger->setLevel(Logger::LogLevel::TRACE);

    auto originalTerminateHandler = std::set_terminate(testTerminate);
    try{
            LOG_ERROR("RuntimeError", "This is an error log message");
    }
    catch (...) {
        std::cout << "wow";
    }

    std::set_terminate(originalTerminateHandler);

    std::this_thread::sleep_for(std::chrono::seconds(1));

    std::string logFile = logger->getLogFilenameForThread(std::this_thread::get_id()) + ".log";
    std::string logContent = readLogFile(logFile);

    assert(logContent.find("This is an error log message") != std::string::npos);

    std::cout << "Error logging test passed!" << std::endl;
}


void testLogRotation() {
    Logger* logger = Logger::getInstance();
    logger->setLogDirectory("test_logs");
    logger->setLevel(Logger::LogLevel::TRACE);
    logger->setLogDirectory("test_logs", 1024); 

    for (int i = 0; i < 100; ++i) {
        LOG_INFO("This is a test log message number ", i);
    }

    std::this_thread::sleep_for(std::chrono::seconds(1));

    std::string logFile = logger->getLogFilenameForThread(std::this_thread::get_id()) + ".log";
    assert(std::filesystem::exists(logFile));

    size_t fileSize = std::filesystem::file_size(logFile);
    assert(fileSize < 1024); 

    std::cout << "Log rotation test passed!" << std::endl;
}

int main() {
    std::filesystem::remove_all("test_logs");

    testBasicLogging();
    // testLogError();
    testLogRotation();

    return 0;
}
