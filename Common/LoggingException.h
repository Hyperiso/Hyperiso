#ifndef LOGGER_EXCEPTION_H
#define LOGGER_EXCEPTION_H

#include <exception>
#include <string>

class LoggingException : public std::exception {
public:
    explicit LoggingException(const std::string& message)
        : message_(message) {}

    virtual const char* what() const noexcept override {
        return message_.c_str();
    }

private:
    std::string message_;
};

#endif // LOGGER_EXCEPTION_H
