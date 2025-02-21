#ifndef ARG_PARSER_H
#define ARG_PARSER_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <sstream>
#include <algorithm>

enum class ArgType {
    STRING,
    INT,
    DOUBLE
};

class IValidator {
protected:
    ArgType type;
public:
    explicit IValidator(ArgType argType) : type(argType) {}
    virtual bool validate(const std::string& value) const = 0;
    virtual std::string errorMessage() const = 0;
    virtual ~IValidator() = default;
};

class RangeValidator : public IValidator {
    double minValue, maxValue;

public:
    RangeValidator(double min, double max, ArgType argType) 
        : IValidator(argType), minValue(min), maxValue(max) {}

    bool validate(const std::string& value) const override {
        try {
            switch (type) {
                case ArgType::INT: {
                    int intValue = std::stoi(value);
                    return intValue >= static_cast<int>(minValue) && intValue <= static_cast<int>(maxValue);
                }
                case ArgType::DOUBLE: {
                    double doubleValue = std::stod(value);
                    return doubleValue >= minValue && doubleValue <= maxValue;
                }
                case ArgType::STRING:
                    // RangeValidator is not logical for strings
                    return false;
                default:
                    return false;
            }
        } catch (...) {
            return false;
        }
    }

    std::string errorMessage() const override {
        return "Value must be between " + std::to_string(minValue) + " and " + std::to_string(maxValue);
    }
};

class AllowedValuesValidator : public IValidator {
    std::vector<std::string> allowedValues;

public:
    AllowedValuesValidator(const std::vector<std::string>& values, ArgType argType) 
        : IValidator(argType), allowedValues(values) {}

    bool validate(const std::string& value) const override {
        switch (type) {
            case ArgType::INT: {
                try {
                    int intValue = std::stoi(value);
                    return std::find(allowedValues.begin(), allowedValues.end(), std::to_string(intValue)) != allowedValues.end();
                } catch (...) {
                    return false;
                }
            }
            case ArgType::DOUBLE: {
                try {
                    double doubleValue = std::stod(value);
                    std::ostringstream oss;
                    oss << doubleValue;
                    return std::find(allowedValues.begin(), allowedValues.end(), oss.str()) != allowedValues.end();
                } catch (...) {
                    return false;
                }
            }
            case ArgType::STRING:
                return std::find(allowedValues.begin(), allowedValues.end(), value) != allowedValues.end();
            default:
                return false;
        }
    }

    std::string errorMessage() const override {
        std::ostringstream oss;
        oss << "Value must be one of: ";
        for (size_t i = 0; i < allowedValues.size(); ++i) {
            oss << allowedValues[i];
            if (i < allowedValues.size() - 1) oss << ", ";
        }
        return oss.str();
    }
};


class Argument {
    std::string longName, shortName, helpText;
    ArgType type;
    bool isRequired;
    bool allowsMultiple;
    std::optional<std::string> defaultValue;
    std::vector<std::shared_ptr<IValidator>> validators;
    bool isPositional;

public:
    Argument(const std::string& longName, const std::string& shortName, const std::string& helpText,
             ArgType type, bool isRequired, bool allowsMultiple,
             std::optional<std::string> defaultValue, bool isPositional);

    void addValidator(std::shared_ptr<IValidator> validator);
    ArgType getType() const;
    bool validate(const std::string& value) const;

    std::string getLongName() const;
    std::string getShortName() const;
    std::string getHelpText() const;
    bool isRequiredArg() const;
    bool allowsMultipleValues() const;
    bool isPositionalArg() const;
    std::optional<std::string> getDefaultValue() const;
};

class ArgumentBuilder {
    std::string longName, shortName, helpText;
    ArgType type = ArgType::STRING;
    bool isRequired = false;
    bool allowsMultiple = false;
    std::optional<std::string> defaultValue = std::nullopt;
    std::vector<std::shared_ptr<IValidator>> validators;
    bool isPositional = false;

public:
    ArgumentBuilder& setType(ArgType argType);
    ArgumentBuilder& setLongName(const std::string& name);
    ArgumentBuilder& setShortName(const std::string& name);
    ArgumentBuilder& setHelpText(const std::string& text);
    ArgumentBuilder& setRequired(bool required);
    ArgumentBuilder& setDefaultValue(const std::string& value);
    ArgumentBuilder& setAllowsMultiple(bool multiple);
    ArgumentBuilder& addValidator(std::shared_ptr<IValidator> validator);
    ArgumentBuilder& setPositional(bool positional);
    Argument build();
};

class ArgParser {
    std::vector<Argument> arguments;
    std::vector<Argument> positional_arguments;
    std::map<std::string, std::vector<std::string>> parsedValues;
    std::vector<std::string> positionalValues;

public:
    void addArgument(const Argument& arg);
    void parse(int argc, char* argv[]);
    std::string getValue(const std::string& name) const;
    std::vector<std::string> getValues(const std::string& name) const;
    std::vector<std::string> getPositionalValues() const;
    void displayHelp() const;
};


#endif // ARG_PARSER_H
