#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <sstream>
#include <functional>
#include <optional>
#include <algorithm>
#include <memory>

enum class ArgType {
    STRING,
    INT,
    DOUBLE
};

class IValidator {
public:
    virtual bool validate(const std::string& value) const = 0;
    virtual std::string errorMessage() const = 0;
    virtual ~IValidator() = default;
};

template <typename T>
class Validator : public IValidator {
public:
    virtual bool validate(const T& value) const = 0;
};

class RangeValidator : public Validator<int> {
    int minValue, maxValue;

public:
    RangeValidator(int min, int max) : minValue(min), maxValue(max) {}

    bool validate(const int& value) const override {
        return value >= minValue && value <= maxValue;
    }

    bool validate(const std::string& value) const override {
        try {
            int intValue = std::stoi(value);
            return validate(intValue);
        } catch (...) {
            return false;
        }
    }

    std::string errorMessage() const override {
        return "Value must be between " + std::to_string(minValue) + " and " + std::to_string(maxValue);
    }
};

template <typename T>
class AllowedValuesValidator : public Validator<T> {
    std::vector<T> allowedValues;

public:
    explicit AllowedValuesValidator(const std::vector<T>& values) : allowedValues(values) {}

    bool validate(const T& value) const override {
        return std::find(allowedValues.begin(), allowedValues.end(), value) != allowedValues.end();
    }

    bool validate(const std::string& value) const override {
        try {
            if constexpr (std::is_same<T, int>::value) {
                int intValue = std::stoi(value);
                return validate(intValue);
            } else if constexpr (std::is_same<T, double>::value) {
                double doubleValue = std::stod(value);
                return validate(doubleValue);
            } else {
                return validate(value);
            }
        } catch (...) {
            return false;
        }
    }

    std::string errorMessage() const override {
        return "Value must be one of: " + joinAllowedValues();
    }

private:
    std::string joinAllowedValues() const {
        std::ostringstream oss;
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
             std::optional<std::string> defaultValue, bool isPositional)
        : longName(longName), shortName(shortName), helpText(helpText), type(type),
          isRequired(isRequired), allowsMultiple(allowsMultiple),
          defaultValue(defaultValue), isPositional(isPositional) {}

    void addValidator(std::shared_ptr<IValidator> validator) {
        validators.push_back(validator);
    }

    ArgType getType() const { return type; }

    bool validate(const std::string& value) const {
        try {
            if (type == ArgType::INT) {
                int intValue = std::stoi(value);
                for (const auto& validator : validators) {
                    auto intValidator = std::dynamic_pointer_cast<Validator<int>>(validator);
                    std::cout << intValidator << std::endl;
                    if (intValidator && !intValidator->validate(intValue)) {
                        throw std::invalid_argument(intValidator->errorMessage());
                    }
                }
            } else if (type == ArgType::DOUBLE) {
                double doubleValue = std::stod(value);
                for (const auto& validator : validators) {
                    auto doubleValidator = std::dynamic_pointer_cast<Validator<double>>(validator);
                    if (doubleValidator && !doubleValidator->validate(doubleValue)) {
                        throw std::invalid_argument(doubleValidator->errorMessage());
                    }
                }
            } else {
                for (const auto& validator : validators) {
                    auto stringValidator = std::dynamic_pointer_cast<Validator<std::string>>(validator);
                    if (stringValidator && !stringValidator->validate(value)) {
                        throw std::invalid_argument(stringValidator->errorMessage());
                    }
                }
            }
        } catch (const std::exception& e) {
            throw std::invalid_argument("Validation failed for argument: " + longName + ". " + e.what());
        }
        return true;
    }

    std::string getLongName() const { return longName; }
    std::string getShortName() const { return shortName; }
    std::string getHelpText() const { return helpText; }
    bool isRequiredArg() const { return isRequired; }
    bool allowsMultipleValues() const { return allowsMultiple; }
    bool isPositionalArg() const { return isPositional; }
    std::optional<std::string> getDefaultValue() const { return defaultValue; }
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
    ArgumentBuilder& setType(ArgType argType) {
        type = argType;
        return *this;
    }
    ArgumentBuilder& setLongName(const std::string& name) {
        longName = name;
        return *this;
    }

    ArgumentBuilder& setShortName(const std::string& name) {
        shortName = name;
        return *this;
    }

    ArgumentBuilder& setHelpText(const std::string& text) {
        helpText = text;
        return *this;
    }

    ArgumentBuilder& setRequired(bool required) {
        isRequired = required;
        return *this;
    }

    ArgumentBuilder& setDefaultValue(const std::string& value) {
        defaultValue = value;
        return *this;
    }

    ArgumentBuilder& setAllowsMultiple(bool multiple) {
        allowsMultiple = multiple;
        return *this;
    }

    ArgumentBuilder& addValidator(std::shared_ptr<IValidator> validator) {
        validators.push_back(validator);
        return *this;
    }

    ArgumentBuilder& setPositional(bool positional) {
        isPositional = positional;
        return *this;
    }

    Argument build() {
        Argument arg = Argument(longName, shortName, helpText, type, isRequired, allowsMultiple, defaultValue, isPositional);
        if (!validators.empty()){
            for (auto &validator : validators) {
                arg.addValidator(validator);
            }
        }
        return arg;
    }
};

class ArgParser {
    std::vector<Argument> arguments;
    std::vector<Argument> positional_arguments;
    std::map<std::string, std::vector<std::string>> parsedValues;
    std::vector<std::string> positionalValues;

public:
    void addArgument(const Argument& arg) {
        if (arg.isPositionalArg()) {
            positional_arguments.push_back(arg);
        }
        arguments.push_back(arg);
    }

    void parse(int argc, char* argv[]) {
        size_t positionalIndex = 0;

        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];

            if (arg.rfind("--", 0) == 0 || arg.rfind("-", 0) == 0) { // Options
                std::string key = arg.substr(arg[1] == '-' ? 2 : 1);
                auto it = std::find_if(arguments.begin(), arguments.end(), [&](const Argument& a) {
                    return a.getLongName() == key || a.getShortName() == key;
                });

                if (it == arguments.end()) {
                    throw std::invalid_argument("Unknown option: " + arg);
                }

                Argument& option = *it;

                if (option.allowsMultipleValues()) {
                    std::vector<std::string> values;
                    while (i + 1 < argc && argv[i + 1][0] != '-') {
                        std::string value = argv[++i];
                        option.validate(value);
                        values.push_back(value);
                    }
                    parsedValues[option.getLongName()] = values;
                } else {
                    if (i + 1 < argc && argv[i + 1][0] != '-') {
                        std::string value = argv[++i];
                        option.validate(value);
                        parsedValues[option.getLongName()] = {value};
                    } else {
                        throw std::invalid_argument("Missing value for option: " + key);
                    }
                }
            } else {
                if (positionalIndex >= positional_arguments.size()) {
                    throw std::invalid_argument("Unexpected positional argument: " + arg);
                }

                Argument& positionalArg = positional_arguments[positionalIndex++];
                positionalArg.validate(arg);
                positionalValues.push_back(arg);
            }
        }

        for (const auto& arg : positional_arguments) {
            if (arg.isRequiredArg() && positionalValues.size() < positional_arguments.size()) {
                throw std::invalid_argument("Missing required positional argument: " + arg.getLongName());
            }
        }
    }

    std::string getValue(const std::string& name) const {
        auto it = parsedValues.find(name);
        if (it != parsedValues.end() && !it->second.empty()) return it->second[0];
        throw std::invalid_argument("Value for argument " + name + " not provided.");
    }

    std::vector<std::string> getValues(const std::string& name) const {
        auto it = parsedValues.find(name);
        if (it != parsedValues.end()) return it->second;
        throw std::invalid_argument("Values for argument " + name + " not provided.");
    }

    std::vector<std::string> getPositionalValues() const {
        return positionalValues;
    }

    void displayHelp() const {
        for (const auto& arg : arguments) {
            if (arg.isPositionalArg()) {
                std::cout << arg.getLongName() << ": " << arg.getHelpText() << " (Positional)\n";
            } else {
                std::cout << "--" << arg.getLongName() << ", -" << arg.getShortName()
                          << ": " << arg.getHelpText();
                if (arg.getDefaultValue()) {
                    std::cout << " (default: " << *arg.getDefaultValue() << ")";
                }
                std::cout << "\n";
            }
        }
    }
};

int main(int argc, char* argv[]) {
    ArgParser parser;
    parser.addArgument(ArgumentBuilder()
                            .setLongName("input")
                            .setHelpText("Input file (positional)")
                            .setType(ArgType::STRING)
                            .setPositional(true)
                            .setRequired(true)
                            .build());

    parser.addArgument(ArgumentBuilder()
                            .setLongName("output")
                            .setHelpText("Output file (positional)")
                            .setType(ArgType::STRING)
                            .setPositional(true)
                            .setRequired(true)
                            .build());

    try {
        parser.parse(argc, argv);
        std::cout << "--input: " << parser.getPositionalValues()[0] << "\n";
        std::cout << "--output: " << parser.getPositionalValues()[1] << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        parser.displayHelp();
        return 1;
    }

    return 0;
}
