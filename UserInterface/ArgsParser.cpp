#include "ArgsParser.h"


Argument::Argument(const std::string& longName, const std::string& shortName, const std::string& helpText,
                   ArgType type, bool isRequired, bool allowsMultiple,
                   std::optional<std::string> defaultValue, bool isPositional)
    : longName(longName), shortName(shortName), helpText(helpText), type(type),
      isRequired(isRequired), allowsMultiple(allowsMultiple),
      defaultValue(defaultValue), isPositional(isPositional) {}

void Argument::addValidator(std::shared_ptr<IValidator> validator) {
    validators.push_back(validator);
}

ArgType Argument::getType() const {
    return type;
}

bool Argument::validate(const std::string& value) const {
    for (const auto& validator : validators) {
        if (!validator->validate(value)) {
            throw std::invalid_argument(validator->errorMessage());
        }
    }
    return true;
}

std::string Argument::getLongName() const {
    return longName;
}

std::string Argument::getShortName() const {
    return shortName;
}

std::string Argument::getHelpText() const {
    return helpText;
}

bool Argument::isRequiredArg() const {
    return isRequired;
}

bool Argument::allowsMultipleValues() const {
    return allowsMultiple;
}

bool Argument::isPositionalArg() const {
    return isPositional;
}

std::optional<std::string> Argument::getDefaultValue() const {
    return defaultValue;
}

ArgumentBuilder& ArgumentBuilder::setType(ArgType argType) {
    type = argType;
    return *this;
}

ArgumentBuilder& ArgumentBuilder::setLongName(const std::string& name) {
    longName = name;
    return *this;
}

ArgumentBuilder& ArgumentBuilder::setShortName(const std::string& name) {
    shortName = name;
    return *this;
}

ArgumentBuilder& ArgumentBuilder::setHelpText(const std::string& text) {
    helpText = text;
    return *this;
}

ArgumentBuilder& ArgumentBuilder::setRequired(bool required) {
    isRequired = required;
    return *this;
}

ArgumentBuilder& ArgumentBuilder::setDefaultValue(const std::string& value) {
    defaultValue = value;
    return *this;
}

ArgumentBuilder& ArgumentBuilder::setAllowsMultiple(bool multiple) {
    allowsMultiple = multiple;
    return *this;
}

ArgumentBuilder& ArgumentBuilder::addValidator(std::shared_ptr<IValidator> validator) {
    validators.push_back(validator);
    return *this;
}

ArgumentBuilder& ArgumentBuilder::setPositional(bool positional) {
    isPositional = positional;
    return *this;
}

Argument ArgumentBuilder::build() {
    Argument arg = Argument(longName, shortName, helpText, type, isRequired, allowsMultiple, defaultValue, isPositional);
    if (!validators.empty()){
        for (auto &validator : validators) {
            arg.addValidator(validator);
        }
    }
    return arg;
}

void ArgParser::addArgument(const Argument& arg) {
    if (arg.isPositionalArg()) {
        positional_arguments.push_back(arg);
    }
    arguments.push_back(arg);
}

void ArgParser::parse(int argc, char* argv[]) {
    size_t positionalIndex = 0;
    positionalValues.clear();
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--", 0) == 0 || arg.rfind("-", 0) == 0) {
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
                if (!values.size()) {
                    throw std::invalid_argument("Please put at least one argument for : " + arg);
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

    for (const auto& arg : arguments) {
        if (parsedValues.find(arg.getLongName()) == parsedValues.end() && arg.getDefaultValue()) {
            parsedValues[arg.getLongName()] = {*arg.getDefaultValue()};
        }
}
}

std::string ArgParser::getValue(const std::string& name) const {
    auto it = parsedValues.find(name);
    if (it != parsedValues.end() && !it->second.empty()) {
        return it->second[0];
    }
    for (const auto& arg : arguments) {
        if (arg.getLongName() == name && arg.getDefaultValue()) {
            return *arg.getDefaultValue();
        }
    }
    throw std::invalid_argument("Value for argument " + name + " not provided.");
}

std::vector<std::string> ArgParser::getValues(const std::string& name) const {
    auto it = parsedValues.find(name);
    if (it != parsedValues.end()) {
        return it->second;
    }
    throw std::invalid_argument("Values for argument " + name + " not provided.");
}

std::vector<std::string> ArgParser::getPositionalValues() const {
    return positionalValues;
}

void ArgParser::displayHelp() const {
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
