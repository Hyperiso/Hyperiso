#include <cassert>
#include <iostream>
#include "ArgsParser.h"

void test_argparser() {
    {
        ArgParser parser;
        parser.addArgument(ArgumentBuilder()
                               .setLongName("level")
                               .setShortName("l")
                               .setHelpText("An integer level between 1 and 5")
                               .setType(ArgType::INT)
                               .addValidator(std::make_shared<RangeValidator>(1, 5, ArgType::INT))
                               .setRequired(true)
                               .build());

        {
            const char* argv[] = {"program", "--level", "3"};
            parser.parse(3, const_cast<char**>(argv));
            assert(parser.getValue("level") == "3");
        }

        try {
            const char* argv[] = {"program", "--level", "6"};
            parser.parse(3, const_cast<char**>(argv));
            assert(false);
        } catch (const std::exception& e) {
            assert(std::string(e.what()).find("Value must be between") != std::string::npos);
        }
    }

    {
        ArgParser parser;
        parser.addArgument(ArgumentBuilder()
                               .setLongName("color")
                               .setShortName("c")
                               .setHelpText("A color (red, green, or blue)")
                               .setType(ArgType::STRING)
                               .addValidator(std::make_shared<AllowedValuesValidator>(
                               std::vector<std::string>{"red", "green", "blue"}, ArgType::STRING))
                               .setRequired(true)
                               .build());

        {
            const char* argv[] = {"program", "--color", "green"};
            parser.parse(3, const_cast<char**>(argv));
            assert(parser.getValue("color") == "green");
        }

        try {
            const char* argv[] = {"program", "--color", "yellow"};
            parser.parse(3, const_cast<char**>(argv));
            assert(false);
        } catch (const std::exception& e) {
            assert(std::string(e.what()).find("Value must be one of") != std::string::npos);
        }
    }

    {
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

        {
            const char* argv[] = {"program", "file1.txt", "file2.txt"};
            parser.parse(3, const_cast<char**>(argv));
            auto positional = parser.getPositionalValues();
            assert(positional.size() == 2);
            assert(positional[0] == "file1.txt");
            assert(positional[1] == "file2.txt");
        }

        try {
            const char* argv[] = {"program", "file1.txt"};
            parser.parse(2, const_cast<char**>(argv));
            assert(false);
        } catch (const std::invalid_argument& e) {
            assert(std::string(e.what()).find("Missing required positional argument") != std::string::npos);
        }
    }

    {
        ArgParser parser;
        parser.addArgument(ArgumentBuilder()
                               .setLongName("numbers")
                               .setShortName("n")
                               .setHelpText("A list of numbers")
                               .setType(ArgType::INT)
                               .setAllowsMultiple(true)
                               .setRequired(true)
                               .build());

        {
            const char* argv[] = {"program", "--numbers", "1", "2", "3"};
            parser.parse(5, const_cast<char**>(argv));
            auto values = parser.getValues("numbers");
            assert(values.size() == 3);
            assert(values[0] == "1");
            assert(values[1] == "2");
            assert(values[2] == "3");
        }

        try {
            const char* argv[] = {"program", "--numbers"};
            parser.parse(2, const_cast<char**>(argv));
            assert(false);
        } catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
            assert(std::string(e.what()).find("Please put at least one argument") != std::string::npos);
        }
    }

    {
        ArgParser parser;

        try {
            const char* argv[] = {"program", "--unknown", "value"};
            parser.parse(3, const_cast<char**>(argv));
            assert(false);
        } catch (const std::exception& e) {
            assert(std::string(e.what()).find("Unknown option") != std::string::npos);
        }
    }

    {
        ArgParser parser;
        parser.addArgument(ArgumentBuilder()
                               .setLongName("help")
                               .setShortName("h")
                               .setHelpText("Show this help message")
                               .setType(ArgType::STRING)
                               .setRequired(false)
                               .build());

        parser.displayHelp();
    }

    std::cout << "All tests passed!\n";
}

int main() {
    test_argparser();
    return 0;
}
