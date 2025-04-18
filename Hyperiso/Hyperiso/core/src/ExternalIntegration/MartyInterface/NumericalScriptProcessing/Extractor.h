#ifndef EXTRACTOR_H
#define EXTRACTOR_H

#include <string>
#include <vector>
#include <regex>
#include <fstream>

class Extractor {
public:
    struct Parameter {
        std::string type;
        std::string name;
        bool complex;
    };

    static std::vector<Parameter> extract(const std::string& filename);
};

#endif // EXTRACTOR_H
