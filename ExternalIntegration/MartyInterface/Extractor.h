#ifndef EXTRACTOR_H
#define EXTRACTOR_H

#include <string>
#include <vector>

class Extractor {
public:
    struct Parameter {
        std::string type;
        std::string name;
    };

    std::vector<Parameter> extract(const std::string& filename);
};

#endif // EXTRACTOR_H
