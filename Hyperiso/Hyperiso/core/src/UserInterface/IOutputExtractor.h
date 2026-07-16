#ifndef IOUTPUT_EXTRACTOR_H
#define IOUTPUT_EXTRACTOR_H

#include <vector>
#include <unordered_map>

#include "ObjectsOutputs.h"

class IOutputExtractor {
public:
    virtual ~IOutputExtractor() = default;

    virtual std::vector<std::string> schema(const OutputSpec& spec) const = 0;

    virtual void extract(std::unordered_map<std::string, Value>& outY,
                         const OutputSpec& spec) const = 0;
};

#endif