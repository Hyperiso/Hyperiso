#ifndef IPARAMETER_RESOLVER_H
#define IPARAMETER_RESOLVER_H

#include <unordered_map>
#include <string>
#include <vector>
#include "Extractor.h"
#include "LhaID.h"

struct ResolvedParam {
    std::string block;
    LhaID code;
    bool is_complex;
    bool is_bsm;
};

class IParameterResolver {
public:
    virtual ~IParameterResolver() = default;

    virtual std::unordered_map<std::string, ResolvedParam>
    resolve(const std::vector<Extractor::Parameter>& params,
            bool modelIsSM) const = 0;

    virtual std::unique_ptr<IParameterResolver> clone() const = 0;
};

#endif
