#ifndef __FREEZER_H__
#define __FREEZER_H__

#include "IFreezer.h"
#include "Include.h"
#include "Parameters.h"

class Freezer : public IFreezer<ParameterType, std::string, ParamId> {
public:
    static void freeze(const ParameterType&, const std::string&);
    static void freeze(const ParamId&);
    static void unfreeze(const ParameterType&, const std::string&);
    static void unfreeze(const ParamId&);
};

#endif // __FREEZER_H__
