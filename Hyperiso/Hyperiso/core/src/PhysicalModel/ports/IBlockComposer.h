#ifndef __IBLOCKCOMPOSER_H__
#define __IBLOCKCOMPOSER_H__

#include "Include.h"
#include "CompositeParamCreator.h"

class IBlockComposer {
public:
    virtual ~IBlockComposer() = default;

    virtual void compose_block(const std::string&, const std::unordered_map<ParameterType, std::vector<std::string>>&, const DepUpdateFunc&) = 0;
    virtual void compose_parameter(const ParamId&, const std::unordered_set<ParamId>&, const DepParamUpdateFunc&) = 0;
    virtual void remove_block(const std::string&) = 0;
    virtual void update(const std::string&) = 0;
};

#endif // __IBLOCKCOMPOSER_H__
