#ifndef __IBLOCKCOMPOSER_H__
#define __IBLOCKCOMPOSER_H__

#include "Include.h"
#include "CompositeParamCreator.h"

class IBlockComposer {
public:
    virtual ~IBlockComposer() = default;

    virtual void compose(const std::string&, const std::unordered_map<ParameterType, std::vector<std::string>>&, const DepUpdateFunc&) = 0;
    virtual void update(const std::string&) = 0;
};

#endif // __IBLOCKCOMPOSER_H__
