#ifndef __WILSONPARAMCOMPOSER_H__
#define __WILSONPARAMCOMPOSER_H__

#include "IBlockComposer.h"
#include "CompositeParamCreator.h"

class WilsonParamComposer : IBlockComposer {
public:
    void compose(const std::string& block_name, const std::unordered_map<ParameterType, std::vector<std::string>>& source_names, const DepUpdateFunc& update_func);
    void update(const std::string& block_name);

};

#endif // __WILSONPARAMCOMPOSER_H__
