#ifndef __WILSONPARAMCOMPOSER_H__
#define __WILSONPARAMCOMPOSER_H__

#include "IBlockComposer.h"
#include "CompositeParamCreator.h"

class WilsonParamComposer : public IBlockComposer {
public:
    void compose_block(const std::string& block_name, const std::unordered_map<ParameterType, std::vector<std::string>>& source_names, const DepUpdateFunc& update_func) override;
    void compose_parameter(const ParamId&, const std::unordered_set<ParamId>&, const DepParamUpdateFunc&) override;
    void remove_block(const std::string&) override;
    void update(const std::string& block_name) override;
    void remove_all_composed_blocks() override;

};

#endif // __WILSONPARAMCOMPOSER_H__
