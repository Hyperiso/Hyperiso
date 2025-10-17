#pragma once
#include <memory>

class IBlockComposer;

class IWilsonParameterHelper {
public:

    IWilsonParameterHelper(std::shared_ptr<IBlockComposer> iblock_c) : iblock_c(iblock_c) {}

    virtual ~IWilsonParameterHelper() = default;

    virtual void init(int gen) = 0;

    virtual void cleanup() = 0;

    bool is_init() const {return initialized;}

protected:
    virtual void init_scale_independent_block(int gen) = 0;
    virtual void init_matching_block() = 0;
    virtual void init_running_block() = 0;

    std::shared_ptr<IBlockComposer> iblock_c;
    bool initialized{false};
};