#ifndef BLOCK_PROVIDER_H
#define BLOCK_PROVIDER_H

#include "Include.h"
#include "Parameters.h"
#include "IBlockProvider.h"

class BlockProvider : public IBlockProvider<ParameterType, const std::string&> {
public:
    bool exists(const std::string& blockname, ParameterType) override;
    void log_all_blocks(ParameterType type) override;
    void log_block(ParameterType type, const std::string& blockname) override;

    
};

#endif