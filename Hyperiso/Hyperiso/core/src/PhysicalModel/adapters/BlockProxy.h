#ifndef BLOCK_PROXY_H
#define BLOCK_PROXY_H

#include "BlockProvider.h"
#include "APIAdapter.h"


class BlockProxy {
public:
    BlockProxy() {bp = BlockProvider();}
    bool exists(const std::string& blockname, ParameterType pt);
    void log_all_blocks(ParameterType pt);
    void log_block(ParameterType pt, const std::string& blockname);
    std::unordered_set<BlockName> get_block_list(ParameterType pt);

private:
    BlockProvider bp;
    APIAdapter api;
};

#endif