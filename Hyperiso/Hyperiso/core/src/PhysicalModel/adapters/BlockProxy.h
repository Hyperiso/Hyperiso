#ifndef BLOCK_PROXY_H
#define BLOCK_PROXY_H

#include "BlockProvider.h"


class BlockProxy {
public:
    BlockProxy() {bp = BlockProvider();}
    bool exists(const std::string& blockname, ParameterType pt);
private:
    BlockProvider bp;
};

#endif