#ifndef BLOCK_PROVIDER_H
#define BLOCK_PROVIDER_H

#include "Include.h"
#include "Parameters.h"

//TODO : port
class BlockProvider {
public:
    bool exists(const std::string& blockname, ParameterType);
    void log_all_blocks(ParameterType type);
};

#endif