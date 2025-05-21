#include "BlockProxy.h"

bool BlockProxy::exists(const std::string& blockname, ParameterType pt) {
    return bp.exists(blockname, pt);
}

void BlockProxy::log_all_blocks(ParameterType pt) {
    bp.log_all_blocks(pt);
}

void BlockProxy::log_block(ParameterType pt, const std::string& blockname) {
    bp.log_block(pt, blockname);
}