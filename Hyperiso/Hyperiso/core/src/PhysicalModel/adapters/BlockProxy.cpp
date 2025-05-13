#include "BlockProxy.h"

bool BlockProxy::exists(const std::string& blockname, ParameterType pt) {
    return bp.exists(blockname, pt);
}