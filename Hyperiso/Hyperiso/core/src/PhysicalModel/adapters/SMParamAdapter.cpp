#include "SMParamAdapter.h"

double SMParamAdapter::operator()(std::string block, int code) {
    return pp(block, code);
}