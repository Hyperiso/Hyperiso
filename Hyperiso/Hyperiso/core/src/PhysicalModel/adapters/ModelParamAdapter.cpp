#include "ModelParamAdapter.h"

double ModelParamAdapter::operator()(std::string block, int code) {
    return pp(block, code);
}