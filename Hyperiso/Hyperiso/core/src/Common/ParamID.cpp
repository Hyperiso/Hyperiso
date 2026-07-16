
#include "ParamID.h"

ParamId::ParamId() : block("NULL"), code(0) {}
ParamId::ParamId(const BlockName& block, const LhaID& code) : block(block), code(code) {}
ParamId::ParamId(ParameterType type, const BlockName& block, const LhaID& code) : type(type), block(block), code(code) {}

void ParamId::set_parameter_type(ParameterType type) { this->type = type; }

bool operator==(const ParamId& lhs, const ParamId& rhs) { 
    return lhs.type == rhs.type && lhs.block == rhs.block && lhs.code == rhs.code;
};

bool ParamId::operator<(const ParamId& other) const {
    if (type != other.type) return type < other.type;
    if (block != other.block) return block < other.block;
    return code < other.code;
}