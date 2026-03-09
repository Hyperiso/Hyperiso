#include "StatDependencyPruner.h"

void StatDependencyPruner::reattach_block(ParameterType tp, const std::string &block_name) {
    dp.reattach_block(tp, block_name);
}

void StatDependencyPruner::detach_block(ParameterType tp, const std::string &block_name) {
    dp.detach_block(tp, block_name);
}

void StatDependencyPruner::reattach_parameter(ParameterType tp, const std::string &block_name, const LhaID &id) {
    dp.reattach_parameter(tp, block_name, id);
}

void StatDependencyPruner::detach_parameter(ParameterType tp, const std::string &block_name, const LhaID &id) {
    dp.detach_parameter(tp, block_name, id);
}
