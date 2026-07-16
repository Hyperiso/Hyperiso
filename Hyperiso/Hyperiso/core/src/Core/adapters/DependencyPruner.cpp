#include "DependencyPruner.h"

void DependencyPruner::reattach_block(ParameterType tp, const BlockName &block_name) {
    Parameters::GetInstance(tp)->reattach_block(block_name);
}

void DependencyPruner::detach_block(ParameterType tp, const BlockName &block_name) {
    Parameters::GetInstance(tp)->detach_block(block_name);
}

void DependencyPruner::reattach_parameter(ParameterType tp, const BlockName &block_name, const LhaID &id) {
    Parameters::GetInstance(tp)->reattach_param(block_name, id);
}

void DependencyPruner::detach_parameter(ParameterType tp, const BlockName &block_name, const LhaID &id) {
    Parameters::GetInstance(tp)->detach_param(block_name, id);
}
