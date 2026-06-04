#include "DependantBlockInfoProvider.h"

bool DependantBlockInfoProvider::is_dependent_block(ParameterType type, const std::string& block_name) const {
    return Parameters::GetInstance(type)->is_dependent_block(BlockName(block_name));
}

std::vector<std::string> DependantBlockInfoProvider::get_source_blocks(
    ParameterType type,
    const std::string& block_name
) const {
    return Parameters::GetInstance(type)->get_source_blocks(BlockName(block_name));
}

std::vector<std::string> DependantBlockInfoProvider::get_dependent_blocks(
    ParameterType type,
    const std::string& block_name
) const {
    return Parameters::GetInstance(type)->get_dependent_blocks(BlockName(block_name));
}

std::vector<std::string> DependantBlockInfoProvider::get_all_source_blocks(
    ParameterType type,
    const std::string& block_name
) const {
    return Parameters::GetInstance(type)->get_all_source_blocks(BlockName(block_name));
}

std::vector<std::string> DependantBlockInfoProvider::get_all_dependent_blocks(
    ParameterType type,
    const std::string& block_name
) const {
    return Parameters::GetInstance(type)->get_all_dependent_blocks(BlockName(block_name));
}
