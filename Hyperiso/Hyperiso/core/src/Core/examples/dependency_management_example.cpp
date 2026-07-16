/**
 * @example examples/dependency_management_example.cpp
 * @brief Example usage of CompositeParamAdapter to manage parameter dependencies.
 *
 * This example demonstrates how to add a dependent block based on multiple source blocks
 * and how to trigger an update manually.
 *
 * @code
 * #include "CompositeParamAdapter.h"
 * #include "Parameters.h"
 * #include <iostream>
 *
 * int main() {
 *     // Assume the Parameters system is already initialized
 *     
 *     // Define the sources: blocks from different parameter types
 *     std::unordered_map<ParameterType, std::vector<std::string>> sources = {
 *         {ParameterType::SM, {"GAUGE", "MASS"}}
 *     };
 *
 *     // Define the recalculation function
 *     auto recalculate = [](const std::unordered_map<std::string, std::shared_ptr<Block>>& src_blocks,
 *                            std::shared_ptr<DependentBlock> dep_block) {
 *         double alpha_em = src_blocks.at("GAUGE")->retrieve(LhaID(1))->get_val();
 *         double mass_w  = src_blocks.at("MASS")->retrieve(LhaID(24))->get_val();
 *         
 *         dep_block->store_or_assign(LhaID(100), std::make_shared<Parameter>(ParamId("NEWBLOCK", 100), alpha_em + mass_w, 0., 0.));
 *     };
 *
 *     // Create the adapter and add the dependent block
 *     CompositeParamAdapter adapter;
 *     adapter.add_block_dependency("NEWBLOCK", sources, ParameterType::SM, recalculate);
 *
 *     // Trigger a manual update
 *     adapter.update_dependency("NEWBLOCK", ParameterType::SM);
 *
 *     // Access and print the new value
 *     auto value = Parameters::GetInstance(ParameterType::SM)->operator()("NEWBLOCK", LhaID(100));
 *     std::cout << "Combined value = " << value << std::endl;
 *
 *     return 0;
 * }
 * @endcode
 */
