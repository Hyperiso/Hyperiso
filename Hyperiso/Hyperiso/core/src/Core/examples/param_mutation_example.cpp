/**
 * @example param_mutation_example.cpp
 * @brief Example usage of ParameterSetter and ParameterShifter.
 *
 * This example demonstrates how to modify parameter values
 * either by setting them directly or by applying a shift.
 *
 * @code
 * #include "ParameterSetter.h"
 * #include "ParameterShifter.h"
 * #include "Parameters.h"
 * 
 * int main() {
 *     // Assume the Parameters system is already initialized elsewhere
 *     
 *     ParamId pid{ParameterType::SM, "MASS", 24}; // W boson mass
 *     
 *     ParameterSetter setter;
 *     setter.mutate(pid, 80.385); // Set the value directly
 *     setter.change_mode(pid, ParameterMode::SHIFTABLE); // Allow shifts
 *     
 *     ParameterShifter shifter;
 *     shifter.mutate(pid, 0.1); // Apply a +0.1 GeV shift to the W mass
 *     
 *     return 0;
 * }
 * @endcode
 */
