cmake_minimum_required(VERSION 3.20)
project(hyperiso_examples_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GSL REQUIRED)
find_package(Hyperiso REQUIRED)  # usually installed in ~/.local/lib/cmake/Hyperiso or /usr/local/lib/cmake/Hyperiso

set(HYPERISO_CORE_LIB Hyperiso::CoreLib)
set(HYPERISO_PHYSICAL_MODEL_LIB Hyperiso::PhysicalModelLib)
set(HYPERISO_BUSINESS_LOGIC_LIB Hyperiso::BusinessLogicLib)

option(BUILD_STATISTIC_EXAMPLES "Build examples requiring Hyperiso::StatisticLib" ON)
set(HYPERISO_STATISTIC_EXAMPLES_AVAILABLE OFF)

if(BUILD_STATISTIC_EXAMPLES AND TARGET Hyperiso::StatisticLib)
    set(_hyperiso_stat_target_broken OFF)
    get_target_property(_hyperiso_stat_includes Hyperiso::StatisticLib INTERFACE_INCLUDE_DIRECTORIES)

    if(_hyperiso_stat_includes)
        foreach(_inc IN LISTS _hyperiso_stat_includes)
            # Generator expressions cannot be checked at configure time.
            if(_inc MATCHES "^\\$<")
                continue()
            endif()

            if(IS_ABSOLUTE "${_inc}" AND NOT EXISTS "${_inc}")
                message(WARNING
                    "Skipping Statistic examples because Hyperiso::StatisticLib exports a missing include directory: ${_inc}")
                set(_hyperiso_stat_target_broken ON)
            endif()
        endforeach()
    endif()

    if(NOT _hyperiso_stat_target_broken)
        set(HYPERISO_STATISTIC_EXAMPLES_AVAILABLE ON)
        set(HYPERISO_STATISTIC_LIB Hyperiso::StatisticLib)
    endif()
elseif(BUILD_STATISTIC_EXAMPLES)
    message(WARNING "Skipping Statistic examples because Hyperiso::StatisticLib was not exported by this Hyperiso installation.")
endif()

function(add_hyperiso_example target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE ${ARGN})
endfunction()

# -------------------------
# Core examples
# -------------------------
add_hyperiso_example(core_base Core/base_core_example.cpp ${HYPERISO_CORE_LIB})
add_hyperiso_example(core_block_logger Core/block_logger_example.cpp ${HYPERISO_CORE_LIB})
add_hyperiso_example(core_dependency_pruner Core/dependency_pruner_example.cpp ${HYPERISO_CORE_LIB})
add_hyperiso_example(core_parameter_provider Core/parameter_provider_example.cpp ${HYPERISO_CORE_LIB})
add_hyperiso_example(core_parameter_setter Core/parameter_setter_example.cpp ${HYPERISO_CORE_LIB})
add_hyperiso_example(core_preinit Core/preinit_core_example.cpp ${HYPERISO_CORE_LIB})
add_hyperiso_example(core_qcd_runner Core/qcd_runner_example.cpp ${HYPERISO_CORE_LIB})
add_hyperiso_example(dev_core_api_diagnostics Core/dev_examples/api_diagnostics_example.cpp ${HYPERISO_CORE_LIB})
add_hyperiso_example(dev_core_composite_param Core/dev_examples/composite_param_creator_example.cpp ${HYPERISO_CORE_LIB})

# -------------------------
# Wilson examples
# -------------------------
add_hyperiso_example(wilson_base Wilson/base_wilson_example.cpp ${HYPERISO_PHYSICAL_MODEL_LIB})
add_hyperiso_example(wilson_marty Wilson/wilson_with_marty_example.cpp ${HYPERISO_PHYSICAL_MODEL_LIB})
add_hyperiso_example(dev_wilson_bulk_dump Wilson/dev_examples/wilson_bulk_dump_example.cpp ${HYPERISO_PHYSICAL_MODEL_LIB})
add_hyperiso_example(dev_wilson_scale_scan Wilson/dev_examples/wilson_scale_scan_example.cpp ${HYPERISO_PHYSICAL_MODEL_LIB})
add_hyperiso_example(dev_wilson_dynamic_ids Wilson/dev_examples/wilson_dynamic_ids_example.cpp ${HYPERISO_PHYSICAL_MODEL_LIB})
add_hyperiso_example(dev_wilson_lambda_custom_group Wilson/dev_examples/wilson_lambda_custom_group_example.cpp ${HYPERISO_PHYSICAL_MODEL_LIB})

# -------------------------
# Observable examples
# -------------------------
add_hyperiso_example(observable_base Observable/base_observable_example.cpp ${HYPERISO_BUSINESS_LOGIC_LIB})
add_hyperiso_example(dev_observable_scan Observable/dev_examples/observable_parameter_scan_example.cpp ${HYPERISO_BUSINESS_LOGIC_LIB})
add_hyperiso_example(dev_observable_dynamic_ids Observable/dev_examples/observable_dynamic_ids_example.cpp ${HYPERISO_BUSINESS_LOGIC_LIB})
add_hyperiso_example(dev_observable_lambda_decay Observable/dev_examples/observable_lambda_decay_example.cpp ${HYPERISO_BUSINESS_LOGIC_LIB})

# -------------------------
# Statistic examples
# -------------------------
if(HYPERISO_STATISTIC_EXAMPLES_AVAILABLE)
    add_hyperiso_example(statistic_base Statistic/base_statistic_example.cpp ${HYPERISO_STATISTIC_LIB})
    add_hyperiso_example(dev_statistic_likelihood_scan Statistic/dev_examples/likelihood_scan_example.cpp ${HYPERISO_STATISTIC_LIB})
    add_hyperiso_example(dev_statistic_custom_observables Statistic/dev_examples/custom_observables_dependency_example.cpp ${HYPERISO_STATISTIC_LIB})
else()
    message(STATUS "Statistic examples are present in the source tree but are not added as CMake targets for this installation.")
endif()
