#ifndef MARTY_INTERFACE_H
#define MARTY_INTERFACE_H

#include <memory>
#include <filesystem>
#include <iostream>
#include "config.hpp"
#include <algorithm>
#include "FileNameManager.h"
#include "GeneralModelModifier.h"
#include "CodeGenerator.h"
#include "GppCompilerStrategy.h"
#include "MakeCompilerStrategy.h"
#include "Logger.h"
#include "GeneralNumModelModifier.h"
#include "ICoreAPI.h"
#include "IMartyParameterProxy.h"

class MartyInterface {
public:

    MartyInterface(std::shared_ptr<ICoreAPI<Model>> core_api, std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> sm_proxy, std::shared_ptr<IInterpreterPortsFactory> ports, std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> bsm_proxy = nullptr) : core_api(core_api), param_proxy_sm(sm_proxy), ports(ports), param_proxy_bsm(bsm_proxy)
     {}
    
    MartyInterface();

    void generate(std::string wilson, std::string model, std::string model_path);
    void compile_run(std::string wilson, std::string model);
    void generate_numlib(std::string wilson, std::string model, double Q_match);
    void compile_run_libs(std::string wilson, std::string model, double Q_match);

    void calculate(std::string wilson, std::string model, double Q_match, std::string model_path, bool new_params=false);

    std::unordered_set<Interpreter::InterpretedParam> get_dependencies(std::string wilson);
    std::set<std::string> get_special_blocks();

private:

    bool already_run(std::string&& outputBinary);
    std::string output_binary_name(std::string& wilson, std::string& model);
    std::set<std::string> specials_block {"KIN", "WEIN", "Finite", "REGPROP", "BETA"};
    std::string num_file_path{};
    std::unordered_map<std::string, std::unordered_set<Interpreter::InterpretedParam>> dependencies;

    std::shared_ptr<ICoreAPI<Model>> core_api;
    std::shared_ptr<IInterpreterPortsFactory> ports;
    std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> param_proxy_sm;
    std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> param_proxy_bsm = nullptr;
};

#endif
