#include <memory>

#include "InterpretedParam.h"
#include "IBlockComposer.h"
#include "ICoreAPI.h"
#include "IMartyWilsonProxy.h"
#include "IParameterProxy.h"
#include "config.hpp"

struct WilsonGroupAdapterConfig {

    WilsonGroupAdapterConfig(std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy, std::shared_ptr<IBlockComposer> ibc,
    std::shared_ptr<ICoreAPI<bool>> use_marty, std::shared_ptr<ICoreAPI<std::string>> marty_model_name, std::shared_ptr<ICoreAPI<fs::path>> marty_model_path,
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy = nullptr) 
    : wilson_proxy(wilson_proxy), iblock_c(ibc), use_marty(use_marty), marty_model_name(marty_model_name), marty_model_path(marty_model_path), marty_proxy(marty_proxy) {}

    std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy;
    std::shared_ptr<IBlockComposer> iblock_c;
    std::shared_ptr<ICoreAPI<bool>> use_marty;
    std::shared_ptr<ICoreAPI<std::string>> marty_model_name;
    std::shared_ptr<ICoreAPI<fs::path>> marty_model_path;
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy = nullptr;

    fs::path sm_path = fs::path(std::string(project_assets_root.data())+"input_files/marty_model/sm.h");
};