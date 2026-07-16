#include "MartyInterface.h"
#include "GeneralNumModelModifier.h"
#include "Extractor.h"
#include "Interpreter.h"
#include "IncludeManager.hpp"
#include "LineProcessor.h"
#include "ModelWriter.h"
#include "SMParamSetter.h"
#include "FileNameManager.h"
#include "MemoryManager.h"
#include "HyperisoMaster.h"
#include "config.hpp"
#include "ModelAPI.h"
#include "MartyParameterProxy.h"
#include "DefaultInterpreterPortsFactory.h"

int main() {
    HyperisoMaster hyp = HyperisoMaster();

    std::string base_path = project_root.data();
    HyperisoConfig config;
    config.model = Model::SM;
    hyp.init(base_path + "Test/InputFiles/testinput_thdm.lha", config);

    std::shared_ptr<ModelAPI> model_api = std::make_shared<ModelAPI>();
    std::shared_ptr<MartyParameterProxy> sm_params = std::make_shared<MartyParameterProxy>(ParameterType::SM);
    std::shared_ptr<MartyParameterProxy> bsm_params = std::make_shared<MartyParameterProxy>(ParameterType::BSM);
    std::shared_ptr<IInterpreterPortsFactory> ports = std::make_shared<DefaultInterpreterPortsFactory>();
    MartyInterface MartyInterface(model_api, sm_params, ports, bsm_params);

    std::cout << sm_params->operator()("MASS", 24) << std::endl;

    MartyInterface.calculate("C7", "SM", 81, project_tp_root.data() + std::string("MARTY/src/MARTY/src/marty/models/sm.h"));
    
    return 0;
}