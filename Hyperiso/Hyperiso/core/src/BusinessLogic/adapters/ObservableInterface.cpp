#include "ObservableInterface.h"

ObservableInterface::ObservableInterface() {
    WilsonBuildConfig config;
    std::shared_ptr<WilsonBuilder> builder_ptr = std::make_shared<WilsonBuilder>(config);
    std::shared_ptr<IObsWilsonBuilder> builder = std::make_shared<ObsWilsonBuilder>(builder_ptr);
    

    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm = std::make_shared<ObsParameterProxy>();

    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_flav = std::make_shared<ObsParameterProxy>(ParameterType::FLAVOR);
    std::shared_ptr<IObsQCDProxy> iobs_qcdp = std::make_shared<ObsQCDProxy>();


    std::shared_ptr<IObsCoreAPI<bool>> iobs_use_marty = std::make_shared<ObsUseMarty>();

    std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer = std::make_shared<WilsonFreezer>(builder);
    
    ports = std::make_shared<ObservablePortsConfig>(builder, iobspp_sm, iobspp_flav, iobs_qcdp, iobs_use_marty, iobs_wfreezer);

    manager = std::make_shared<ObsManager>(*ports);
}

void ObservableInterface::add_custom_decay(DecayId id, std::shared_ptr<DecayParent> ptr) {
    manager->add_custom_decay(id, ptr);
}