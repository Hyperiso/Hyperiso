#ifndef OBS_PORTS_CONFIG_H
#define OBS_PORTS_CONFIG_H

#include "IObsParameterProxy.h"
#include "IObsWilsonBuilder.h"
#include "IObsCoreAPI.h"
#include "IWilsonFreezer.h"
#include "IObsQCDProxy.h"

struct ObservablePortsConfig {

    ObservablePortsConfig(std::shared_ptr<IObsWilsonBuilder> iobswb, std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm, std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_flav, std::shared_ptr<IObsQCDProxy> iobs_qcdp, std::shared_ptr<IObsCoreAPI<bool>> iobs_use_marty, std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer) :
        iobswb(iobswb), 
        iobspp_sm(iobspp_sm), iobspp_flav(iobspp_flav), iobs_qcdp(iobs_qcdp),
        iobs_use_marty(iobs_use_marty), iobs_wfreezer(iobs_wfreezer) {}


    std::shared_ptr<IObsWilsonBuilder> iobswb;

    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm;
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_flav;

    std::shared_ptr<IObsQCDProxy> iobs_qcdp;

    std::shared_ptr<IObsCoreAPI<bool>> iobs_use_marty;

    std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer;

};

#endif