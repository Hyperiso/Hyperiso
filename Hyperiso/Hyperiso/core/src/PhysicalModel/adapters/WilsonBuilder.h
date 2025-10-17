#ifndef __WILSONBUILDER_H__
#define __WILSONBUILDER_H__

#include "IWilsonBuilder.h"
#include "WilsonProvider.h"
#include "Include.h"
#include "WilsonManager.h"
#include "WilsonGroupFactory.h"
#include "Configs.h"
#include "WilsonParamComposer.h"
#include "ScaleSetter.h"
#include "Wilson_parameters.h"
#include "ModelAPI.h"
#include "UseMarty.h"
#include "susy_parameters.h"
#include "thdm_parameters.h"
#include "MartyWilsonProxy.h"

class WilsonBuilder : public IWilsonBuilder<WilsonBuildConfig, WilsonProvider>, public std::enable_shared_from_this<WilsonBuilder> { 
public:

    WilsonBuilder(WilsonBuildConfig config);
    WilsonBuilder(std::shared_ptr<CoefficientManager> manager);

    void build(WilsonBuildConfig config) override;
    void add(WilsonBuildConfig config) override;

    std::shared_ptr<WilsonProvider> get_wilson_provider();
    std::shared_ptr<CoefficientManager> get_coefficient_manager();

private:

    std::shared_ptr<CoefficientManager> cm;
};

#endif // __IWILSONBUILDER_H__
