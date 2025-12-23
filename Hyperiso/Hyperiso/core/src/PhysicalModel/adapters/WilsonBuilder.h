#ifndef WILSONBUILDER_H
#define WILSONBUILDER_H

#include "IWilsonBuilder.h"
#include "WilsonProvider.h"
#include "Include.h"
#include "WilsonManager.h"
// #include "WilsonGroupFactory.h"
#include "Configs.h"
#include "WilsonParamComposer.h"
#include "ScaleSetter.h"
#include "WilsonParametersHelper.h"
#include "ModelAPI.h"
#include "UseMarty.h"
#include "HasWilsonAPI.h"
#include "SUSYParametersHelper.h"
#include "THDMParametersHelper.h"
#include "MartyWilsonProxy.h"
#include "MartyModelNameAPI.h"
#include "MartyModelPathAPI.h"

class WilsonBuilder : public IWilsonBuilder<WilsonBuildConfig, WilsonProvider>, public std::enable_shared_from_this<WilsonBuilder> { 
public:

    WilsonBuilder(WilsonBuildConfig config);
    WilsonBuilder(std::shared_ptr<CoefficientManager> manager);

    void build(WilsonBuildConfig config) override;
    void add(WilsonBuildConfig config) override;

    std::shared_ptr<WilsonProvider> get_wilson_provider();
    std::shared_ptr<CoefficientManager> get_coefficient_manager();

private:
    std::map<Model, std::shared_ptr<IWilsonParameterHelper>> wilson_param_helpers;
    std::shared_ptr<CoefficientManager> cm;
};

#endif // IWILSONBUILDER_H
