#ifndef __WILSONBUILDER_H__
#define __WILSONBUILDER_H__

#include "IWilsonBuilder.h"
#include "IWilsonProvider.h"
#include "WilsonProvider.h"
#include "Include.h"
#include "WilsonManager.h"
#include "WilsonGroupFactory.h"

struct WilsonBuildConfig : public AbstractConfig {
    std::unordered_set<WGroup> groups;
    double matching_scale;
    double hadronic_scale;
    QCDOrder order {QCDOrder::LO}; 
};

class WilsonBuilder : public IWilsonBuilder<WilsonBuildConfig, WGroup> { 
public:
    WilsonBuilder(WilsonBuildConfig config);
    WilsonBuilder(std::shared_ptr<CoefficientManager> manager);

    void build(WilsonBuildConfig config) override;
    void add(WilsonBuildConfig config) override;
    void switch_basis(WGroup group_id) override;
    
    std::shared_ptr<IWilsonProvider> get_wilson_provider();

private:
    std::shared_ptr<CoefficientManager> cm;
};

#endif // __IWILSONBUILDER_H__
