#ifndef __WILSONPROVIDER_H__
#define __WILSONPROVIDER_H__

#include "IWilsonProvider.h"
#include "Include.h"
#include "WilsonManager.h"
#include "Configs.h"

class WilsonBuilder;

class WilsonProvider : public IWilsonProvider<WilsonBuilder> {
public:
    WilsonProvider(std::shared_ptr<WilsonBuilder> builder);

    scalar_t get(std::shared_ptr<AbstractConfig> config) override;
    std::shared_ptr<WilsonBuilder> get_builder() override;
    std::unordered_set<WilsonBasis> get_bases(WGroupId group) override;

private:
    std::shared_ptr<WilsonBuilder> builder;
    std::shared_ptr<CoefficientManager> cm;
};

#endif // __WILSONPROVIDER_H__
