#ifndef __OBSWILSONHELPER_H__
#define __OBSWILSONHELPER_H__

#include "Include.h"
#include "WilsonFreezer.h"
#include "IObsWilsonBuilder.h"
#include "Configs.h"
#include "ObsWilsonProxy.h"

class ObsWilsonHelper {
public:
    static void build(WilsonBuildConfig config, std::shared_ptr<IObsWilsonBuilder<ObsWilsonProxy, WGroup>> wil_builder);

private:
    static std::unordered_set<WGroup> get_all_groups(const std::unordered_set<WGroup>& needed);
    static std::unordered_set<WGroup> update_state(const std::unordered_set<WGroup>& needed, std::shared_ptr<IObsWilsonBuilder<ObsWilsonProxy, WGroup>> wil_builder);
    static inline std::unordered_map<WGroup, bool> state;
};

#endif // __OBSWILSONHELPER_H__
