#ifndef __OBSWILSONHELPER_H__
#define __OBSWILSONHELPER_H__

#include "Include.h"
#include "WilsonFreezer.h"
#include "ObsWilsonBuilder.h"
#include "Configs.h"
#include "ObsWilsonProxy.h"

class ObsWilsonHelper {
public:
    static void build(WilsonBuildConfig config, std::shared_ptr<ObsWilsonBuilder>& wil_builder);
    ObsWilsonHelper(bool reset) {reset ? state = {} : state;}
    ObsWilsonHelper() = default;
private:
    static std::unordered_set<WGroupId> get_all_groups(const std::unordered_set<WGroupId>& needed);
    static std::unordered_set<WGroupId> update_state(const std::unordered_set<WGroupId>& needed, std::shared_ptr<ObsWilsonBuilder>& wil_builder);
    static inline std::unordered_map<WGroupId, bool> state;
};

#endif // __OBSWILSONHELPER_H__
