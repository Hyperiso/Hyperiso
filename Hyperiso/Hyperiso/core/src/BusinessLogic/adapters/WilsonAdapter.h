#ifndef WILSON_ADAPTER_H
#define WILSON_ADAPTER_H

#include "WilsonInterface.h"

class WilsonAdapter {
public:
    WilsonAdapter() = default;

    void build(std::vector<WGroup> groups, double mu_W, double mu_h, QCDOrder order) {wi.build(groups, mu_W, mu_h, order); built = true;}

    void addWilsonGroup(WGroup group_name) {wi.addWilsonGroup(group_name);}

    void switchbasis(WGroup group_name) {wi.switchbasis(group_name);}

private:
    static inline WilsonInterface wi {};
    static inline bool built {false};

friend class ObsWilsonProxy;
};

#endif