#ifndef WILSON_ADAPTER_H
#define WILSON_ADAPTER_H

#include "WilsonInterface.h"

class WilsonAdapter {
public:
    WilsonAdapter();

    void build(WilsonConfig config) {wi->build(config); built = true;}

    void addWilsonGroup(WGroup group_name) {wi->addWilsonGroup(group_name);}

    void switchbasis(WGroup group_name) {wi->switchbasis(group_name);}

private:
    static inline std::shared_ptr<WilsonInterface> wi;
    static inline bool built {false};

friend class ObsWilsonProxy;
};

#endif