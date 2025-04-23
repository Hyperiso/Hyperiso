#include "WilsonInterface.h"

class WilsonAdapter {
    WilsonAdapter() = default;

    void build(std::vector<WGroup> groups, double mu_W, double mu_h, QCDOrder order) {wi.build(groups, mu_W, mu_h, order);}

    void addWilsonGroup(WGroup group_name) {wi.addWilsonGroup(group_name);}

    void switchbasis(WGroup group_name) {wi.switchbasis(group_name);}
    
private:
    WilsonInterface wi = WilsonInterface();
};