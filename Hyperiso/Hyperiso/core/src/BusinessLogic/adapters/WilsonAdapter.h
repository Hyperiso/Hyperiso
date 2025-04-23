#include "WilsonInterface.h"

class WilsonAdapter {
    WilsonAdapter() = default;

    void build(std::vector<WGroup> groups, double mu_W, double mu_h, QCDOrder order) {wi.build(groups, mu_W, mu_h, order);}
private:
    WilsonInterface wi = WilsonInterface();
};