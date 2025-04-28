#include "WilsonAdapter.h"


WilsonAdapter::WilsonAdapter() {wi = std::make_shared<WilsonInterface>();}

void WilsonAdapter::build(WilsonConfig config) {
    wi->build(config);
    built = true;
}
