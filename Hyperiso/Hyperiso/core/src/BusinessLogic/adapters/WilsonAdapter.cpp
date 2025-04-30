#include "WilsonAdapter.h"


WilsonAdapter::WilsonAdapter() {wi = std::make_shared<WilsonInterface>();}

void WilsonAdapter::build(WilsonBuildConfig config) {
    wi->build(config);
    built = true;
}
