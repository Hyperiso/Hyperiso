#include "QEDProvider.h"

double QEDProvider::operator()(AlphasConfig config) {
    return EWHelper::alpha_em(config.scale);
}
