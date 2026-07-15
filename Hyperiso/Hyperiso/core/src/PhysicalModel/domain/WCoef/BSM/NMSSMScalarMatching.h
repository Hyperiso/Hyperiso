#ifndef NMSSM_SCALAR_MATCHING_H
#define NMSSM_SCALAR_MATCHING_H

#include "SourcesView.h"
#include "scalar.h"

namespace nmssm_scalar_matching {

struct Result {
    bool active {false};
    scalar_t cq1 {0.0};
    scalar_t cq2 {0.0};
};

bool is_active();
Result compute(const ParamSrc& src, int lepton_mass_slot);

} // namespace nmssm_scalar_matching

#endif
