#include "WilsonProvider.h"
#include "WilsonBuilder.h" // Théo is reponsible for this horror

WilsonProvider::WilsonProvider(std::shared_ptr<WilsonBuilder> builder) : builder(builder), cm(builder->get_coefficient_manager()) {
    if (!cm) {
        LOG_ERROR("NullPointerError", "(WilsonProvider) CoefficientManager is null.");
    }
    
    if (cm->getGroups().empty()) {
        LOG_WARN("(WilsonProvider) CoefficientManager does not contain any coefficient groups. Pretty sus.");
    }
}

scalar_t WilsonProvider::get(std::shared_ptr<AbstractConfig> config) {
    WilsonRequest request = *static_cast<WilsonRequest*>(config.get());

    // if (UseMarty().get() && request.order > QCDOrder::LO) {
    //     LOG_WARN("Trying to access coefficient at ", OrderMapper::str(request.order), "but coefficients were only calculated at LO (MARTY). Defaulting to LO.");
    //     request.order = QCDOrder::LO;
    // }

    ScaleType scale_type = request.scale_type;

    if (scale_type == ScaleType::MATCHING) {
        return request.sum_qcd_orders ? 
            cm->getFullMatchingCoefficient(
                GroupMapper::str(request.group),
                WCoefMapper::str(request.coefficient),
                OrderMapper::str(request.order),
                request.contribution) :
            cm->getMatchingCoefficient(
                GroupMapper::str(request.group),
                WCoefMapper::str(request.coefficient),
                OrderMapper::str(request.order),
                request.contribution);
    } else if (scale_type == ScaleType::HADRONIC) {
        return request.sum_qcd_orders ? 
            cm->getFullRunCoefficient(
                GroupMapper::str(request.group),
                WCoefMapper::str(request.coefficient),
                OrderMapper::str(request.order),
                request.contribution,
                request.basis.value_or(WilsonBasis::B_STANDARD)) :
            cm->getRunCoefficient(
                GroupMapper::str(request.group),
                WCoefMapper::str(request.coefficient),
                OrderMapper::str(request.order),
                request.contribution,
                request.basis.value_or(WilsonBasis::B_STANDARD));
    } else {
        LOG_ERROR("Invalid argument", "(WilsonProvider) Invalid scale type.");
    }
}

std::shared_ptr<WilsonBuilder> WilsonProvider::get_builder() {
    return builder;
}

std::unordered_set<WilsonBasis> WilsonProvider::get_bases(WGroupId group) {
    return this->cm->getGroupBases(group);
}
