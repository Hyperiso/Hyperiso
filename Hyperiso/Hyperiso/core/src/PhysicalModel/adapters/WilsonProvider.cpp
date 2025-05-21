#include "WilsonProvider.h"
#include "WilsonBuilder.h" // Théo is reponsible for this horror

WilsonProvider::WilsonProvider(std::shared_ptr<CoefficientManager> manager) : cm(manager) {
    if (!cm) {
        LOG_ERROR("NullPointerError", "(WilsonProvider) CoefficientManager is null.");
    }
    
    if (cm->getGroups().empty()) {
        LOG_WARN("(WilsonProvider) CoefficientManager does not contain any coefficient groups. Pretty sus.");
    }
}

scalar_t WilsonProvider::get(std::shared_ptr<AbstractConfig> config) {
    WilsonRequest request = *static_cast<WilsonRequest*>(config.get());

    ScaleType scale_type = request.scale_type;
    bool sm_only = request.contribution == ContributionType::SM; // TODO : manage properly

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
                request.contribution) :
            cm->getRunCoefficient(
                GroupMapper::str(request.group),
                WCoefMapper::str(request.coefficient),
                OrderMapper::str(request.order),
                request.contribution);
    } else {
        LOG_ERROR("Invalid argument", "(WilsonProvider) Invalid scale type.");
    }
}

std::shared_ptr<WilsonBuilder> WilsonProvider::get_builder() {
    return std::make_shared<WilsonBuilder>(this->cm);
}
