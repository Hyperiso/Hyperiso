#include "Wilson.h"


bool WilsonCoefficient::fill_from_flha() {
    if (!from_lha && HAS_WILSON_API().get()) {
        auto wc = Parameters::GetInstance(ParameterType::WILSON);
        WCoef id = WCoefMapper::enum_elt(this->get_name());
        Model m = ModelAPI().get();
        int w_type = (*wc)("REWCOEF", -2);
        if (m != Model::SM && w_type == 0) {
            LOG_ERROR("Value", "SM Wilsons coefficients were given, but the selected model is not SM.");
        }

        this->set_CoefficientMatchingValue("LO", complex_t((*wc)("REWCOEF", (int)id * 10 + (int)QCDOrder::LO-1), 
                                                            (*wc)("IMWCOEF", (int)id * 10 + (int)QCDOrder::LO-1)));
        this->set_CoefficientMatchingValue("NLO", complex_t((*wc)("REWCOEF", (int)id * 10 + (int)QCDOrder::NLO-1), 
                                                             (*wc)("IMWCOEF", (int)id * 10 + (int)QCDOrder::NLO-1)));
        this->set_CoefficientMatchingValue("NNLO", complex_t((*wc)("REWCOEF", (int)id * 10 + (int)QCDOrder::NNLO-1), 
                                                              (*wc)("IMWCOEF", (int)id * 10 + (int)QCDOrder::NNLO-1)));         
        from_lha = true;
    }
    return from_lha;
}

void WilsonCoefficient::set_owned(bool owned) {
    if (this->is_owned && owned) {
        LOG_ERROR("LogicError", "WilsonCoefficient is already owned by a WilsonGroup and cannot be shared.");
    }

    this->is_owned = owned;
}
