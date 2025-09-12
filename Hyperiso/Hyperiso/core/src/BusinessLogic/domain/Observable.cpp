#include "Observable.h"

scalar_t Observable::get_exp_val() const {
    ObsParameterProxy opp {ParameterType::OBSERVABLE};
    return opp("FOBS", ObservableMapper::flha_of(ObservableMapper::to_id(this->id)).value(), DataType::VALUE);
}

scalar_t Observable::get_exp_uncertainty(UncertaintyType u_type) const {
    ObsParameterProxy opp {ParameterType::OBSERVABLE};
    return opp("FOBS", ObservableMapper::flha_of(ObservableMapper::to_id(this->id)).value(), UncertaintyTypeMapper::d_type(u_type));
}

Estimate Observable::get_exp() const {
    ObsParameterProxy opp {ParameterType::OBSERVABLE};
    Estimate est;
    est.central_value = opp("FOBS", ObservableMapper::flha_of(ObservableMapper::to_id(this->id)).value(), DataType::VALUE);
    est.stat_std = opp("FOBS", ObservableMapper::flha_of(ObservableMapper::to_id(this->id)).value(), DataType::STD_STAT);
    est.syst_std = opp("FOBS", ObservableMapper::flha_of(ObservableMapper::to_id(this->id)).value(), DataType::STD_SYST);

    return est;
}

scalar_t Observable::eval() const {
    decay_parent->enable();
    return decay_parent->compute_observable(id);
}