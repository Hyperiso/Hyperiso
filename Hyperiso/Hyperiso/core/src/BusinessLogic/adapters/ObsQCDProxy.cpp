#include "ObsQCDProxy.h"

double ObsQCDProxy::operator()(AlphasConfig config) {
    return QCDHelper::alpha_s(config.scale, config.m_b_type, config.m_t_type);
}

double ObsQCDProxy::operator()(MassConfig config) {
    return QCDHelper::msbar_mass(config.pdg_id, config.scale, config.m_b_type, config.m_t_type);
}

QCDConstants *ObsQCDProxy::get_constants() {
    return QCDHelper::constants;
}
