#include "QCDProvider.h"

double QCDProvider::operator()(AlphasConfig config) {
    return QCDHelper::alpha_s(config.scale, config.m_b_type, config.m_t_type);
}

double QCDProvider::operator()(MassConfig config) {
    return QCDHelper::msbar_mass(config.pdg_id, config.scale, config.m_b_type, config.m_t_type);
}

double QCDProvider::operator()(QCDMass mass_id) {
    switch (mass_id) {
    case QCDMass::C_POLE:
        return QCDHelper::mass_c_pole();
        break;
    case QCDMass::B_POLE:
        return QCDHelper::mass_b_pole();
        break;
    case QCDMass::B_MSBAR:
        return QCDHelper::mass_b_msbar();
        break;
    case QCDMass::B_1S:
        return QCDHelper::mass_b_1S();
        break;
    case QCDMass::T_POLE:
        return QCDHelper::mass_t_pole();
        break;
    case QCDMass::T_MSBAR:
        return QCDHelper::mass_t_msbar();
        break;
    default:
        break;
    }

    LOG_ERROR("ValueError", "Unknown special QCD mass.");
}

QCDConstants *QCDProvider::get_constants() {
    return QCDHelper::constants;
}
