#include "Wilson_THDMv2.h"

std::complex<double> C4_THDM::NLO_calculation() {
    double coeff_temp = EH(thdm_params->yt,thdm_params->lu);
    return this->double_to_complex_save("NLO", coeff_temp);
}

std::complex<double> C7_THDM::LO_calculation() {
    double coeff_temp = 1./3.*thdm_params->lu*thdm_params->lu*F7_1(thdm_params->yt) - thdm_params->lu*thdm_params->ld*F7_2(thdm_params->yt);
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> C7_THDM::NLO_calculation() {
    double coeff_temp = G7H(thdm_params->yt,thdm_params->lu,thdm_params->ld)+Delta7H(thdm_params->yt,thdm_params->lu,thdm_params->ld)*log(pow(this->get_Q_match()/thdm_params->m_H,2.))
    -4./9.*EH(thdm_params->yt,thdm_params->lu);
    return this->double_to_complex_save("NLO", coeff_temp);
}

std::complex<double> C8_THDM::LO_calculation() {
    double coeff_temp = 1./3.*thdm_params->lu*thdm_params->lu*F8_1(thdm_params->yt) - thdm_params->lu*thdm_params->ld*F8_2(thdm_params->yt);
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> C8_THDM::NLO_calculation() {
    double coeff_temp = G8H(thdm_params->yt,thdm_params->lu,thdm_params->ld)+Delta8H(thdm_params->yt,thdm_params->lu,thdm_params->ld)*log(pow(this->get_Q_match()/thdm_params->m_H,2.))
    -1./6.*EH(thdm_params->yt,thdm_params->lu);
    return this->double_to_complex_save("NLO", coeff_temp);
}

std::complex<double> C9_THDM::LO_calculation() {
    double coeff_temp = (1.-4.*thdm_params->sw2)/thdm_params->sw2*C9llH0(thdm_params->xt,thdm_params->yt,thdm_params->lu)-D9H0(thdm_params->yt,thdm_params->lu);
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> C9_THDM::NLO_calculation() {
    double coeff_temp = (1.-4.*thdm_params->sw2)/thdm_params->sw2*C9llH1(thdm_params->xt,thdm_params->yt,thdm_params->lu,log(pow(this->get_Q_match()/thdm_params->m_H,2.)))
    -D9H1(thdm_params->yt,thdm_params->lu,log(pow(this->get_Q_match()/thdm_params->m_H,2.)));
    return this->double_to_complex_save("NLO", coeff_temp);
}

std::complex<double> C10_THDM::LO_calculation() {
    double coeff_temp = -C9llH0(thdm_params->xt,thdm_params->yt,thdm_params->lu)/thdm_params->sw2;
    return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> C10_THDM::NLO_calculation() {
    double coeff_temp = -C9llH1(thdm_params->xt,thdm_params->yt,thdm_params->lu,log(pow(this->get_Q_match()/thdm_params->m_H,2.)))/thdm_params->sw2;
    return this->double_to_complex_save("NLO", coeff_temp);
}