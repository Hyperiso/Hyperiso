#include "Wilson_THDMv2.h"

std::complex<double> C3_THDM::NNLO_calculation() {
    double coeff_temp = G3H(thdm_params->yt,thdm_params->lu)+Delta3H(thdm_params->yt,thdm_params->lu)*log(pow(this->get_Q_match()/(*mod)("MASS",37),2.));
    return this->double_to_complex_save("NNLO", coeff_temp);
}

std::complex<double> C4_THDM::NLO_calculation() {
    double coeff_temp = EH(thdm_params->yt,thdm_params->lu);
    return this->double_to_complex_save("NLO", coeff_temp);
}

std::complex<double> C4_THDM::NNLO_calculation() {
    double coeff_temp=G4H(thdm_params->yt,thdm_params->lu)+Delta4H(thdm_params->yt,thdm_params->lu)*log(pow(this->get_Q_match()/(*mod)("MASS",37),2.));
    return this->double_to_complex_save("NNLO", coeff_temp);
}

std::complex<double> C5_THDM::NNLO_calculation() {
    double C4H_1=EH(thdm_params->yt,thdm_params->lu);
    double C3H_2=G3H(thdm_params->yt,thdm_params->lu)+Delta3H(thdm_params->yt,thdm_params->lu)*log(pow(this->get_Q_match()/(*mod)("MASS",37),2.));
	double coeff_temp=-C3H_2/10.+2./15.*C4H_1;
    return this->double_to_complex_save("NNLO", coeff_temp);
}

std::complex<double> C6_THDM::NNLO_calculation() {
    double C4H_1=EH(thdm_params->yt,thdm_params->lu);
    double C3H_2=G3H(thdm_params->yt,thdm_params->lu)+Delta3H(thdm_params->yt,thdm_params->lu)*log(pow(this->get_Q_match()/(*mod)("MASS",37),2.));
	double coeff_temp=-3./16.*C3H_2+1./4.*C4H_1;

    return this->double_to_complex_save("NNLO", coeff_temp);
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

std::complex<double> C7_THDM::NNLO_calculation() {
	double coeff_temp =C7H2(thdm_params->yt,thdm_params->lu,thdm_params->ld,log(pow(this->get_Q_match()/thdm_params->mass_top_muW, 2.)));
    return this->double_to_complex_save("NNLO", coeff_temp);
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

std::complex<double> C8_THDM::NNLO_calculation() {
	double coeff_temp =C8H2(thdm_params->yt,thdm_params->lu,thdm_params->ld,log(pow(this->get_Q_match()/thdm_params->mass_top_muW, 2.)));
    return this->double_to_complex_save("NNLO", coeff_temp);
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

std::complex<double> CQ1_THDM::LO_calculation() {
    double le = (*mod)("YL",10*(gen-1)+gen-1);
	double G1=-3./4.+thdm_params->ld*thdm_params->lu*F4SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*thdm_params->lu*F5SP(thdm_params->xt,thdm_params->xH);
	double G2=thdm_params->ld*(thdm_params->ld*thdm_params->lu+1.)*F6SP(thdm_params->xt,thdm_params->xH)-thdm_params->ld*thdm_params->lu*thdm_params->lu*F7SP(thdm_params->xt,thdm_params->xH)
	+thdm_params->lu*thdm_params->lu*(thdm_params->ld*F8SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F9SP(thdm_params->xt,thdm_params->xH)-thdm_params->lu*F10SP(thdm_params->xt,thdm_params->xH))+thdm_params->lu*F11SP(thdm_params->xt,thdm_params->xH)-thdm_params->lu*F12SP(thdm_params->xt,thdm_params->xH);

	double CSn_2HDM=thdm_params->xt*(F0SP(thdm_params->xt)+le*(thdm_params->ld*F1SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F2SP(thdm_params->xt,thdm_params->xH))+le*thdm_params->lu*F3SP(thdm_params->xt,thdm_params->xH))
	+thdm_params->xt/2./thdm_params->xh*(sin(thdm_params->alpha-thdm_params->beta)+cos(thdm_params->alpha-thdm_params->beta)*le)*(sin(thdm_params->alpha-thdm_params->beta)*G1+cos(thdm_params->alpha-thdm_params->beta)*G2)
	+thdm_params->xt/2./thdm_params->xH0*(cos(thdm_params->alpha-thdm_params->beta)-sin(thdm_params->alpha-thdm_params->beta)*le)*(cos(thdm_params->alpha-thdm_params->beta)*G1-sin(thdm_params->alpha-thdm_params->beta)*G2);
	double coeff_temp=CSc_2HDM(thdm_params->xH,thdm_params->xt,thdm_params->lu,thdm_params->ld,le)+CSn_2HDM;
	coeff_temp*=(W_param->ml*thdm_params->mass_b_muW/(*sm)("MASS",24)/(*sm)("MASS",24))/thdm_params->sw2;
	return this->double_to_complex_save("LO", coeff_temp);
}

std::complex<double> CQ2_THDM::LO_calculation() {
    double le = (*mod)("YL",10*(gen-1)+gen-1);
    double G3=thdm_params->ld*(thdm_params->ld*thdm_params->lu+1.)*F6SP(thdm_params->xt,thdm_params->xH)+thdm_params->ld*thdm_params->lu*thdm_params->lu*F7SP(thdm_params->xt,thdm_params->xH)
	+thdm_params->lu*thdm_params->lu*(thdm_params->ld*F8SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F9SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F10SP(thdm_params->xt,thdm_params->xH))+thdm_params->lu*F11SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F12SP(thdm_params->xt,thdm_params->xH);
    double CPn_2HDM=thdm_params->xt*(-le*(thdm_params->ld*F1SP(thdm_params->xt,thdm_params->xH)+thdm_params->lu*F2SP(thdm_params->xt,thdm_params->xH))+le*thdm_params->lu*F3SP(thdm_params->xt,thdm_params->xH))+thdm_params->xt/2./thdm_params->xA*(le)*G3;

    double coeff_temp=CPc_2HDM(thdm_params->xH,thdm_params->xt,thdm_params->lu,thdm_params->ld,le,thdm_params->sw2)+CPn_2HDM;
    coeff_temp*=(W_param->ml*thdm_params->mass_b_muW/(*sm)("MASS",24)/(*sm)("MASS",24))/thdm_params->sw2;
    return this->double_to_complex_save("LO", coeff_temp);
}