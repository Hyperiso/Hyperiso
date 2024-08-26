#include "thdm_parameters.h"

void thdm_parameters::set_params(double Q_match) {
    this->scale = Q_match;

    mass_top_muW=(*sm).running_mass((*sm)("MASS",6), (*sm)("MASS",6),Q_match, "running", "pole");
	mass_b_muW=(*sm).running_mass((*sm)("MASS",5), (*sm)("MASS",5), Q_match);

    sw2=pow(sin(atan((*sm)("GAUGE",1)/(*sm)("GAUGE",2))),2.); //1 = param-> gp and 2 = param->g2
    xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	yt= pow(mass_top_muW/(*mod)("MASS",37),2.); // param->mass_H (25)
    xh=pow((*mod)("MASS",25)/(*sm)("MASS",24),2.);

    m_H = (*mod)("MASS", 37);
}