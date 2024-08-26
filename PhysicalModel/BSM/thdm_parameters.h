#include <algorithm>
#include <array>
#include <functional>
#include "Parameters.h"
#include "Logger.h"

class thdm_parameters {

    static thdm_parameters* instance;
    double scale{81};
	bool is_PrimeCQG = false;

    explicit thdm_parameters();
    thdm_parameters(const thdm_parameters&) = delete;
    void operator=(const thdm_parameters&) = delete;


	Parameters* mod = Parameters::GetInstance(2);
    Parameters* sm = Parameters::GetInstance();

public:
    static thdm_parameters* GetInstance() {
        if (!thdm_parameters::instance) {
            thdm_parameters::instance = new thdm_parameters();
        }
        return thdm_parameters::instance;
    }
    void set_lu(double lu) {this->lu = lu;}
    void set_ld(double ld) {this->ld = ld;}
    void set_le(double le) {this->le = le;}

    void set_sm_parameters(Parameters* sm) {this->sm = sm;}
    void set_mod_parameters(Parameters* sm) {this->mod = mod;}
    void set_params(double Q_match);
    void set_gen(int new_gen) {this->gen = new_gen; le = (*mod)("YL",10*(gen-1)+gen-1); ml = (*sm)("MASS", 13+2*(this->gen-2));}

    double mass_top_muW=(*sm).running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale, "running", "pole");
	double mass_b_muW=(*sm).running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale);

    double sw2=pow(sin(atan((*sm)("GAUGE",1)/(*sm)("GAUGE",2))),2.); //1 = param-> gp and 2 = param->g2
    double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*mod)("MASS",37),2.); // param->mass_H (25)
    double xh=pow((*mod)("MASS",25)/(*sm)("MASS",24),2.);

	double L; // scale -> mu_W
	double alphas_muW;
	double z;
    double m_H = (*mod)("MASS", 37);
    int gen{3};
    double ml = (*sm)("MASS", 13+2*(this->gen-2));
    double lu = (*mod)("YU", 22);;
	double ld = (*mod)("YD", 22);
    double le = (*mod)("YL",10*(gen-1)+gen-1);
};