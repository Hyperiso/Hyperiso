#pragma once
#include <vector>
#include <cmath>
#include "ports/IModel.h"


// Simple model mirroring your Python formula for quick testing
class BToMuMuToy final : public IModel {
public:
BToMuMuToy() = default;


std::size_t n_observables() const override { return 2; }


Vec predict(const Vec& p, const Vec& eta) override {
// p: [C10, Cp10]
const double C10 = p.at(0);
const double Cp10 = p.at(1);
// eta: [f_Bd, f_Bs, y_s, |Vts|, |Vtd|]
const double f_d = eta.at(0);
const double f_s = eta.at(1);
const double ys = eta.at(2);
const double Vts = eta.at(3);
const double Vtd = eta.at(4);


// constants
const double HBAR = 6.58211889e-25;
const double G_F = 1.1663787e-5;
const double alpha_em = 1.0/1.27930e+02;
const double m_mu = 0.1056583755;
const double m_Bd = 5.27972;
const double m_Bs = 5.36693;
const double tau_Bd = 1.517e-12;
const double tau_Bs = 1.527e-12;
const double abs_Vtb = 0.999118;


const double k = (G_F*G_F)*(alpha_em*alpha_em)/(16.0*M_PI*M_PI*M_PI);
const double ps_s = std::sqrt(1.0 - 4.0*m_mu*m_mu/(m_Bs*m_Bs));
const double ps_d = std::sqrt(1.0 - 4.0*m_mu*m_mu/(m_Bd*m_Bd));


const double BR_s_untag = (1.0 + ys)/(1.0 - ys*ys) * k * f_s*f_s * tau_Bs * m_Bs * m_mu*m_mu * Vts*Vts * abs_Vtb*abs_Vtb * ps_s * (C10 - Cp10)*(C10 - Cp10) / HBAR;
const double BR_d = k * f_d*f_d * tau_Bd * m_Bd * m_mu*m_mu * Vtd*Vtd * abs_Vtb*abs_Vtb * ps_d * C10*C10 / HBAR;
return Vec{BR_s_untag, BR_d};
}

std::map<ObservableId, double> predict(const std::map<ParamId, double>& p, const std::map<ParamId, double>& eta) override {
// p: [C10, Cp10]
// const double C10 = p.at(0);
// const double Cp10 = p.at(1);
// // eta: [f_Bd, f_Bs, y_s, |Vts|, |Vtd|]
// const double f_d = eta.at(0);
// const double f_s = eta.at(1);
// const double ys = eta.at(2);
// const double Vts = eta.at(3);
// const double Vtd = eta.at(4);


// // constants
// const double HBAR = 6.58211889e-25;
// const double G_F = 1.1663787e-5;
// const double alpha_em = 1.0/1.27930e+02;
// const double m_mu = 0.1056583755;
// const double m_Bd = 5.27972;
// const double m_Bs = 5.36693;
// const double tau_Bd = 1.517e-12;
// const double tau_Bs = 1.527e-12;
// const double abs_Vtb = 0.999118;


// const double k = (G_F*G_F)*(alpha_em*alpha_em)/(16.0*M_PI*M_PI*M_PI);
// const double ps_s = std::sqrt(1.0 - 4.0*m_mu*m_mu/(m_Bs*m_Bs));
// const double ps_d = std::sqrt(1.0 - 4.0*m_mu*m_mu/(m_Bd*m_Bd));


// const double BR_s_untag = (1.0 + ys)/(1.0 - ys*ys) * k * f_s*f_s * tau_Bs * m_Bs * m_mu*m_mu * Vts*Vts * abs_Vtb*abs_Vtb * ps_s * (C10 - Cp10)*(C10 - Cp10) / HBAR;
// const double BR_d = k * f_d*f_d * tau_Bd * m_Bd * m_mu*m_mu * Vtd*Vtd * abs_Vtb*abs_Vtb * ps_d * C10*C10 / HBAR;
// return Vec{BR_s_untag, BR_d};
return std::map<ObservableId, double>{};
}

void add_observables(std::map<ObservableId, QCDOrder> obs_ids) override {};
std::unordered_set<ParamId> get_obs_deps(ObservableId id) override {return std::unordered_set<ParamId>{};}
};