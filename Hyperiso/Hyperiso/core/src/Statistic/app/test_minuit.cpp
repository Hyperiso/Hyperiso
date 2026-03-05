#include <cmath>
#include <iostream>
#include <random>
#include <vector>

  #include "minuit-cpp/FCNBase.hh"
  #include "minuit-cpp/FunctionMinimum.hh"
  #include "minuit-cpp/MnMigrad.hh"
  #include "minuit-cpp/MnUserParameters.hh"
  namespace M2 = MinuitCpp;
// #ifdef USE_MINUITCPP
//   #include "minuit-cpp/FCNBase.hh"
//   #include "minuit-cpp/FunctionMinimum.hh"
//   #include "minuit-cpp/MnMigrad.hh"
//   #include "minuit-cpp/MnUserParameters.hh"
//   namespace M2 = MinuitCpp;
// #else
//   #include "Minuit2/FCNBase.h"
//   #include "Minuit2/FunctionMinimum.h"
//   #include "Minuit2/MnMigrad.h"
//   #include "Minuit2/MnUserParameters.h"
//   namespace M2 = ROOT::Minuit2;
// #endif

// NLL gaussien : somme_i [ 0.5*((x_i-mu)/sigma)^2 + log(sigma) ] + constante
// On paramètre sigma via logSigma pour garantir sigma>0 sans limites.
class GaussianNLL final : public M2::FCNBase {
public:
  explicit GaussianNLL(std::vector<double> data) : data_(std::move(data)) {}

  double operator()(const std::vector<double>& p) const override {
    const double mu = p.at(0);
    const double logSigma = p.at(1);
    const double sigma = std::exp(logSigma);

    // garde-fou (normalement inutile avec exp, mais évite les NaN si logSigma très négatif)
    if (!std::isfinite(sigma) || sigma <= 0.0) return 1e300;

    double nll = 0.0;
    for (double x : data_) {
      const double z = (x - mu) / sigma;
      nll += 0.5 * z * z + logSigma; // + 0.5*log(2*pi) est une constante -> optionnelle
    }
    return nll;
  }

  // Convention Minuit : Up=1 pour chi2, Up=0.5 pour -log(L)
  double Up() const override { return 0.5; }

private:
  std::vector<double> data_;
};

int main() {
  // --- Génère un petit dataset gaussien (juste pour le test) ---
  constexpr double trueMu = 1.5;
  constexpr double trueSigma = 0.8;

  std::mt19937 rng(42);
  std::normal_distribution<double> gauss(trueMu, trueSigma);

  std::vector<double> data;
  data.reserve(200);
  for (int i = 0; i < 200; ++i) data.push_back(gauss(rng));

  // --- FCN (likelihood) ---
  GaussianNLL fcn(data);

  // --- Paramètres (valeurs initiales + pas) ---
  M2::MnUserParameters upar;
  upar.Add("mu", 0.0, 0.1);
  upar.Add("logSigma", std::log(1.0), 0.05);

  // --- Minimisation (MIGRAD) ---
  M2::MnMigrad migrad(fcn, upar);
  M2::FunctionMinimum min = migrad();

//   std::cout << "=== FunctionMinimum ===\n" << min << "\n";

  if (!min.IsValid()) {
    std::cerr << "Minimisation invalide (min.IsValid()==false)\n";
    return 2;
  }

  const auto& st = min.UserState();
  const double muHat = st.Value("mu");
  const double muErr = st.Error("mu");

  const double logSigmaHat = st.Value("logSigma");
  const double logSigmaErr = st.Error("logSigma");
  const double sigmaHat = std::exp(logSigmaHat);

  std::cout << "=== Fit result (MLE) ===\n";
  std::cout << "mu       = " << muHat << " +- " << muErr << "\n";
  std::cout << "sigma    = " << sigmaHat << "   (via logSigma=" << logSigmaHat
            << " +- " << logSigmaErr << ")\n";
  std::cout << "NLL(min) = " << min.Fval() << "\n";
  std::cout << "EDM      = " << min.Edm() << "\n";
  std::cout << "NFcn     = " << min.NFcn() << "\n";

  return 0;
}