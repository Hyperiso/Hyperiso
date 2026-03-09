#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>

#include "minuit-cpp/FCNBase.hh"
#include "minuit-cpp/FunctionMinimum.hh"
#include "minuit-cpp/MnContours.hh"
#include "minuit-cpp/MnHesse.hh"
#include "minuit-cpp/MnMigrad.hh"
#include "minuit-cpp/MnUserParameters.hh"

namespace M2 = MinuitCpp;

// NLL gaussien : sum_i [ 0.5*((x_i-mu)/sigma)^2 + log(sigma) ] (+const)
class GaussianNLL final : public M2::FCNBase {
public:
  GaussianNLL(std::vector<double> data, double up) : data_(std::move(data)), up_(up) {}

  double operator()(const std::vector<double>& p) const override {
    const double mu = p.at(0);
    const double logSigma = p.at(1);
    const double sigma = std::exp(logSigma);
    if (!std::isfinite(sigma) || sigma <= 0.0) return 1e300;

    double nll = 0.0;
    for (double x : data_) {
      const double z = (x - mu) / sigma;
      nll += 0.5 * z * z + logSigma;
    }
    return nll;
  }

  // Up (= error definition) :
  // - pour -logL : 0.5 correspond aux erreurs 1D à 1σ
  // - pour contours 2D, on met Up = ΔNLL voulu (ex: 1.15, 2.995)
  double Up() const override { return up_; }

private:
  std::vector<double> data_;
  double up_;
};

static void write_bestfit_csv(const std::string& path,
                              double mu, double muErr,
                              double logSig, double logSigErr) {
  const double sig = std::exp(logSig);
  const double sigErr = std::fabs(sig * logSigErr); // propagation approx

  std::ofstream out(path);
  out << "name,value,error\n";
  out << "mu," << std::setprecision(17) << mu << "," << muErr << "\n";
  out << "sigma," << std::setprecision(17) << sig << "," << sigErr << "\n";
  out << "logSigma," << std::setprecision(17) << logSig << "," << logSigErr << "\n";
}

static void write_contours_csv(const std::string& path,
                               const std::vector<std::pair<double,double>>& c68,
                               const std::vector<std::pair<double,double>>& c95) {
  std::ofstream out(path);
  out << "cl,mu,logSigma,sigma\n";
  auto dump = [&](double cl, const std::vector<std::pair<double,double>>& c) {
    for (auto [mu, logSig] : c) {
      out << cl << "," << std::setprecision(17) << mu << "," << logSig << "," << std::exp(logSig) << "\n";
    }
  };
  dump(0.683, c68);
  dump(0.95,  c95);
}

int main() {
  // ---- Data synthétique ----
  constexpr double trueMu = 1.5;
  constexpr double trueSigma = 0.8;

  std::mt19937 rng(42);
  std::normal_distribution<double> gauss(trueMu, trueSigma);

  std::vector<double> data;
  data.reserve(400);
  for (int i = 0; i < 400; ++i) data.push_back(gauss(rng));

  // ---- Fit : Up=0.5 (NLL) ----
  GaussianNLL fcn_fit(data, /*up=*/0.5);

  M2::MnUserParameters upar;
  upar.Add("mu", 0.0, 0.1);
  upar.Add("logSigma", std::log(1.0), 0.05);

  // (Optionnel) limites pour éviter logSigma extrême
  upar.SetLimits("logSigma", std::log(1e-6), std::log(1e6));

  M2::MnMigrad migrad(fcn_fit, upar);
  M2::FunctionMinimum min = migrad(/*maxfcn=*/20000, /*tolerance=*/1e-8);

  // HESSE (améliore la covariance/erreurs)
  M2::MnHesse hesse;
  hesse(fcn_fit, min);

  if (!min.IsValid()) {
    std::cerr << "Minimisation invalide (min.IsValid()==false)\n";
    // std::cerr << min << "\n";
    return 2;
  }

  const auto& st = min.UserState();
  const double muHat = st.Value("mu");
  const double muErr = st.Error("mu");
  const double logSigmaHat = st.Value("logSigma");
  const double logSigmaErr = st.Error("logSigma");

  std::cout << "=== Fit result (MLE) ===\n";
  std::cout << "mu       = " << muHat << " +- " << muErr << "\n";
  std::cout << "sigma    = " << std::exp(logSigmaHat) << "  (logSigma=" << logSigmaHat
            << " +- " << logSigmaErr << ")\n";
  std::cout << "NLL(min) = " << min.Fval() << "\n";
  std::cout << "EDM      = " << min.Edm() << "\n";
  std::cout << "NFcn     = " << min.NFcn() << "\n";

  // ---- Contours (2 paramètres => indices 0 et 1) ----
  // Wilks 2D : Δχ²(68.3%)=2.30, Δχ²(95%)=5.99, donc ΔNLL = Δχ²/2 :contentReference[oaicite:2]{index=2}
  const double up68 = 2.30 / 2.0;   // 1.15
  const double up95 = 5.99 / 2.0;   // 2.995
  const unsigned px = 0;
  const unsigned py = 1;
  const unsigned npoints = 80;

  GaussianNLL fcn68(data, up68);
  GaussianNLL fcn95(data, up95);

  M2::MnContours contours68(fcn68, min);
  M2::MnContours contours95(fcn95, min);

  auto c68 = contours68(px, py, npoints);
  auto c95 = contours95(px, py, npoints);

  std::cout << "Contour sizes: 68%=" << c68.size() << "  95%=" << c95.size() << "\n";

  // ---- CSV ----
  write_bestfit_csv("bestfit.csv", muHat, muErr, logSigmaHat, logSigmaErr);
  write_contours_csv("contours.csv", c68, c95);
  std::cout << "Wrote bestfit.csv and contours.csv\n";

  return 0;
}