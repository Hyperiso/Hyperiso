// #include <algorithm>
// #include <cassert>
// #include <chrono>
// #include <cmath>
// #include <functional>
// #include <iomanip>
// #include <iostream>
// #include <memory>
// #include <stdexcept>
// #include <string>
// #include <vector>

// #include "StatisticManager.h"
// #include "ObservableInterfaceProxy.h"
// #include "StatCorrelationProxy.h"
// #include "StatParameterProxy.h"
// #include "ObservableInterface.h"
// #include "StatParamSourcesProxy.h"
// #include "StatDependencyPruner.h"
// #include "BaseLikelihood.h"
// #include "DefaultNuisancePathsProvider.h"

// #include "minuit-cpp/FCNBase.hh"
// #include "minuit-cpp/FunctionMinimum.hh"
// #include "minuit-cpp/MnHesse.hh"
// #include "minuit-cpp/MnMigrad.hh"
// #include "minuit-cpp/MnUserCovariance.hh"
// #include "minuit-cpp/MnUserParameters.hh"
// #include "minuit-cpp/MnUserParameterState.hh"

// #include "BlockProxy.h"
// #include "NuisanceReader.h"

// namespace M2 = MinuitCpp;

// static void print_vec(const std::vector<double>& vec) {
//     std::cout << "[ ";
//     for (size_t i = 0; i < vec.size(); i++) {
//         std::cout << vec.at(i) << (i == vec.size() - 1 ? " " : ", ");
//     }
//     std::cout << "]\n";
// }



// static double safe_step(double scale) {
//     double s = std::fabs(scale);
//     if (!std::isfinite(s) || s <= 0.0) return 1.0;
//     return 0.1 * s;
// }

// struct MinuitJointFit {
//     Vector x_hat;                 // [p, eta]
//     std::vector<double> x_err; 
//     RealMatrix cov;          
//     double fmin = 0.0;
//     bool ok = false;
// };

// class GenericFCN final : public M2::FCNBase {
// public:
//     explicit GenericFCN(std::function<double(const std::vector<double>&)> f) : f_(std::move(f)) {}

//     double operator()(const std::vector<double>& x) const override {
//         const double v = f_(x);
//         return std::isfinite(v) ? v : 1e300;
//     }

//     // NLL => Up=0.5
//     double Up() const override { return 0.5; }

// private:
//     std::function<double(const std::vector<double>&)> f_;
// };

// static MinuitJointFit minuit_migrad_hesse(
//     const std::function<double(const std::vector<double>&)>& f,
//     const std::vector<double>& x0,
//     const std::vector<double>& scales,
//     std::size_t max_fcn,
//     double tol_edm
// ) {
//     if (x0.size() != scales.size()) {
//         throw std::invalid_argument(
//             "x0 and scales must have same size (x0=" + std::to_string(x0.size()) +
//             ", scales=" + std::to_string(scales.size()) + ")"
//         );
//     }

//     GenericFCN fcn(f);

//     M2::MnUserParameters upar;
//     for (std::size_t i = 0; i < x0.size(); ++i) {
//         upar.Add("x" + std::to_string(i), x0[i], safe_step(scales[i]));
//     }

//     M2::MnMigrad migrad(fcn, upar);
//     M2::FunctionMinimum min = migrad(max_fcn, tol_edm);

//     M2::MnHesse hesse;
//     hesse(fcn, min);

//     MinuitJointFit out;
//     out.fmin = min.Fval();
//     out.ok = min.IsValid();

//     const auto& st = min.UserState();
//     const std::size_t n = x0.size();

//     out.x_hat.resize(n);
//     out.x_err.assign(n, 0.0);

//     for (std::size_t i = 0; i < n; ++i) {
//         out.x_hat[i] = st.Value("x" + std::to_string(i));
//         out.x_err[i] = st.Error("x" + std::to_string(i));
//     }

//     out.cov = RealMatrix(n, n);
//     const auto& cov = st.Covariance();
//     for (std::size_t i = 0; i < n; ++i)
//         for (std::size_t j = 0; j < n; ++j)
//             out.cov.at(i, j) = cov(i, j);

//     return out;
// }


// class MinuitMLEstimatorLocal {
// public:
//     // using ModelFn = ProfiledLikelihood::ModelFn;

//     MinuitMLEstimatorLocal(LikelihoodContext ctx, ModelFn model, std::size_t max_fcn, double tol_edm)
//         : like_(std::move(ctx)), model_(std::move(model)), max_fcn_(max_fcn), tol_edm_(tol_edm) {}

//     FitResult fit_joint(const Vector& p0) const {
//         const std::size_t p_dim = p0.size();

//         // Vector eta0 = like_.nuisance_central_values;
//         Vector eta0;
//         for (auto elem : like_.nuis_defs) {
//             eta0.push_back(elem.value);
//         }
//         Vector eta_scales = like_.nuisance_dist->get_stds();

//         if (eta0.size() != eta_scales.size()) {
//             std::cerr << "[ERROR] eta0.size() != eta_scales.size(): "
//                       << eta0.size() << " vs " << eta_scales.size() << "\n";
//             throw std::runtime_error("Nuisance central values and nuisance stds dimensions mismatch");
//         }

//         std::vector<double> x0;
//         x0.reserve(p_dim + eta0.size());
//         x0.insert(x0.end(), p0.begin(), p0.end());
//         x0.insert(x0.end(), eta0.begin(), eta0.end());

//         Vector p_scales = p0;
//         for (auto& v : p_scales) v = std::fabs(v);

//         std::vector<double> scales;
//         scales.reserve(p_scales.size() + eta_scales.size());
//         scales.insert(scales.end(), p_scales.begin(), p_scales.end());
//         scales.insert(scales.end(), eta_scales.begin(), eta_scales.end());

//         auto f = [this, p_dim](const std::vector<double>& x) -> double {
//             Vector p(x.begin(), x.begin() + p_dim);
//             Vector eta(x.begin() + p_dim, x.end());
//             return nll(p, eta);
//         };

//         MinuitJointFit mj = minuit_migrad_hesse(f, x0, scales, max_fcn_, tol_edm_);

//         FitResult fr;
//         fr.ell_hat = mj.fmin;

//         if (!mj.ok) {
//             fr.p_hat = p0;
//             fr.eta_hat = eta0;
//             fr.p_hat_std = Vector(p_dim, 0.0);
//             fr.p_hat_correlations = RealMatrix(p_dim, p_dim);
//             return fr;
//         }

//         fr.p_hat.assign(mj.x_hat.begin(), mj.x_hat.begin() + p_dim);
//         fr.eta_hat.assign(mj.x_hat.begin() + p_dim, mj.x_hat.end());

//         RealMatrix cov_p(p_dim, p_dim);
//         for (std::size_t i = 0; i < p_dim; ++i)
//             for (std::size_t j = 0; j < p_dim; ++j)
//                 cov_p.at(i, j) = mj.cov.at(i, j);

//         fr.p_hat_std.assign(p_dim, 0.0);
//         for (std::size_t i = 0; i < p_dim; ++i)
//             fr.p_hat_std[i] = std::sqrt(std::max(0.0, cov_p.at(i, i)));

//         fr.p_hat_correlations = RealMatrix(p_dim, p_dim);
//         for (std::size_t i = 0; i < p_dim; ++i) {
//             for (std::size_t j = 0; j < p_dim; ++j) {
//                 double denom = fr.p_hat_std[i] * fr.p_hat_std[j];
//                 fr.p_hat_correlations.at(i, j) = (denom > 0.0) ? (cov_p.at(i, j) / denom) : 0.0;
//             }
//         }

//         return fr;
//     }

// private:
//     double nll(const Vector& p, const Vector& eta) const {
//         Vector pred = model_(p, eta);

//         Vector r(pred.size());
//         for (std::size_t i = 0; i < pred.size(); ++i)
//             r[i] = pred[i] - like_.exp_obs_values[i];

//         double ell_obs = like_.exp_obs_dist->logpdf(r);
//         double ell_eta = like_.nuisance_dist->logpdf(eta);

//         return -(ell_obs + ell_eta);
//     }
//     LikelihoodContext like_;
//     ModelFn model_;
//     std::size_t max_fcn_;
//     double tol_edm_;
// };


int main(int argc, char** argv) {

//     auto F_gauss_nll = [](const Vector& p) -> double {
//         double mu0 = 1e-12, mu1 = 1e3,  mu2 = 1e-4;
//         double s0  = 1e-13, s1  = 1e2,  s2  = 1e-6;

//         RealMatrix corr({
//             {1,    0.2, -0.5},
//             {0.2,  1,    0.7},
//             {-0.5, 0.7,  1}
//         });

//         RealMatrix z({
//             Vector{(p[0] - mu0) / s0},
//             Vector{(p[1] - mu1) / s1},
//             Vector{(p[2] - mu2) / s2}
//         });

//         return 0.5 * (2 * PI * corr.slogdet().logdet + (z.transpose() * corr.inv() * z).at(0, 0));
//     };

//     {
//         std::vector<double> x0     = {3e-12, 250.0, 0.0002};
//         std::vector<double> scales = {1e-13, 1e2,   1e-6};

//         auto f = [&](const std::vector<double>& x) -> double {
//             Vector p = {x[0], x[1], x[2]};
//             return F_gauss_nll(p);
//         };

//         auto r = minuit_migrad_hesse(f, x0, scales, /*max_fcn*/ 50000, /*tol_edm*/ 1e-10);

//         std::cout << "=== Minuit test (gauss corr) ===\n";
//         std::cout << "ok=" << r.ok << " fmin=" << std::setprecision(17) << r.fmin << "\n";
//         std::cout << "xhat = "; print_vec(r.x_hat);
//         std::cout << "xerr = "; print_vec(r.x_err);
//         std::cout << "\n";
//     }

//     HyperisoMaster hyp;
//     HyperisoConfig config_hyp;
//     config_hyp.model = Model::SM;
//     hyp.init("lha/si_input.flha", config_hyp);

//     std::shared_ptr<ObservableInterface> oint = std::make_shared<ObservableInterface>();
//     oint->add_observable(ObservableMapper::to_id(Observables::BR_BS_MUMU_UNTAG), QCDOrder::LO, true)
//         .add_observable(ObservableMapper::to_id(Observables::BR_BD_MUMU), QCDOrder::LO, true);

//     StatisticConfig config;
//     config.MC_draws     = 100;
//     config.MLE_max_iter = 10000;
//     config.MLE_tol      = 1e-6; 

//     std::vector<ParamId> p_specs = {
//         ParamId(ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::MATCHING), WCoefMapper::flha_full(WCoef::C10, QCDOrder::LO, ContributionType::SM))
//     };

//     std::shared_ptr<IStatParamOptimizerProxy> spop = std::make_shared<StatParamOptimizerProxy>();
//     auto model = std::make_shared<ObservableInterfaceProxy>(oint, spop);

//     std::shared_ptr<INuisancePathsProvider> npp = std::make_shared<DefaultNuisancePathsProvider>();

//     StatisticManager stat(
//         config,
//         model,
//         std::make_shared<StatCorrelationProxy>(),
//         std::make_shared<StatParameterProxy>(),
//         std::make_shared<StatParamSourcesProxy>(),
//         std::make_shared<StatDependencyPruner>(),
//         std::make_shared<NuisanceReader>(npp),
//         spop
//     );

//     LOG_INFO("YO1");
//     stat.update_cache(p_specs);

//     BlockProxy bp;
//     bp.log_all_blocks(ParameterType::WILSON);

//     auto start_u = std::chrono::steady_clock::now();
//     stat.compute_uncertainties();
//     auto stop_u  = std::chrono::steady_clock::now();
//     LOG_INFO("YO2");
//     auto us_u = std::chrono::duration_cast<std::chrono::microseconds>(stop_u - start_u).count();
//     std::cout << "Uncertainty estimation time : " << us_u << " µs\n";


//     stat.update_cache(p_specs);

//     auto p_specs_map = stat.get_p_specs(p_specs);
//     auto eta_specs_real = stat.get_all_obss_deps();
//     for (const auto& [pid, _] : p_specs_map) eta_specs_real.erase(pid);
//     auto exp_obs_map = stat.get_obs_exp();

//     auto unz_p   = unzip(p_specs_map);
//     auto unz_eta = unzip(eta_specs_real);
//     auto unz_obs = unzip(exp_obs_map);

//     std::vector<ParamId> p_ids = unz_p.ids;
//     std::vector<ParamId> eta_ids = unz_eta.ids;
//     std::vector<ExperimentObs> obs_ids = unz_obs.ids;

//     auto nuisance_dist = stat.build_nuisance_distribution();
//     auto exp_obs_dist  = stat.build_exp_data_distribution();

//     if (nuisance_dist->get_stds().size() != unz_eta.vals.size()) {
//         std::cerr << "[ERROR] nuisance_dist->get_stds().size()=" << nuisance_dist->get_stds().size()
//                   << " but eta central values size=" << unz_eta.vals.size() << "\n";
//         std::cerr << "=> appelle fill_cache() juste avant de construire ctx et évite de recalculer eta ailleurs.\n";
//         return 3;
//     }
//     if (exp_obs_dist->dim() != unz_obs.vals.size()) {
//         std::cerr << "[ERROR] exp_obs_dist->dim()=" << exp_obs_dist->dim()
//                   << " but exp obs values size=" << unz_obs.vals.size() << "\n";
//         return 4;
//     }

//     LikelihoodContext ctx;
//     ctx.nuisance_dist = std::move(nuisance_dist);
//     ctx.exp_obs_dist  = std::move(exp_obs_dist);
//     Vec _ = unz_eta.vals;
//     // ctx.nuisance_central_values = unz_eta.vals;
//     ctx.exp_obs_values          = unz_obs.vals;

//     auto model_fn = [model, obs_ids, p_ids, eta_ids](const Vec& p_vec, const Vec& eta_vec) -> Vec {
//         auto pred_map = model->predict_optimized(zip(p_ids, p_vec), zip(eta_ids, eta_vec));

//         Vec out;
//         out.reserve(obs_ids.size());

//         for (const auto& bid : obs_ids) {
//             const auto& vec = pred_map.at(bid.obs.s);

//             auto it = std::find_if(vec.begin(), vec.end(), [&](const ObservableValue& ov){
//                 auto bin = ov.bin.value_or(std::pair<double,double>{0.,0.});
//                 return bin == bid.obs.p;
//             });
//             if (it == vec.end()) throw std::runtime_error("Missing predicted observable/bin");
//             out.push_back(it->value);
//         }
//         return out;
//     };

//     auto start_m = std::chrono::steady_clock::now();
//     MinuitMLEstimatorLocal est(std::move(ctx), model_fn, config.MLE_max_iter, config.MLE_tol);
//     FitResult fr = est.fit_joint(unz_p.vals);
//     auto stop_m  = std::chrono::steady_clock::now();
//     auto us_m = std::chrono::duration_cast<std::chrono::microseconds>(stop_m - start_m).count();
//     std::cout << "MLE (Minuit) fitting time : " << us_m << " µs\n";

//     std::cout << "ell_hat = " << std::setprecision(17) << fr.ell_hat << "\n";

//     std::cout << "p_hat = ";
//     for (auto v : fr.p_hat) std::cout << v << " ";
//     std::cout << "\n";

//     std::cout << "p_hat_std = ";
//     for (auto v : fr.p_hat_std) std::cout << v << " ";
//     std::cout << "\n";

//     std::cout << "p_hat_correlations:\n" << fr.p_hat_correlations << "\n";

//     // std::cout << *StatParameterProxy(ParameterType::OBSERVABLE)
//     //                  .get_param("FOBS", LhaID("511_1_0_0_2_13_-13"))
//     //           << std::endl;

    return 0;
}