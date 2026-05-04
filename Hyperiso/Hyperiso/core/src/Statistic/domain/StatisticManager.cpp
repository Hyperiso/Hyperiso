// #include "StatisticManager.h"

// namespace {

// std::string param_name(const ParamId& pid) {
//     std::ostringstream oss;
//     oss << pid;
//     return oss.str();
// }

// fit_app::ParameterDefinition make_fit_param_def(const ParamId& pid, double value, double sigma_hint) {
//     fit_app::ParameterDefinition out;
//     out.name = param_name(pid);
//     out.value = value;
//     out.step_hint = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
//         ? sigma_hint
//         : std::max(1e-3, 0.01 * std::abs(value));

//     const std::string& nm = out.name;

//     if (nm.find("FCONST") != std::string::npos) {
//         out.limits = std::make_pair(0.05, 0.35);
//     }

//     return out;
// }

// fit_app::ParameterDefinition make_nuisance_param_def(const ParamId& pid, double value, double sigma_hint) {
//     fit_app::ParameterDefinition out;
//     out.name = param_name(pid);
//     out.value = value;

//     const double s = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
//         ? std::abs(sigma_hint)
//         : std::max(1e-3, 0.01 * std::abs(value));

//     out.step_hint = s;

//     const std::string& nm = out.name;

//     // TODO : Very ugly
//     if (nm.find("SMINPUTS:3") != std::string::npos) {
//         out.limits = std::make_pair(0.05, 0.30);
//     } else if (nm.find("MASS:") != std::string::npos ||
//                nm.find("FLIFE:") != std::string::npos ||
//                nm.find("FCONST:") != std::string::npos ||
//                nm.find("FMASS:") != std::string::npos ||
//                nm.find("SMINPUTS:5") != std::string::npos ||
//                nm.find("SMINPUTS:6") != std::string::npos) {
//         out.limits = std::make_pair(std::max(1e-12, value - 5.0 * s), value + 5.0 * s);
//     }

//     return out;
// }

// template <class PredMapT>
// std::vector<double> ordered_prediction_vector(
//     const std::vector<ExperimentObs>& obs_ids,
//     const PredMapT& pred_map)
// {
//     std::vector<double> out;
//     out.reserve(obs_ids.size());

//     for (const auto& bid : obs_ids) {
//         const auto& vec = pred_map.at(bid.obs.s);

//         auto it = std::find_if(vec.begin(), vec.end(), [&](const auto& ov) {
//             auto bin = ov.bin.value_or(std::pair<double, double>{0.0, 0.0});
//             return bin == bid.obs.p;
//         });

//         if (it == vec.end()) {
//             throw std::runtime_error("Missing predicted observable/bin.");
//         }

//         out.push_back(it->value);
//     }

//     return out;
// }

// }

// StatisticManager::StatisticManager(StatisticConfig config,
//                                    std::shared_ptr<IModel> obs_int,
//                                    std::shared_ptr<IStatCorrelationProxy> pscp,
//                                    std::shared_ptr<IStatParameterProxy> pspp,
//                                    std::shared_ptr<IStatSourcesProxy> sp,
//                                    std::shared_ptr<IStatDependencyPruner> dp,
//                                    std::shared_ptr<INuisanceReader> nuisance_reader)
//     : obs_int(std::move(obs_int)),
//       pscp(std::move(pscp)),
//       pspp(std::move(pspp)),
//       sp(std::move(sp)),
//       dp(std::move(dp)),
//       nuisance_reader_(std::move(nuisance_reader)),
//       config(std::move(config))
// {
//     if (!nuisance_reader_) {
//         throw std::invalid_argument("StatisticManager: nuisance_reader is null");
//     }

//     this->obs_int->compute_observables();
//     reload_nuisance_specs();
//     invalidate_fit_state();
// }

// // std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
// //     // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
// //     unsigned int seed = 123456u;
// //     std::map<ParamId, MarginalType> nuisance_marginals;

// //     for (auto& [pid, v] : cache.eta_specs_real) {
// //         nuisance_marginals.emplace(pid, MarginalType::GAUSSIAN);
// //     }

// //     for (auto& [pid, mt] : config.override_nuisance_marginals) {
// //         for (auto elem : nuisance_marginals) {
// //             std::cout << elem.first << std::endl;
// //             std::cout << ParameterTypeMapper::str(elem.first.type.value_or(ParameterType::WILSON)) << std::endl;
// //         }
// //         if (!nuisance_marginals.contains(pid)) {
// //             LOG_WARN("Parameter", pid, "is not a nuisance for the selected observables");
// //             continue;
// //         }
// //         nuisance_marginals.at(pid) = mt;
// //     }

// //     auto unzipped = unzip(nuisance_marginals);
// //     std::vector<ParamId> nuisance_ids = unzipped.ids;
// //     std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

// //     for (auto& [pid, mt] : nuisance_marginals) {
// //         MarginalConfig cfg = MarginalConfigFactory().create(pid, mt);
// //         auto m_ptr = DistributionFactory::create(mt, cfg, seed);
// //         marginals.emplace_back(std::move(m_ptr));
// //     }

// //     return marginals;
// // }

// std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
//     unsigned int seed = 123456u;
//     std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

//     marginals.reserve(cache.eta_specs_real.size());

//     for (const auto& [pid, _] : cache.eta_specs_real) {
//         const MarginalType mt = resolve_nuisance_marginal_type(pid);
//         MarginalConfig cfg = make_nuisance_marginal_config(pid, mt);
//         marginals.emplace_back(DistributionFactory::create(mt, cfg, seed));
//     }

//     return marginals;
// }

// std::unique_ptr<JointDistribution> StatisticManager::build_nuisance_distribution() {
//     // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
//     unsigned int seed = 123456u;

//     std::unique_ptr<ICopula> copula;
//     if (config.nuisance_copula_type == CopulaType::GAUSSIAN) {
//         GaussianCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
//         copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
//     } else if (config.nuisance_copula_type == CopulaType::STUDENT_T) {
//         StudentTCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
//         copula_cfg.nu = cache.eta_specs_real.size();
//         copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
//     }

//     return std::make_unique<JointDistribution>(
//         std::move(build_nuisance_marginal_distributions()), 
//         std::move(copula)
//     );
// }

// std::unique_ptr<JointDistribution> StatisticManager::build_exp_data_distribution() {
//     unsigned int seed = std::random_device{}();
//     std::map<ExperimentObs, MarginalType> exp_data_marginals;

//     for (auto& [oid, v] : cache.exp_obs) {
//         exp_data_marginals.emplace(oid, MarginalType::GAUSSIAN);
//     }

//     for (auto& [oid, mt] : config.override_exp_data_marginals) {
//         if (!exp_data_marginals.contains(oid)) {
//             LOG_WARN("Observable", ObservableMapper::str(oid.obs.s), "is not a taken into account in this fit.");
//             continue;
//         }
//         exp_data_marginals.at(oid) = mt;
//     }

//     auto unzipped = unzip(exp_data_marginals);
//     std::vector<ExperimentObs> obs_ids = unzipped.ids;
//     std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

//     for (auto& [oid, mt] : exp_data_marginals) {//TODO : checkkkkkk
//         MarginalConfig cfg = MarginalConfigFactory().create(oid, mt);
//         auto m_ptr = DistributionFactory::create(mt, cfg, seed);
//         marginals.emplace_back(std::move(m_ptr));
//     }

//     std::unique_ptr<ICopula> copula;
//     if (config.exp_data_copula_type == CopulaType::GAUSSIAN) {
//         GaussianCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
//         copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
//     } else if (config.exp_data_copula_type == CopulaType::STUDENT_T) {
//         StudentTCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
//         copula_cfg.nu = obs_int->n_observables() - 1;
//         copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
//     }

//     return std::make_unique<JointDistribution>(std::move(marginals), std::move(copula));
// }

// FitResultWithMaps StatisticManager::compute_MLE(const std::vector<ParamId>& p_specs) {
//     update_cache(p_specs);

//     if (cache.p_specs.empty()) {
//         throw std::invalid_argument("compute_MLE called with an empty fit parameter list.");
//     }

//     auto unzipped_fit_params = unzip(cache.p_specs);
//     auto unzipped_nuisances  = unzip(cache.eta_specs_real);
//     auto unzipped_exp_obs    = unzip(cache.exp_obs);

//     const std::vector<ParamId> p_ids = unzipped_fit_params.ids;
//     const std::vector<double> p0 = unzipped_fit_params.vals;

//     const std::vector<ParamId> eta_ids = unzipped_nuisances.ids;
//     const std::vector<double> eta0 = unzipped_nuisances.vals;

//     const std::vector<ExperimentObs> obs_ids = unzipped_exp_obs.ids;
//     const std::vector<double> exp_obs_vals = unzipped_exp_obs.vals;

//     auto ctx = std::make_shared<LikelihoodContext>();
//     ctx->nuisance_dist = build_nuisance_distribution();
//     ctx->exp_obs_dist = build_exp_data_distribution();
//     ctx->exp_obs_values = exp_obs_vals;

//     // ctx->fp_defs.reserve(p_ids.size());
//     // for (std::size_t i = 0; i < p_ids.size(); ++i) {
//     //     double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
//     //     ctx->fp_defs.emplace_back(make_param_def(p_ids[i], p0[i], sigma));
//     // }

//     // ctx->nuis_defs.reserve(eta_ids.size());
//     // for (std::size_t i = 0; i < eta_ids.size(); ++i) {
//     //     double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
//     //     ctx->nuis_defs.emplace_back(make_param_def(eta_ids[i], eta0[i], sigma));
//     // }

//     ctx->fp_defs.reserve(p_ids.size());
//     for (std::size_t i = 0; i < p_ids.size(); ++i) {
//         double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
//         ctx->fp_defs.emplace_back(make_fit_param_def(p_ids[i], p0[i], sigma));
//     }

//     // ctx->nuis_defs.reserve(eta_ids.size());
//     // for (std::size_t i = 0; i < eta_ids.size(); ++i) {
//     //     double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
//     //     ctx->nuis_defs.emplace_back(make_nuisance_param_def(eta_ids[i], eta0[i], sigma));
//     // }
//     ctx->nuis_defs.reserve(eta_ids.size());
//     for (std::size_t i = 0; i < eta_ids.size(); ++i) {
//         double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
//         ctx->nuis_defs.emplace_back(
//             make_nuisance_parameter_definition(eta_ids[i], eta0[i], sigma)
//         );
//     }

//     auto model_fn = [this, obs_ids, p_ids, eta_ids](
//         const std::vector<double>& p_vec,
//         const std::vector<double>& eta_vec) -> std::vector<double>
//     {
//         auto pred_map = this->obs_int->predict_optimized(
//             zip(p_ids, p_vec),
//             zip(eta_ids, eta_vec)
//         );

//         return ordered_prediction_vector(obs_ids, pred_map);
//     };

//     last_ctx_ = ctx;
//     last_like_ = std::make_shared<BaseLikelihood>(model_fn, ctx, p_ids.size());
//     last_fitter_ = std::make_shared<MLFitter>(ctx, model_fn);

//     last_fit_raw_ = last_fitter_->maximum_likelihood_fit(p0);

//     last_fit_param_ids_ = p_ids;
//     last_nuisance_ids_ = eta_ids;
//     last_fit_param_index_.clear();
//     for (std::size_t i = 0; i < p_ids.size(); ++i) {
//         last_fit_param_index_[p_ids[i]] = i;
//     }

//     FitResultWithMaps out;
//     out.fit_ok = !last_fit_raw_.p_hat.empty();
//     out.ell_hat = last_fit_raw_.ell_hat;
//     out.p_hat = zip(p_ids, last_fit_raw_.p_hat);
//     out.eta_hat = zip(eta_ids, last_fit_raw_.eta_hat);
//     out.p_hat_std = zip(p_ids, last_fit_raw_.p_hat_std);
//     out.p_correlations = zip(p_ids, last_fit_raw_.p_hat_correlations);

//     cache.mle_result = out;
//     return out;
// }

// Contour StatisticManager::confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, ContourOptions options) {
//     if (!cache.mle_result.fit_ok || !last_fitter_) {
//         throw std::runtime_error("Please run compute_MLE before requesting a confidence contour.");
//     }

//     if (!last_fit_param_index_.contains(p1) || !last_fit_param_index_.contains(p2)) {
//         throw std::invalid_argument("Contour requested for parameters that are not in the last fitted parameter set.");
//     }

//     if (p1 == p2) {
//         throw std::invalid_argument("Contour requires two distinct parameters.");
//     }

//     const std::size_t x_id = last_fit_param_index_.at(p1);
//     const std::size_t y_id = last_fit_param_index_.at(p2);

//     Contour cl = last_fitter_->contour(
//         x_id,
//         y_id,
//         z,
//         bounds,
//         options
//     );

//     return cl;
// }

// void StatisticManager::print_cache()
// {
//     for (auto elem : cache.eta_specs_real) {
//         std::cout << " eta_specs_real : " << elem.first << " = " << elem.second;
//     }
//     std::cout << std::endl;

//     for (auto elem : cache.SigmaEta) {
//         for (auto elem2 : elem.second) {
//             std::cout << " SigmaEta : " << elem.first << " | " << elem2.first << " = " << elem2.second;
//         }
//         std::cout << std::endl;
//     }

//     for (auto elem : cache.exp_obs) {
//         std::cout << " exp_obs : " << elem.first.str() << " = " << elem.second;
//     }
//     std::cout << std::endl;

//     for (auto elem : cache.SigmaObs) {
//         for (auto elem2 : elem.second) {
//             std::cout << " SigmaObs : " << elem.first.str() << " | " << elem2.first.str() << " = " << elem2.second;
//         }
//         std::cout << std::endl;
//     }

//     for (auto elem : cache.p_specs) {
//         std::cout << " p_specs : " << elem.first << " = " << elem.second;
//     }
//     std::cout << std::endl;

//     std::cout << "eta size : " << this->cache.eta_specs_real.size() << std::endl;
//     std::cout << "etasigma size : " << this->cache.eta_specs_real.size() << " | " << this->cache.SigmaEta.begin()->second.size() << std::endl;
//     std::cout << "p_specs size : " << this->cache.p_specs.size() << std::endl;
//     std::cout << "exp_obs : " << this->cache.exp_obs.size() << std::endl;
//     std::cout << "SigmaObs size : " << this->cache.SigmaObs.size() << " | " << this->cache.SigmaObs.begin()->second.size() << std::endl;
// }

// void StatisticManager::reload_nuisance_specs() {
//     default_nuisance_specs_ = nuisance_reader_->load_default();

//     if (current_user_nuisance_file_.has_value()) {
//         user_nuisance_specs_ = nuisance_reader_->load_user(*current_user_nuisance_file_);
//     } else {
//         user_nuisance_specs_ = nuisance_reader_->load_user();
//     }

//     rebuild_merged_nuisance_specs();
//     invalidate_fit_state();
// }

// void StatisticManager::set_nuisance_user_file(const fs::path& user_yaml_path) {
//     current_user_nuisance_file_ = user_yaml_path;
//     reload_nuisance_specs();
// }

// void StatisticManager::clear_nuisance_user_file() {
//     current_user_nuisance_file_.reset();
//     reload_nuisance_specs();
// }

// void StatisticManager::rebuild_merged_nuisance_specs() {
//     merged_nuisance_specs_ = default_nuisance_specs_;
//     for (const auto& [pid, spec] : user_nuisance_specs_) {
//         merged_nuisance_specs_[pid] = spec;
//     }
// }

// void StatisticManager::invalidate_fit_state() {
//     cache.mle_result = FitResultWithMaps{};

//     last_ctx_.reset();
//     last_like_.reset();
//     last_fitter_.reset();

//     last_fit_param_ids_.clear();
//     last_nuisance_ids_.clear();
//     last_fit_param_index_.clear();
// }

// const NuisanceSpec* StatisticManager::find_nuisance_spec(const ParamId& pid) const {
//     auto it = merged_nuisance_specs_.find(pid);
//     if (it == merged_nuisance_specs_.end()) {
//         return nullptr;
//     }
//     return &it->second;
// }

// MarginalType StatisticManager::resolve_nuisance_marginal_type(const ParamId& pid) const {
//     MarginalType mt = MarginalType::GAUSSIAN;

//     if (const auto* spec = find_nuisance_spec(pid)) {
//         mt = spec->marginal;
//     }

//     if (config.override_nuisance_marginals.contains(pid)) {
//         mt = config.override_nuisance_marginals.at(pid);
//     }

//     return mt;
// }

// fit_app::ParameterDefinition StatisticManager::make_nuisance_parameter_definition(const ParamId& pid,
//                                                                                   double value,
//                                                                                   double sigma_hint) const
// {
//     fit_app::ParameterDefinition out;
//     out.name = param_name(pid);
//     out.value = value;

//     const double s = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
//         ? std::abs(sigma_hint)
//         : std::max(1e-3, 0.01 * std::abs(value));

//     out.step_hint = s;

//     if (const auto* spec = find_nuisance_spec(pid)) {
//         out.limits = spec->bounds;
//         return out;
//     }

//     const std::string& nm = out.name;
//     if (nm.find("SMINPUTS:3") != std::string::npos) {
//         out.limits = std::make_pair(0.05, 0.30);
//     } else if (nm.find("MASS:") != std::string::npos ||
//                nm.find("FLIFE:") != std::string::npos ||
//                nm.find("FCONST:") != std::string::npos ||
//                nm.find("FMASS:") != std::string::npos ||
//                nm.find("SMINPUTS:5") != std::string::npos ||
//                nm.find("SMINPUTS:6") != std::string::npos) {
//         out.limits = std::make_pair(std::max(1e-12, value - 5.0 * s), value + 5.0 * s);
//     }

//     return out;
// }

// MarginalConfig StatisticManager::make_nuisance_marginal_config(const ParamId& pid,
//                                                                MarginalType mt) const
// {
//     // if (const auto* spec = find_nuisance_spec(pid)) {
//     //     return MarginalConfigFactory().create(pid, mt, *spec);
//     // }
//     return MarginalConfigFactory().create(pid, mt);
// }

// std::map<ParamId, double> StatisticManager::get_all_obss_deps() {
//     std::unordered_set<ParamId> eta_infos;

//     for (const auto& obsId : obs_int->get_obs_ids()) {
//         for (auto paramId : obs_int->get_obs_deps(obsId.s)) {
//             eta_infos.insert(paramId);
//         }
//     }

//     std::unordered_set<ParamId> eta_infos_leaf = this->sp->get_all_leaf_sources(eta_infos);
//     std::map<ParamId, double> eta_specs_real_leaf;
//     std::map<ParamId, double> delta_rel;

//     for (auto& paramId : eta_infos_leaf) {
//         double u = pspp->get_param(paramId)->get_combined_std().real();
//         double val = pspp->get_param(paramId)->get_val().real();
//         if (fpeq(val, 0.0) && !fpeq(u, 0.0)) {
//             eta_specs_real_leaf[paramId] = val;
//             continue;
//         }
            
//         if (!std::isfinite(u) || fpeq(std::abs(u), 0.0))
//             delta_rel[paramId] = 0;
//         else
//             delta_rel[paramId] = std::abs(u / val);
//     }
        
    
//     using T = std::pair<ParamId, double>;
//     double delta_rel_max = std::max_element(
//         delta_rel.begin(), delta_rel.end(), 
//         [] (const T& p, const T& q) { return p.second < q.second; }
//     )->second;

//     for (auto& paramId : eta_infos_leaf) {
//         LOG_INFO("Compared to max relative uncertainty", paramId, delta_rel[paramId] / delta_rel_max);
//     }

//     for (auto& paramId : eta_infos_leaf) {
//         if (delta_rel[paramId] / delta_rel_max > config.nuisance_relevance_cutoff)
//             eta_specs_real_leaf[paramId] = pspp->get_param(paramId)->get_val();
//     }

//     LOG_INFO("Significant nuisances");
//     for (const auto& [pid, val] : eta_specs_real_leaf) {
//         LOG_INFO(pid, val);
//     }

//     return eta_specs_real_leaf;
// }

// std::map<ParamId, double> StatisticManager::get_p_specs(const std::vector<ParamId>& p_specs) {
//     std::map<ParamId, double> out;
//     for (auto elem : p_specs) {            
//         out[elem] = pspp->get_param(elem)->get_val();
//     }
//     return out;
// }
// std::map<ParamId, std::map<ParamId, double>> StatisticManager::get_all_correlations() {
//     std::map<ParamId, std::map<ParamId, double>> res;
//     CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);
//     res = ct.transform(this->cache.eta_specs_real);
//     return res;
// }

// std::map<ExperimentObs, std::map<ExperimentObs, double>> StatisticManager::get_all_obs_correlations() {
//     std::map<ExperimentObs, std::map<ExperimentObs, double>> res;
//     CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);

//     res = ct.transform(this->cache.exp_obs);
//     return res;
// }

// std::map<ExperimentObs, double> StatisticManager::get_obs_exp() {
//     std::map<ExperimentObs, double> out;

//     for (const auto& obsId : obs_int->get_obs_ids()) {
//         auto _ = pspp->get_obs_param(obsId);
//         for (auto _2 : _) {
//             out[_2.first] = _2.second->get_val(); 
//         }
//         // out[obsId] = pspp->get_obs_param(obsId)->get_val();
//     }
//     return out;
// }

// #include "StatisticManager.h"

// #include <limits>
// #include <numeric>
// #include <set>


// bool starts_with(const std::string& s, const std::string& prefix) {
//     return s.size() >= prefix.size() &&
//            std::equal(prefix.begin(), prefix.end(), s.begin());
// }

// std::string param_name(const ParamId& pid) {
//     std::ostringstream oss;
//     oss << pid;
//     return oss.str();
// }
// std::vector<std::size_t> select_indices_by_prefix(const std::vector<ParamId>& ids,
//                                                   const std::string& prefix)
// {
//     std::vector<std::size_t> out;
//     for (std::size_t i = 0; i < ids.size(); ++i) {
//         if (starts_with(param_name(ids[i]), prefix)) {
//             out.push_back(i);
//         }
//     }
//     return out;
// }

// std::vector<std::size_t> select_indices_by_names(const std::vector<ParamId>& ids,
//                                                  const std::set<std::string>& names)
// {
//     std::vector<std::size_t> out;
//     for (std::size_t i = 0; i < ids.size(); ++i) {
//         if (names.contains(param_name(ids[i]))) {
//             out.push_back(i);
//         }
//     }
//     return out;
// }

// RealMatrix principal_submatrix(const RealMatrix& M,
//                                const std::vector<std::size_t>& idx)
// {
//     RealMatrix out(idx.size(), idx.size());
//     for (std::size_t i = 0; i < idx.size(); ++i) {
//         for (std::size_t j = 0; j < idx.size(); ++j) {
//             out.at(i, j) = M.at(idx[i], idx[j]);
//         }
//     }
//     return out;
// }

// void log_smallest_mode(const RealMatrix& M,
//                        const std::vector<ParamId>& ids,
//                        const std::vector<std::size_t>& idx,
//                        const std::string& tag)
// {
//     std::cout << "[FIT] ---- " << tag << " ----\n";
//     std::cout << "[FIT] size=" << idx.size() << "\n";
//     if (idx.empty()) {
//         std::cout << "[FIT] empty block\n";
//         return;
//     }

//     RealMatrix S = principal_submatrix(M, idx);

//     double diag_min = std::numeric_limits<double>::infinity();
//     double diag_max = -std::numeric_limits<double>::infinity();
//     double max_abs_offdiag = 0.0;

//     for (std::size_t i = 0; i < S.rows(); ++i) {
//         diag_min = std::min(diag_min, S.at(i, i));
//         diag_max = std::max(diag_max, S.at(i, i));
//         for (std::size_t j = i + 1; j < S.cols(); ++j) {
//             max_abs_offdiag = std::max(max_abs_offdiag, std::abs(S.at(i, j)));
//         }
//     }

//     std::cout << "[FIT] diag_min=" << diag_min
//               << " diag_max=" << diag_max
//               << " max_abs_offdiag=" << max_abs_offdiag << "\n";

//     EigenSystem e = S.eig();

//     std::size_t imin = 0;
//     double evmin = e.D.at(0, 0);
//     for (std::size_t i = 1; i < S.rows(); ++i) {
//         if (e.D.at(i, i) < evmin) {
//             evmin = e.D.at(i, i);
//             imin = i;
//         }
//     }

//     std::cout << "[FIT] smallest_eigenvalue=" << evmin << "\n";

//     std::vector<std::pair<double, std::size_t>> comps;
//     for (std::size_t i = 0; i < idx.size(); ++i) {
//         comps.push_back({std::abs(e.P.at(i, imin)), i});
//     }

//     std::sort(comps.begin(), comps.end(),
//               [](const auto& a, const auto& b) { return a.first > b.first; });

//     std::cout << "[FIT] dominant entries of smallest mode:\n";
//     for (std::size_t k = 0; k < std::min<std::size_t>(10, comps.size()); ++k) {
//         const std::size_t local_i = comps[k].second;
//         const std::size_t global_i = idx[local_i];
//         std::cout << "[FIT]   " << ids[global_i]
//                   << " coeff=" << e.P.at(local_i, imin) << "\n";
//     }

//     std::cout << "[FIT] matrix entries:\n";
//     for (std::size_t i = 0; i < idx.size(); ++i) {
//         for (std::size_t j = 0; j < idx.size(); ++j) {
//             std::cout << "[FIT]   "
//                       << ids[idx[i]] << " | " << ids[idx[j]]
//                       << " = " << S.at(i, j) << "\n";
//         }
//     }
// }


// static void dump_matrix_sanity(const RealMatrix& M,
//                                const std::vector<ParamId>& ids,
//                                const std::string& name)
// {
//     std::cout << "[FIT] Matrix " << name
//               << " shape=" << M.rows() << "x" << M.cols()
//               << " symmetric=" << M.is_symmetric() << "\n";

//     double diag_min = std::numeric_limits<double>::infinity();
//     double diag_max = -std::numeric_limits<double>::infinity();
//     std::size_t nonpos_diag = 0;

//     double max_abs_offdiag = 0.0;
//     std::size_t n_abs_gt_1 = 0;
//     std::pair<std::size_t,std::size_t> argmax_offdiag{0,0};

//     for (std::size_t i = 0; i < M.rows(); ++i) {
//         const double d = M.at(i,i);
//         diag_min = std::min(diag_min, d);
//         diag_max = std::max(diag_max, d);
//         if (!(d > 0.0)) ++nonpos_diag;

//         for (std::size_t j = i + 1; j < M.cols(); ++j) {
//             const double a = std::abs(M.at(i,j));
//             if (a > max_abs_offdiag) {
//                 max_abs_offdiag = a;
//                 argmax_offdiag = {i,j};
//             }
//             if (a > 1.0 + 1e-10) ++n_abs_gt_1;
//         }
//     }

//     std::cout << "[FIT] " << name
//               << " diag_min=" << diag_min
//               << " diag_max=" << diag_max
//               << " nonpos_diag=" << nonpos_diag
//               << " max_abs_offdiag=" << max_abs_offdiag
//               << " n_abs_offdiag_gt_1=" << n_abs_gt_1
//               << "\n";

//     if (M.rows() == M.cols() && M.is_symmetric()) {
//         EigenSystem e = M.eig();

//         std::size_t imin = 0;
//         double evmin = e.D.at(0,0);
//         for (std::size_t i = 1; i < M.rows(); ++i) {
//             if (e.D.at(i,i) < evmin) {
//                 evmin = e.D.at(i,i);
//                 imin = i;
//             }
//         }

//         std::vector<std::pair<double,std::size_t>> comps;
//         for (std::size_t i = 0; i < M.rows(); ++i) {
//             comps.push_back({std::abs(e.P.at(i, imin)), i});
//         }
//         std::sort(comps.begin(), comps.end(),
//                   [](const auto& a, const auto& b) { return a.first > b.first; });

//         std::cout << "[FIT] " << name << " smallest_eigenvalue=" << evmin << "\n";
//         std::cout << "[FIT] " << name << " dominant entries of smallest mode:\n";
//         for (std::size_t k = 0; k < std::min<std::size_t>(10, comps.size()); ++k) {
//             std::size_t i = comps[k].second;
//             std::cout << "[FIT]   " << ids[i]
//                       << " coeff=" << e.P.at(i, imin) << "\n";
//         }
//     }
// }

// namespace {

// // std::string param_name(const ParamId& pid) {
// //     std::ostringstream oss;
// //     oss << pid;
// //     return oss.str();
// // }

// fit_app::ParameterDefinition make_fit_param_def(const ParamId& pid, double value, double sigma_hint) {
//     fit_app::ParameterDefinition out;
//     out.name = param_name(pid);
//     out.value = value;
//     out.step_hint = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
//         ? sigma_hint
//         : std::max(1e-3, 0.01 * std::abs(value));

//     const std::string& nm = out.name;

//     if (nm.find("FCONST") != std::string::npos) {
//         out.limits = std::make_pair(0.05, 0.35);
//     }

//     return out;
// }

// fit_app::ParameterDefinition make_nuisance_param_def(const ParamId& pid, double value, double sigma_hint) {
//     fit_app::ParameterDefinition out;
//     out.name = param_name(pid);
//     out.value = value;

//     const double s = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
//         ? std::abs(sigma_hint)
//         : std::max(1e-3, 0.01 * std::abs(value));

//     out.step_hint = s;

//     const std::string& nm = out.name;

//     // TODO : Very ugly
//     if (nm.find("SMINPUTS:3") != std::string::npos) {
//         out.limits = std::make_pair(0.05, 0.30);
//     } else if (nm.find("MASS:") != std::string::npos ||
//                nm.find("FLIFE:") != std::string::npos ||
//                nm.find("FCONST:") != std::string::npos ||
//                nm.find("FMASS:") != std::string::npos ||
//                nm.find("SMINPUTS:5") != std::string::npos ||
//                nm.find("SMINPUTS:6") != std::string::npos) {
//         out.limits = std::make_pair(std::max(1e-12, value - 5.0 * s), value + 5.0 * s);
//     }

//     return out;
// }

// template <class PredMapT>
// std::vector<double> ordered_prediction_vector(
//     const std::vector<ExperimentObs>& obs_ids,
//     const PredMapT& pred_map)
// {
//     std::vector<double> out;
//     out.reserve(obs_ids.size());

//     for (const auto& bid : obs_ids) {
//         const auto& vec = pred_map.at(bid.obs.s);

//         auto it = std::find_if(vec.begin(), vec.end(), [&](const auto& ov) {
//             auto bin = ov.bin.value_or(std::pair<double, double>{0.0, 0.0});
//             return bin == bid.obs.p;
//         });

//         if (it == vec.end()) {
//             throw std::runtime_error("Missing predicted observable/bin.");
//         }

//         out.push_back(it->value);
//     }

//     return out;
// }


// void log_selected_nuisance_diagnostics(const std::vector<ParamId>& eta_ids,
//                                       const std::map<ParamId, std::map<ParamId, double>>& sigma_eta) {
//     std::cout << "[FIT] Selected nuisance parameters: " << eta_ids.size() << std::endl;
//     for (const auto& pid : eta_ids) {
//         std::cout << "[FIT]   nuisance " << pid << std::endl;
//     }

//     if (sigma_eta.empty()) {
//         std::cout << "[FIT] SigmaEta is empty" << std::endl;
//         return;
//     }

//     try {
//         RealMatrix sigma = RealMatrix(unzip(sigma_eta).vals);
//         std::cout << "[FIT] SigmaEta shape=" << sigma.rows() << "x" << sigma.cols()
//                   << ", symmetric=" << sigma.is_symmetric() << std::endl;

//         RealMatrix sym = 0.5 * (sigma + sigma.transpose());
//         EigenSystem eig = sym.eig();

//         double min_eig = std::numeric_limits<double>::infinity();
//         double max_eig = -std::numeric_limits<double>::infinity();
//         std::size_t non_pos = 0;
//         std::size_t tiny = 0;

//         for (std::size_t i = 0; i < eig.D.rows(); ++i) {
//             const double lambda = eig.D.at(i, i);
//             min_eig = std::min(min_eig, lambda);
//             max_eig = std::max(max_eig, lambda);
//             if (!(lambda > 0.0)) ++non_pos;
//             if (std::abs(lambda) < 1e-12) ++tiny;
//         }

//         std::cout << "[FIT] SigmaEta eig_min=" << min_eig
//                   << ", eig_max=" << max_eig
//                   << ", non_positive=" << non_pos
//                   << ", tiny(|eig|<1e-12)=" << tiny
//                   << std::endl;
//     } catch (const std::exception& e) {
//         std::cout << "[FIT] SigmaEta diagnostics failed: " << e.what() << std::endl;
//     }
// }

// }

// StatisticManager::StatisticManager(StatisticConfig config,
//                                    std::shared_ptr<IModel> obs_int,
//                                    std::shared_ptr<IStatCorrelationProxy> pscp,
//                                    std::shared_ptr<IStatParameterProxy> pspp,
//                                    std::shared_ptr<IStatSourcesProxy> sp,
//                                    std::shared_ptr<IStatDependencyPruner> dp,
//                                    std::shared_ptr<INuisanceReader> nuisance_reader)
//     : obs_int(std::move(obs_int)),
//       pscp(std::move(pscp)),
//       pspp(std::move(pspp)),
//       sp(std::move(sp)),
//       dp(std::move(dp)),
//       nuisance_reader_(std::move(nuisance_reader)),
//       config(std::move(config))
// {
//     if (!nuisance_reader_) {
//         throw std::invalid_argument("StatisticManager: nuisance_reader is null");
//     }

//     this->obs_int->compute_observables();
//     reload_nuisance_specs();
//     invalidate_fit_state();
// }

// // std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
// //     // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
// //     unsigned int seed = 123456u;
// //     std::map<ParamId, MarginalType> nuisance_marginals;

// //     for (auto& [pid, v] : cache.eta_specs_real) {
// //         nuisance_marginals.emplace(pid, MarginalType::GAUSSIAN);
// //     }

// //     for (auto& [pid, mt] : config.override_nuisance_marginals) {
// //         for (auto elem : nuisance_marginals) {
// //             std::cout << elem.first << std::endl;
// //             std::cout << ParameterTypeMapper::str(elem.first.type.value_or(ParameterType::WILSON)) << std::endl;
// //         }
// //         if (!nuisance_marginals.contains(pid)) {
// //             LOG_WARN("Parameter", pid, "is not a nuisance for the selected observables");
// //             continue;
// //         }
// //         nuisance_marginals.at(pid) = mt;
// //     }

// //     auto unzipped = unzip(nuisance_marginals);
// //     std::vector<ParamId> nuisance_ids = unzipped.ids;
// //     std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

// //     for (auto& [pid, mt] : nuisance_marginals) {
// //         MarginalConfig cfg = MarginalConfigFactory().create(pid, mt);
// //         auto m_ptr = DistributionFactory::create(mt, cfg, seed);
// //         marginals.emplace_back(std::move(m_ptr));
// //     }

// //     return marginals;
// // }

// std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
//     unsigned int seed = 123456u;
//     std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

//     marginals.reserve(cache.eta_specs_real.size());

//     for (const auto& [pid, _] : cache.eta_specs_real) {
//         const MarginalType mt = resolve_nuisance_marginal_type(pid);
//         MarginalConfig cfg = make_nuisance_marginal_config(pid, mt);
//         marginals.emplace_back(DistributionFactory::create(mt, cfg, seed));
//     }

//     return marginals;
// }

// std::unique_ptr<JointDistribution> StatisticManager::build_nuisance_distribution() {
//     // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
//     unsigned int seed = 123456u;

//     std::unique_ptr<ICopula> copula;
//     if (config.nuisance_copula_type == CopulaType::GAUSSIAN) {
//         GaussianCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
//         copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
//     } else if (config.nuisance_copula_type == CopulaType::STUDENT_T) {
//         StudentTCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
//         copula_cfg.nu = cache.eta_specs_real.size();
//         copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
//     }

//     return std::make_unique<JointDistribution>(
//         std::move(build_nuisance_marginal_distributions()), 
//         std::move(copula)
//     );
// }

// std::unique_ptr<JointDistribution> StatisticManager::build_exp_data_distribution() {
//     unsigned int seed = std::random_device{}();
//     std::map<ExperimentObs, MarginalType> exp_data_marginals;

//     for (auto& [oid, v] : cache.exp_obs) {
//         exp_data_marginals.emplace(oid, MarginalType::GAUSSIAN);
//     }

//     for (auto& [oid, mt] : config.override_exp_data_marginals) {
//         if (!exp_data_marginals.contains(oid)) {
//             LOG_WARN("Observable", ObservableMapper::str(oid.obs.s), "is not a taken into account in this fit.");
//             continue;
//         }
//         exp_data_marginals.at(oid) = mt;
//     }

//     auto unzipped = unzip(exp_data_marginals);
//     std::vector<ExperimentObs> obs_ids = unzipped.ids;
//     std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

//     for (auto& [oid, mt] : exp_data_marginals) {//TODO : checkkkkkk
//         MarginalConfig cfg = MarginalConfigFactory().create(oid, mt);
//         auto m_ptr = DistributionFactory::create(mt, cfg, seed);
//         marginals.emplace_back(std::move(m_ptr));
//     }

//     std::unique_ptr<ICopula> copula;
//     if (config.exp_data_copula_type == CopulaType::GAUSSIAN) {
//         GaussianCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
//         copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
//     } else if (config.exp_data_copula_type == CopulaType::STUDENT_T) {
//         StudentTCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
//         copula_cfg.nu = obs_int->n_observables() - 1;
//         copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
//     }

//     return std::make_unique<JointDistribution>(std::move(marginals), std::move(copula));
// }

// static void scan_leave_one_out(const RealMatrix& M,
//                                const std::vector<ParamId>& ids,
//                                const std::vector<std::size_t>& idx,
//                                const std::string& tag)
// {
//     std::cout << "[FIT] ---- leave-one-out " << tag << " ----\n";
//     if (idx.size() <= 1) return;

//     for (std::size_t k = 0; k < idx.size(); ++k) {
//         std::vector<std::size_t> sub;
//         sub.reserve(idx.size() - 1);
//         for (std::size_t i = 0; i < idx.size(); ++i) {
//             if (i != k) sub.push_back(idx[i]);
//         }

//         RealMatrix S = principal_submatrix(M, sub);
//         EigenSystem e = S.eig();

//         double evmin = e.D.at(0,0);
//         for (std::size_t i = 1; i < S.rows(); ++i) {
//             evmin = std::min(evmin, e.D.at(i,i));
//         }

//         std::cout << "[FIT] remove " << ids[idx[k]]
//                   << " => smallest_eigenvalue = " << evmin << "\n";
//     }
// }

// FitResultWithMaps StatisticManager::compute_MLE(const std::vector<ParamId>& p_specs) {
//     update_cache(p_specs);

//     if (cache.p_specs.empty()) {
//         throw std::invalid_argument("compute_MLE called with an empty fit parameter list.");
//     }

//     auto unzipped_fit_params = unzip(cache.p_specs);
//     auto unzipped_nuisances  = unzip(cache.eta_specs_real);
//     auto unzipped_exp_obs    = unzip(cache.exp_obs);

//     {auto eta_ids = unzip(cache.eta_specs_real).ids;
//     RealMatrix Reta(unzip(cache.SigmaEta).vals);
//     dump_matrix_sanity(Reta, eta_ids, "SigmaEta");

// }

//     const std::vector<ParamId> p_ids = unzipped_fit_params.ids;
//     const std::vector<double> p0 = unzipped_fit_params.vals;

//     const std::vector<ParamId> eta_ids = unzipped_nuisances.ids;
//     const std::vector<double> eta0 = unzipped_nuisances.vals;
    
//         RealMatrix SigmaEta(unzip(cache.SigmaEta).vals);

//     // bloc complet
//     {
//         std::vector<std::size_t> all_idx(eta_ids.size());
//         std::iota(all_idx.begin(), all_idx.end(), 0);
//         log_smallest_mode(SigmaEta, eta_ids, all_idx, "SigmaEta full");
//     }

//     // bloc B_Ks:1_*
//     {
//         auto idx = select_indices_by_prefix(eta_ids, "B_Ks:1_");
//         log_smallest_mode(SigmaEta, eta_ids, idx, "SigmaEta B_Ks:1_*");
//     }

//     // bloc B_Ks:18_*
//     {
//         auto idx = select_indices_by_prefix(eta_ids, "B_Ks:18_");
//         log_smallest_mode(SigmaEta, eta_ids, idx, "SigmaEta B_Ks:18_*");
//     }

//     // bloc mixte 1_* + {7_*, 8_*, 14}
//     {
//         std::set<std::string> names = {
//             "B_Ks:7_1", "B_Ks:7_2",
//             "B_Ks:8_1", "B_Ks:8_2",
//             "B_Ks:14"
//         };

//         auto idx1 = select_indices_by_prefix(eta_ids, "B_Ks:1_");
//         auto idx2 = select_indices_by_names(eta_ids, names);

//         std::vector<std::size_t> idx = idx1;
//         idx.insert(idx.end(), idx2.begin(), idx2.end());
//         std::sort(idx.begin(), idx.end());
//         idx.erase(std::unique(idx.begin(), idx.end()), idx.end());

//         log_smallest_mode(SigmaEta, eta_ids, idx, "SigmaEta B_Ks:1_* + {7,8,14}");
//     }

//     {
//         std::set<std::string> names = {
//             "B_Ks:1_2_1",
//             "B_Ks:1_5_1",
//             "B_Ks:1_5_2",
//             "B_Ks:1_6_2",
//             "B_Ks:1_2_2",
//             "B_Ks:1_4_1",
//             "B_Ks:1_3_1",
//             "B_Ks:1_1_1",
//             "B_Ks:1_4_2",
//             "B_Ks:1_3_2"
//         };
//         auto idx = select_indices_by_names(eta_ids, names);
//         log_smallest_mode(SigmaEta, eta_ids, idx, "SigmaEta suspects");
//         scan_leave_one_out(SigmaEta, eta_ids, idx, "SigmaEta suspects");
//     }

//     log_selected_nuisance_diagnostics(eta_ids, cache.SigmaEta);

//     const std::vector<ExperimentObs> obs_ids = unzipped_exp_obs.ids;
//     const std::vector<double> exp_obs_vals = unzipped_exp_obs.vals;

//     auto ctx = std::make_shared<LikelihoodContext>();
//     ctx->nuisance_dist = build_nuisance_distribution();
//     ctx->exp_obs_dist = build_exp_data_distribution();
//     ctx->exp_obs_values = exp_obs_vals;

//     // ctx->fp_defs.reserve(p_ids.size());
//     // for (std::size_t i = 0; i < p_ids.size(); ++i) {
//     //     double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
//     //     ctx->fp_defs.emplace_back(make_param_def(p_ids[i], p0[i], sigma));
//     // }

//     // ctx->nuis_defs.reserve(eta_ids.size());
//     // for (std::size_t i = 0; i < eta_ids.size(); ++i) {
//     //     double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
//     //     ctx->nuis_defs.emplace_back(make_param_def(eta_ids[i], eta0[i], sigma));
//     // }

//     ctx->fp_defs.reserve(p_ids.size());
//     for (std::size_t i = 0; i < p_ids.size(); ++i) {
//         double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
//         ctx->fp_defs.emplace_back(make_fit_param_def(p_ids[i], p0[i], sigma));
//     }

//     // ctx->nuis_defs.reserve(eta_ids.size());
//     // for (std::size_t i = 0; i < eta_ids.size(); ++i) {
//     //     double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
//     //     ctx->nuis_defs.emplace_back(make_nuisance_param_def(eta_ids[i], eta0[i], sigma));
//     // }
//     ctx->nuis_defs.reserve(eta_ids.size());
//     for (std::size_t i = 0; i < eta_ids.size(); ++i) {
//         double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
//         ctx->nuis_defs.emplace_back(
//             make_nuisance_parameter_definition(eta_ids[i], eta0[i], sigma)
//         );
//     }

//     auto model_fn = [this, obs_ids, p_ids, eta_ids](
//         const std::vector<double>& p_vec,
//         const std::vector<double>& eta_vec) -> std::vector<double>
//     {
//         auto pred_map = this->obs_int->predict_optimized(
//             zip(p_ids, p_vec),
//             zip(eta_ids, eta_vec)
//         );

//         return ordered_prediction_vector(obs_ids, pred_map);
//     };

//     last_ctx_ = ctx;
//     last_like_ = std::make_shared<BaseLikelihood>(model_fn, ctx, p_ids.size());
//     last_fitter_ = std::make_shared<MLFitter>(ctx, model_fn);

//     last_fit_raw_ = last_fitter_->maximum_likelihood_fit(p0);

//     last_fit_param_ids_ = p_ids;
//     last_nuisance_ids_ = eta_ids;
//     last_fit_param_index_.clear();
//     for (std::size_t i = 0; i < p_ids.size(); ++i) {
//         last_fit_param_index_[p_ids[i]] = i;
//     }

//     FitResultWithMaps out;
//     out.fit_ok = !last_fit_raw_.p_hat.empty();
//     out.ell_hat = last_fit_raw_.ell_hat;
//     out.p_hat = zip(p_ids, last_fit_raw_.p_hat);
//     out.eta_hat = zip(eta_ids, last_fit_raw_.eta_hat);
//     out.p_hat_std = zip(p_ids, last_fit_raw_.p_hat_std);
//     out.p_correlations = zip(p_ids, last_fit_raw_.p_hat_correlations);

//     cache.mle_result = out;
//     return out;
// }

// Contour StatisticManager::confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, ContourOptions options) {
//     if (!cache.mle_result.fit_ok || !last_fitter_) {
//         throw std::runtime_error("Please run compute_MLE before requesting a confidence contour.");
//     }

//     if (!last_fit_param_index_.contains(p1) || !last_fit_param_index_.contains(p2)) {
//         throw std::invalid_argument("Contour requested for parameters that are not in the last fitted parameter set.");
//     }

//     if (p1 == p2) {
//         throw std::invalid_argument("Contour requires two distinct parameters.");
//     }

//     const std::size_t x_id = last_fit_param_index_.at(p1);
//     const std::size_t y_id = last_fit_param_index_.at(p2);

//     Contour cl = last_fitter_->contour(
//         x_id,
//         y_id,
//         z,
//         bounds,
//         options
//     );

//     return cl;
// }

// void StatisticManager::print_cache()
// {
//     for (auto elem : cache.eta_specs_real) {
//         std::cout << " eta_specs_real : " << elem.first << " = " << elem.second;
//     }
//     std::cout << std::endl;

//     for (auto elem : cache.SigmaEta) {
//         for (auto elem2 : elem.second) {
//             std::cout << " SigmaEta : " << elem.first << " | " << elem2.first << " = " << elem2.second;
//         }
//         std::cout << std::endl;
//     }

//     for (auto elem : cache.exp_obs) {
//         std::cout << " exp_obs : " << elem.first.str() << " = " << elem.second;
//     }
//     std::cout << std::endl;

//     for (auto elem : cache.SigmaObs) {
//         for (auto elem2 : elem.second) {
//             std::cout << " SigmaObs : " << elem.first.str() << " | " << elem2.first.str() << " = " << elem2.second;
//         }
//         std::cout << std::endl;
//     }

//     for (auto elem : cache.p_specs) {
//         std::cout << " p_specs : " << elem.first << " = " << elem.second;
//     }
//     std::cout << std::endl;

//     std::cout << "eta size : " << this->cache.eta_specs_real.size() << std::endl;
//     std::cout << "etasigma size : " << this->cache.eta_specs_real.size() << " | " << this->cache.SigmaEta.begin()->second.size() << std::endl;
//     std::cout << "p_specs size : " << this->cache.p_specs.size() << std::endl;
//     std::cout << "exp_obs : " << this->cache.exp_obs.size() << std::endl;
//     std::cout << "SigmaObs size : " << this->cache.SigmaObs.size() << " | " << this->cache.SigmaObs.begin()->second.size() << std::endl;
// }

// void StatisticManager::reload_nuisance_specs() {
//     default_nuisance_specs_ = nuisance_reader_->load_default();

//     if (current_user_nuisance_file_.has_value()) {
//         user_nuisance_specs_ = nuisance_reader_->load_user(*current_user_nuisance_file_);
//     } else {
//         user_nuisance_specs_ = nuisance_reader_->load_user();
//     }

//     rebuild_merged_nuisance_specs();
//     invalidate_fit_state();
// }

// void StatisticManager::set_nuisance_user_file(const fs::path& user_yaml_path) {
//     current_user_nuisance_file_ = user_yaml_path;
//     reload_nuisance_specs();
// }

// void StatisticManager::clear_nuisance_user_file() {
//     current_user_nuisance_file_.reset();
//     reload_nuisance_specs();
// }

// void StatisticManager::rebuild_merged_nuisance_specs() {
//     merged_nuisance_specs_ = default_nuisance_specs_;
//     for (const auto& [pid, spec] : user_nuisance_specs_) {
//         merged_nuisance_specs_[pid] = spec;
//     }
// }

// void StatisticManager::invalidate_fit_state() {
//     cache.mle_result = FitResultWithMaps{};

//     last_ctx_.reset();
//     last_like_.reset();
//     last_fitter_.reset();

//     last_fit_param_ids_.clear();
//     last_nuisance_ids_.clear();
//     last_fit_param_index_.clear();
// }

// const NuisanceSpec* StatisticManager::find_nuisance_spec(const ParamId& pid) const {
//     auto it = merged_nuisance_specs_.find(pid);
//     if (it == merged_nuisance_specs_.end()) {
//         return nullptr;
//     }
//     return &it->second;
// }

// MarginalType StatisticManager::resolve_nuisance_marginal_type(const ParamId& pid) const {
//     MarginalType mt = MarginalType::GAUSSIAN;

//     if (const auto* spec = find_nuisance_spec(pid)) {
//         mt = spec->marginal;
//     }

//     if (config.override_nuisance_marginals.contains(pid)) {
//         mt = config.override_nuisance_marginals.at(pid);
//     }

//     return mt;
// }

// fit_app::ParameterDefinition StatisticManager::make_nuisance_parameter_definition(const ParamId& pid,
//                                                                                   double value,
//                                                                                   double sigma_hint) const
// {
//     fit_app::ParameterDefinition out;
//     out.name = param_name(pid);
//     out.value = value;

//     const double s = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
//         ? std::abs(sigma_hint)
//         : std::max(1e-3, 0.01 * std::abs(value));

//     out.step_hint = s;

//     if (const auto* spec = find_nuisance_spec(pid)) {
//         out.limits = spec->bounds;
//         return out;
//     }

//     const std::string& nm = out.name;
//     if (nm.find("SMINPUTS:3") != std::string::npos) {
//         out.limits = std::make_pair(0.05, 0.30);
//     } else if (nm.find("MASS:") != std::string::npos ||
//                nm.find("FLIFE:") != std::string::npos ||
//                nm.find("FCONST:") != std::string::npos ||
//                nm.find("FMASS:") != std::string::npos ||
//                nm.find("SMINPUTS:5") != std::string::npos ||
//                nm.find("SMINPUTS:6") != std::string::npos) {
//         out.limits = std::make_pair(std::max(1e-12, value - 5.0 * s), value + 5.0 * s);
//     }

//     return out;
// }

// MarginalConfig StatisticManager::make_nuisance_marginal_config(const ParamId& pid,
//                                                                MarginalType mt) const
// {
//     // if (const auto* spec = find_nuisance_spec(pid)) {
//     //     return MarginalConfigFactory().create(pid, mt, *spec);
//     // }
//     return MarginalConfigFactory().create(pid, mt);
// }

// std::map<ParamId, double> StatisticManager::get_all_obss_deps() {
//     std::unordered_set<ParamId> eta_infos;

//     for (const auto& obsId : obs_int->get_obs_ids()) {
//         for (auto paramId : obs_int->get_obs_deps(obsId.s)) {
//             eta_infos.insert(paramId);
//         }
//     }

//     std::unordered_set<ParamId> eta_infos_leaf = this->sp->get_all_leaf_sources(eta_infos);
//     std::map<ParamId, double> eta_specs_real_leaf;
//     std::map<ParamId, double> delta_rel;

//     for (auto& paramId : eta_infos_leaf) {
//         double u = pspp->get_param(paramId)->get_combined_std().real();
//         double val = pspp->get_param(paramId)->get_val().real();
//         if (fpeq(val, 0.0) && !fpeq(u, 0.0)) {
//             eta_specs_real_leaf[paramId] = val;
//             continue;
//         }
            
//         if (!std::isfinite(u) || fpeq(std::abs(u), 0.0))
//             delta_rel[paramId] = 0;
//         else
//             delta_rel[paramId] = std::abs(u / val);
//     }
        
    
//     using T = std::pair<ParamId, double>;
//     double delta_rel_max = std::max_element(
//         delta_rel.begin(), delta_rel.end(), 
//         [] (const T& p, const T& q) { return p.second < q.second; }
//     )->second;

//     for (auto& paramId : eta_infos_leaf) {
//         LOG_INFO("Compared to max relative uncertainty", paramId, delta_rel[paramId] / delta_rel_max);
//     }

//     for (auto& paramId : eta_infos_leaf) {
//         if (delta_rel[paramId] / delta_rel_max > config.nuisance_relevance_cutoff)
//             eta_specs_real_leaf[paramId] = pspp->get_param(paramId)->get_val();
//     }

//     LOG_INFO("Significant nuisances");
//     for (const auto& [pid, val] : eta_specs_real_leaf) {
//         LOG_INFO(pid, val);
//     }

//     // for (auto it = eta_specs_real_leaf.begin(); it != eta_specs_real_leaf.end(); ) {
//     //     std::ostringstream oss;
//     //     oss << it->first;
//     //     if (oss.str() == "B_Ks:1_5_2") {
//     //         std::cout << "[FIT] DEBUG removing nuisance " << oss.str() << "\n";
//     //         it = eta_specs_real_leaf.erase(it);
//     //     } else {
//     //         ++it;
//     //     }
//     // }
//     return eta_specs_real_leaf;
// }

// std::map<ParamId, double> StatisticManager::get_p_specs(const std::vector<ParamId>& p_specs) {
//     std::map<ParamId, double> out;
//     for (auto elem : p_specs) {            
//         out[elem] = pspp->get_param(elem)->get_val();
//     }
//     return out;
// }
// std::map<ParamId, std::map<ParamId, double>> StatisticManager::get_all_correlations() {
//     std::map<ParamId, std::map<ParamId, double>> res;
//     CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);
//     res = ct.transform(this->cache.eta_specs_real);
//     return res;
// }

// std::map<ExperimentObs, std::map<ExperimentObs, double>> StatisticManager::get_all_obs_correlations() {
//     std::map<ExperimentObs, std::map<ExperimentObs, double>> res;
//     CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);

//     res = ct.transform(this->cache.exp_obs);
//     return res;
// }

// std::map<ExperimentObs, double> StatisticManager::get_obs_exp() {
//     std::map<ExperimentObs, double> out;

//     for (const auto& obsId : obs_int->get_obs_ids()) {
//         auto _ = pspp->get_obs_param(obsId);
//         for (auto _2 : _) {
//             out[_2.first] = _2.second->get_val(); 
//         }
//         // out[obsId] = pspp->get_obs_param(obsId)->get_val();
//     }
//     return out;
// }


#include "StatisticManager.h"

#include <limits>
#include <numeric>
#include <set>


bool starts_with(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() &&
           std::equal(prefix.begin(), prefix.end(), s.begin());
}

std::string param_name(const ParamId& pid) {
    std::ostringstream oss;
    oss << pid;
    return oss.str();
}
std::vector<std::size_t> select_indices_by_prefix(const std::vector<ParamId>& ids,
                                                  const std::string& prefix)
{
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (starts_with(param_name(ids[i]), prefix)) {
            out.push_back(i);
        }
    }
    return out;
}

std::vector<std::size_t> select_indices_by_names(const std::vector<ParamId>& ids,
                                                 const std::set<std::string>& names)
{
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (names.contains(param_name(ids[i]))) {
            out.push_back(i);
        }
    }
    return out;
}

RealMatrix principal_submatrix(const RealMatrix& M,
                               const std::vector<std::size_t>& idx)
{
    RealMatrix out(idx.size(), idx.size());
    for (std::size_t i = 0; i < idx.size(); ++i) {
        for (std::size_t j = 0; j < idx.size(); ++j) {
            out.at(i, j) = M.at(idx[i], idx[j]);
        }
    }
    return out;
}

void log_smallest_mode(const RealMatrix& M,
                       const std::vector<ParamId>& ids,
                       const std::vector<std::size_t>& idx,
                       const std::string& tag)
{
    std::cout << "[FIT] ---- " << tag << " ----\n";
    std::cout << "[FIT] size=" << idx.size() << "\n";
    if (idx.empty()) {
        std::cout << "[FIT] empty block\n";
        return;
    }

    RealMatrix S = principal_submatrix(M, idx);

    double diag_min = std::numeric_limits<double>::infinity();
    double diag_max = -std::numeric_limits<double>::infinity();
    double max_abs_offdiag = 0.0;

    for (std::size_t i = 0; i < S.rows(); ++i) {
        diag_min = std::min(diag_min, S.at(i, i));
        diag_max = std::max(diag_max, S.at(i, i));
        for (std::size_t j = i + 1; j < S.cols(); ++j) {
            max_abs_offdiag = std::max(max_abs_offdiag, std::abs(S.at(i, j)));
        }
    }

    std::cout << "[FIT] diag_min=" << diag_min
              << " diag_max=" << diag_max
              << " max_abs_offdiag=" << max_abs_offdiag << "\n";

    EigenSystem e = S.eig();

    std::size_t imin = 0;
    double evmin = e.D.at(0, 0);
    for (std::size_t i = 1; i < S.rows(); ++i) {
        if (e.D.at(i, i) < evmin) {
            evmin = e.D.at(i, i);
            imin = i;
        }
    }

    std::cout << "[FIT] smallest_eigenvalue=" << evmin << "\n";

    std::vector<std::pair<double, std::size_t>> comps;
    for (std::size_t i = 0; i < idx.size(); ++i) {
        comps.push_back({std::abs(e.P.at(i, imin)), i});
    }

    std::sort(comps.begin(), comps.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    std::cout << "[FIT] dominant entries of smallest mode:\n";
    for (std::size_t k = 0; k < std::min<std::size_t>(10, comps.size()); ++k) {
        const std::size_t local_i = comps[k].second;
        const std::size_t global_i = idx[local_i];
        std::cout << "[FIT]   " << ids[global_i]
                  << " coeff=" << e.P.at(local_i, imin) << "\n";
    }

    std::cout << "[FIT] matrix entries:\n";
    for (std::size_t i = 0; i < idx.size(); ++i) {
        for (std::size_t j = 0; j < idx.size(); ++j) {
            std::cout << "[FIT]   "
                      << ids[idx[i]] << " | " << ids[idx[j]]
                      << " = " << S.at(i, j) << "\n";
        }
    }
}


static void dump_matrix_sanity(const RealMatrix& M,
                               const std::vector<ParamId>& ids,
                               const std::string& name)
{
    std::cout << "[FIT] Matrix " << name
              << " shape=" << M.rows() << "x" << M.cols()
              << " symmetric=" << M.is_symmetric() << "\n";

    double diag_min = std::numeric_limits<double>::infinity();
    double diag_max = -std::numeric_limits<double>::infinity();
    std::size_t nonpos_diag = 0;

    double max_abs_offdiag = 0.0;
    std::size_t n_abs_gt_1 = 0;
    std::pair<std::size_t,std::size_t> argmax_offdiag{0,0};

    for (std::size_t i = 0; i < M.rows(); ++i) {
        const double d = M.at(i,i);
        diag_min = std::min(diag_min, d);
        diag_max = std::max(diag_max, d);
        if (!(d > 0.0)) ++nonpos_diag;

        for (std::size_t j = i + 1; j < M.cols(); ++j) {
            const double a = std::abs(M.at(i,j));
            if (a > max_abs_offdiag) {
                max_abs_offdiag = a;
                argmax_offdiag = {i,j};
            }
            if (a > 1.0 + 1e-10) ++n_abs_gt_1;
        }
    }

    std::cout << "[FIT] " << name
              << " diag_min=" << diag_min
              << " diag_max=" << diag_max
              << " nonpos_diag=" << nonpos_diag
              << " max_abs_offdiag=" << max_abs_offdiag
              << " n_abs_offdiag_gt_1=" << n_abs_gt_1
              << "\n";

    if (M.rows() == M.cols() && M.is_symmetric()) {
        EigenSystem e = M.eig();

        std::size_t imin = 0;
        double evmin = e.D.at(0,0);
        for (std::size_t i = 1; i < M.rows(); ++i) {
            if (e.D.at(i,i) < evmin) {
                evmin = e.D.at(i,i);
                imin = i;
            }
        }

        std::vector<std::pair<double,std::size_t>> comps;
        for (std::size_t i = 0; i < M.rows(); ++i) {
            comps.push_back({std::abs(e.P.at(i, imin)), i});
        }
        std::sort(comps.begin(), comps.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        std::cout << "[FIT] " << name << " smallest_eigenvalue=" << evmin << "\n";
        std::cout << "[FIT] " << name << " dominant entries of smallest mode:\n";
        for (std::size_t k = 0; k < std::min<std::size_t>(10, comps.size()); ++k) {
            std::size_t i = comps[k].second;
            std::cout << "[FIT]   " << ids[i]
                      << " coeff=" << e.P.at(i, imin) << "\n";
        }
    }
}

namespace {

// std::string param_name(const ParamId& pid) {
//     std::ostringstream oss;
//     oss << pid;
//     return oss.str();
// }

fit_app::ParameterDefinition make_fit_param_def(const ParamId& pid, double value, double sigma_hint) {
    fit_app::ParameterDefinition out;
    out.name = param_name(pid);
    out.value = value;
    out.step_hint = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
        ? sigma_hint
        : std::max(1e-3, 0.01 * std::abs(value));

    const std::string& nm = out.name;

    if (nm.find("FCONST") != std::string::npos) {
        out.limits = std::make_pair(0.05, 0.35);
    }

    return out;
}

fit_app::ParameterDefinition make_nuisance_param_def(const ParamId& pid, double value, double sigma_hint) {
    fit_app::ParameterDefinition out;
    out.name = param_name(pid);
    out.value = value;

    const double s = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
        ? std::abs(sigma_hint)
        : std::max(1e-3, 0.01 * std::abs(value));

    out.step_hint = s;

    const std::string& nm = out.name;

    // TODO : Very ugly
    if (nm.find("SMINPUTS:3") != std::string::npos) {
        out.limits = std::make_pair(0.05, 0.30);
    } else if (nm.find("MASS:") != std::string::npos ||
            //    nm.find("FLIFE:") != std::string::npos ||
               nm.find("FCONST:") != std::string::npos ||
               nm.find("FMASS:") != std::string::npos ||
               nm.find("SMINPUTS:5") != std::string::npos ||
               nm.find("SMINPUTS:6") != std::string::npos) {
        out.limits = std::make_pair(std::max(1e-12, value - 5.0 * s), value + 5.0 * s);
    }

    return out;
}

template <class PredMapT>
std::vector<double> ordered_prediction_vector(
    const std::vector<ExperimentObs>& obs_ids,
    const PredMapT& pred_map)
{
    std::vector<double> out;
    out.reserve(obs_ids.size());

    for (const auto& bid : obs_ids) {
        const auto& vec = pred_map.at(bid.obs.s);

        auto it = std::find_if(vec.begin(), vec.end(), [&](const auto& ov) {
            auto bin = ov.bin.value_or(std::pair<double, double>{0.0, 0.0});
            return bin == bid.obs.p;
        });

        if (it == vec.end()) {
            throw std::runtime_error("Missing predicted observable/bin.");
        }

        out.push_back(it->value);
    }

    return out;
}


void log_selected_nuisance_diagnostics(const std::vector<ParamId>& eta_ids,
                                      const std::map<ParamId, std::map<ParamId, double>>& sigma_eta) {
    std::cout << "[FIT] Selected nuisance parameters: " << eta_ids.size() << std::endl;
    for (const auto& pid : eta_ids) {
        std::cout << "[FIT]   nuisance " << pid << std::endl;
    }

    if (sigma_eta.empty()) {
        std::cout << "[FIT] SigmaEta is empty" << std::endl;
        return;
    }

    try {
        RealMatrix sigma = RealMatrix(unzip(sigma_eta).vals);
        std::cout << "[FIT] SigmaEta shape=" << sigma.rows() << "x" << sigma.cols()
                  << ", symmetric=" << sigma.is_symmetric() << std::endl;

        RealMatrix sym = 0.5 * (sigma + sigma.transpose());
        EigenSystem eig = sym.eig();

        double min_eig = std::numeric_limits<double>::infinity();
        double max_eig = -std::numeric_limits<double>::infinity();
        std::size_t non_pos = 0;
        std::size_t tiny = 0;

        for (std::size_t i = 0; i < eig.D.rows(); ++i) {
            const double lambda = eig.D.at(i, i);
            min_eig = std::min(min_eig, lambda);
            max_eig = std::max(max_eig, lambda);
            if (!(lambda > 0.0)) ++non_pos;
            if (std::abs(lambda) < 1e-12) ++tiny;
        }

        std::cout << "[FIT] SigmaEta eig_min=" << min_eig
                  << ", eig_max=" << max_eig
                  << ", non_positive=" << non_pos
                  << ", tiny(|eig|<1e-12)=" << tiny
                  << std::endl;
    } catch (const std::exception& e) {
        std::cout << "[FIT] SigmaEta diagnostics failed: " << e.what() << std::endl;
    }
}

}

StatisticManager::StatisticManager(StatisticConfig config,
                                   std::shared_ptr<IModel> obs_int,
                                   std::shared_ptr<IStatCorrelationProxy> pscp,
                                   std::shared_ptr<IStatParameterProxy> pspp,
                                   std::shared_ptr<IStatSourcesProxy> sp,
                                   std::shared_ptr<IStatDependencyPruner> dp,
                                   std::shared_ptr<INuisanceReader> nuisance_reader,
                                   std::shared_ptr<IStatParamOptimizerProxy> spop)
    : obs_int(std::move(obs_int)),
      pscp(std::move(pscp)),
      pspp(std::move(pspp)),
      sp(std::move(sp)),
      dp(std::move(dp)),
      nuisance_reader_(std::move(nuisance_reader)),
      spop(std::move(spop)),
      config(std::move(config))
{
    if (!nuisance_reader_) {
        throw std::invalid_argument("StatisticManager: nuisance_reader is null");
    }

    this->obs_int->compute_observables();
    reload_nuisance_specs();
    invalidate_fit_state();
}

// std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
//     // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
//     unsigned int seed = 123456u;
//     std::map<ParamId, MarginalType> nuisance_marginals;

//     for (auto& [pid, v] : cache.eta_specs_real) {
//         nuisance_marginals.emplace(pid, MarginalType::GAUSSIAN);
//     }

//     for (auto& [pid, mt] : config.override_nuisance_marginals) {
//         for (auto elem : nuisance_marginals) {
//             std::cout << elem.first << std::endl;
//             std::cout << ParameterTypeMapper::str(elem.first.type.value_or(ParameterType::WILSON)) << std::endl;
//         }
//         if (!nuisance_marginals.contains(pid)) {
//             LOG_WARN("Parameter", pid, "is not a nuisance for the selected observables");
//             continue;
//         }
//         nuisance_marginals.at(pid) = mt;
//     }

//     auto unzipped = unzip(nuisance_marginals);
//     std::vector<ParamId> nuisance_ids = unzipped.ids;
//     std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

//     for (auto& [pid, mt] : nuisance_marginals) {
//         MarginalConfig cfg = MarginalConfigFactory().create(pid, mt);
//         auto m_ptr = DistributionFactory::create(mt, cfg, seed);
//         marginals.emplace_back(std::move(m_ptr));
//     }

//     return marginals;
// }

#include "BlockProxy.h"
std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
    unsigned int seed = 123456u;
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

    marginals.reserve(cache.eta_specs_real.size());

    for (const auto& [pid, _] : cache.eta_specs_real) {
        const MarginalType mt = resolve_nuisance_marginal_type(pid);
        MarginalConfig cfg = make_nuisance_marginal_config(pid, mt);
        std::ostringstream oss;
        oss << pid;
        const std::string name = oss.str();
        if (name == "B_SCALE:1" || name == "EW_SCALE:1") {
            BlockProxy().log_block(ParameterType::WILSON, "B_SCALE");
            std::visit([&](const auto& c) {
                using T = std::decay_t<decltype(c)>;
                if constexpr (std::is_same_v<T, FlatMarginalCfg>) {
                    std::cout << "[MC CFG] " << pid
                            << " FLAT a=" << c.a
                            << " b=" << c.b << '\n';
                } else if constexpr (std::is_same_v<T, GaussianMarginalCfg>) {
                    std::cout << "[MC CFG] " << pid
                            << " GAUSSIAN mu=" << c.mu
                            << " sigma=" << c.sigma << '\n';
                } else if constexpr (std::is_same_v<T, SplitGaussianMarginalCfg>) {
                    std::cout << "[MC CFG] " << pid
                            << " SPLIT mu=" << c.mu
                            << " sigma_p=" << c.sigma_p
                            << " sigma_m=" << c.sigma_m << '\n';
                } else {
                    std::cout << "[MC CFG] " << pid << " other cfg\n";
                }
            }, cfg);
        }

        marginals.emplace_back(MarginalFactory::create(mt, cfg, seed));
    }

    return marginals;
}

std::unique_ptr<JointDistribution> StatisticManager::build_nuisance_distribution() {
    // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
    unsigned int seed = 123456u;

    std::unique_ptr<ICopula> copula;
    if (config.nuisance_copula_type == CopulaType::GAUSSIAN) {
        GaussianCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
        copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
    } else if (config.nuisance_copula_type == CopulaType::STUDENT_T) {
        StudentTCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
        copula_cfg.nu = cache.eta_specs_real.size();
        copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
    }

    return std::make_unique<JointDistribution>(
        std::move(build_nuisance_marginal_distributions()), 
        std::move(copula)
    );
}

std::unique_ptr<JointDistribution> StatisticManager::build_exp_data_distribution() {
    unsigned int seed = std::random_device{}();
    std::map<ExperimentObs, MarginalType> exp_data_marginals;

    for (auto& [oid, v] : cache.exp_obs) {
        exp_data_marginals.emplace(oid, MarginalType::GAUSSIAN);
    }

    for (auto& [oid, mt] : config.override_exp_data_marginals) {
        if (!exp_data_marginals.contains(oid)) {
            LOG_WARN("Observable", ObservableMapper::str(oid.obs.s), "is not a taken into account in this fit.");
            continue;
        }
        exp_data_marginals.at(oid) = mt;
    }

    auto unzipped = unzip(exp_data_marginals);
    std::vector<ExperimentObs> obs_ids = unzipped.ids;
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

    for (auto& [oid, mt] : exp_data_marginals) {//TODO : checkkkkkk
        MarginalConfig cfg = MarginalConfigFactory().create(oid, mt);
        auto m_ptr = MarginalFactory::create(mt, cfg, seed);
        marginals.emplace_back(std::move(m_ptr));
    }

    std::unique_ptr<ICopula> copula;
    if (config.exp_data_copula_type == CopulaType::GAUSSIAN) {
        GaussianCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
        copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
    } else if (config.exp_data_copula_type == CopulaType::STUDENT_T) {
        StudentTCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
        copula_cfg.nu = obs_int->n_observables() - 1;
        copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
    }

    return std::make_unique<JointDistribution>(std::move(marginals), std::move(copula));
}

static void scan_leave_one_out(const RealMatrix& M,
                               const std::vector<ParamId>& ids,
                               const std::vector<std::size_t>& idx,
                               const std::string& tag)
{
    std::cout << "[FIT] ---- leave-one-out " << tag << " ----\n";
    if (idx.size() <= 1) return;

    for (std::size_t k = 0; k < idx.size(); ++k) {
        std::vector<std::size_t> sub;
        sub.reserve(idx.size() - 1);
        for (std::size_t i = 0; i < idx.size(); ++i) {
            if (i != k) sub.push_back(idx[i]);
        }

        RealMatrix S = principal_submatrix(M, sub);
        EigenSystem e = S.eig();

        double evmin = e.D.at(0,0);
        for (std::size_t i = 1; i < S.rows(); ++i) {
            evmin = std::min(evmin, e.D.at(i,i));
        }

        std::cout << "[FIT] remove " << ids[idx[k]]
                  << " => smallest_eigenvalue = " << evmin << "\n";
    }
}

FitResultWithMaps StatisticManager::compute_MLE(const std::vector<ParamId>& p_specs) {
    update_cache(p_specs);

    if (cache.p_specs.empty()) {
        throw std::invalid_argument("compute_MLE called with an empty fit parameter list.");
    }

    auto unzipped_fit_params = unzip(cache.p_specs);
    auto unzipped_nuisances  = unzip(cache.eta_specs_real);
    auto unzipped_exp_obs    = unzip(cache.exp_obs);

    {auto eta_ids = unzip(cache.eta_specs_real).ids;
    RealMatrix Reta(unzip(cache.SigmaEta).vals);
    dump_matrix_sanity(Reta, eta_ids, "SigmaEta");

}

    const std::vector<ParamId> p_ids = unzipped_fit_params.ids;
    const std::vector<double> p0 = unzipped_fit_params.vals;

    const std::vector<ParamId> eta_ids = unzipped_nuisances.ids;
    const std::vector<double> eta0 = unzipped_nuisances.vals;
    
        RealMatrix SigmaEta(unzip(cache.SigmaEta).vals);

    // bloc complet
    // {
    //     std::vector<std::size_t> all_idx(eta_ids.size());
    //     std::iota(all_idx.begin(), all_idx.end(), 0);
    //     log_smallest_mode(SigmaEta, eta_ids, all_idx, "SigmaEta full");
    // }

    // bloc B_Ks:1_*
    // {
    //     auto idx = select_indices_by_prefix(eta_ids, "B_Ks:1_");
    //     log_smallest_mode(SigmaEta, eta_ids, idx, "SigmaEta B_Ks:1_*");
    // }

    // bloc B_Ks:18_*
    // {
    //     auto idx = select_indices_by_prefix(eta_ids, "B_Ks:18_");
    //     log_smallest_mode(SigmaEta, eta_ids, idx, "SigmaEta B_Ks:18_*");
    // }

    // bloc mixte 1_* + {7_*, 8_*, 14}
    // {
    //     std::set<std::string> names = {
    //         "B_Ks:7_1", "B_Ks:7_2",
    //         "B_Ks:8_1", "B_Ks:8_2",
    //         "B_Ks:14"
    //     };

    //     auto idx1 = select_indices_by_prefix(eta_ids, "B_Ks:1_");
    //     auto idx2 = select_indices_by_names(eta_ids, names);

    //     std::vector<std::size_t> idx = idx1;
    //     idx.insert(idx.end(), idx2.begin(), idx2.end());
    //     std::sort(idx.begin(), idx.end());
    //     idx.erase(std::unique(idx.begin(), idx.end()), idx.end());

    //     log_smallest_mode(SigmaEta, eta_ids, idx, "SigmaEta B_Ks:1_* + {7,8,14}");
    // }

    // {
    //     std::set<std::string> names = {
    //         "B_Ks:1_2_1",
    //         "B_Ks:1_5_1",
    //         "B_Ks:1_5_2",
    //         "B_Ks:1_6_2",
    //         "B_Ks:1_2_2",
    //         "B_Ks:1_4_1",
    //         "B_Ks:1_3_1",
    //         "B_Ks:1_1_1",
    //         "B_Ks:1_4_2",
    //         "B_Ks:1_3_2"
    //     };
    //     auto idx = select_indices_by_names(eta_ids, names);
    //     log_smallest_mode(SigmaEta, eta_ids, idx, "SigmaEta suspects");
    //     scan_leave_one_out(SigmaEta, eta_ids, idx, "SigmaEta suspects");
    // }

    // log_selected_nuisance_diagnostics(eta_ids, cache.SigmaEta);

    const std::vector<ExperimentObs> obs_ids = unzipped_exp_obs.ids;
    const std::vector<double> exp_obs_vals = unzipped_exp_obs.vals;

    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->nuisance_dist = build_nuisance_distribution();
    ctx->exp_obs_dist = build_exp_data_distribution();
    ctx->exp_obs_values = exp_obs_vals;

    // ctx->fp_defs.reserve(p_ids.size());
    // for (std::size_t i = 0; i < p_ids.size(); ++i) {
    //     double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
    //     ctx->fp_defs.emplace_back(make_param_def(p_ids[i], p0[i], sigma));
    // }

    // ctx->nuis_defs.reserve(eta_ids.size());
    // for (std::size_t i = 0; i < eta_ids.size(); ++i) {
    //     double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
    //     ctx->nuis_defs.emplace_back(make_param_def(eta_ids[i], eta0[i], sigma));
    // }

    ctx->fp_defs.reserve(p_ids.size());
    for (std::size_t i = 0; i < p_ids.size(); ++i) {
        double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
        ctx->fp_defs.emplace_back(make_fit_param_def(p_ids[i], p0[i], sigma));
    }

    // ctx->nuis_defs.reserve(eta_ids.size());
    // for (std::size_t i = 0; i < eta_ids.size(); ++i) {
    //     double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
    //     ctx->nuis_defs.emplace_back(make_nuisance_param_def(eta_ids[i], eta0[i], sigma));
    // }
    ctx->nuis_defs.reserve(eta_ids.size());
    for (std::size_t i = 0; i < eta_ids.size(); ++i) {
        double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
        ctx->nuis_defs.emplace_back(
            make_nuisance_parameter_definition(eta_ids[i], eta0[i], sigma)
        );
    }

    auto model_fn = [this, obs_ids, p_ids, eta_ids](
        const std::vector<double>& p_vec,
        const std::vector<double>& eta_vec) -> std::vector<double>
    {
        auto pred_map = this->obs_int->predict_optimized(
            zip(p_ids, p_vec),
            zip(eta_ids, eta_vec)
        );

        return ordered_prediction_vector(obs_ids, pred_map);
    };

    // last_ctx_ = ctx;
    // last_like_ = std::make_shared<BaseLikelihood>(model_fn, ctx, p_ids.size());
    // last_fitter_ = std::make_shared<MLFitter>(ctx, model_fn);

    // last_fit_raw_ = last_fitter_->maximum_likelihood_fit(p0);

    last_ctx_ = ctx;
    last_like_ = std::make_shared<BaseLikelihood>(model_fn, ctx, p_ids.size());

    MLFitOptions fitopt;
    fitopt.run_hesse = config.MLE_run_hesse;
    fitopt.request_minos = config.MLE_request_minos;
    fitopt.verbose = config.MLE_verbose;
    fitopt.strategy = config.MLE_strategy;
    fitopt.max_fcn = static_cast<unsigned>(config.MLE_max_iter);
    fitopt.tolerance = config.MLE_tol;

    fitopt.allow_profile_hessian_fallback = config.MLE_allow_profile_hessian_fallback;
    fitopt.profile_hessian_step_scale = config.MLE_profile_hessian_step_scale;
    fitopt.profile_hessian_eig_floor_rel = config.MLE_profile_hessian_eig_floor_rel;

    fitopt.trace_first_evals = config.MLE_trace_first_evals;
    fitopt.trace_max_evals = config.MLE_trace_max_evals;

    last_fitter_ = std::make_shared<MLFitter>(ctx, model_fn, fitopt);
    last_fit_raw_ = last_fitter_->maximum_likelihood_fit(p0);

    last_fit_param_ids_ = p_ids;
    last_nuisance_ids_ = eta_ids;
    last_fit_param_index_.clear();
    for (std::size_t i = 0; i < p_ids.size(); ++i) {
        last_fit_param_index_[p_ids[i]] = i;
    }

    FitResultWithMaps out;
    out.fit_ok = !last_fit_raw_.p_hat.empty();
    out.ell_hat = last_fit_raw_.ell_hat;
    out.p_hat = zip(p_ids, last_fit_raw_.p_hat);
    out.eta_hat = zip(eta_ids, last_fit_raw_.eta_hat);
    out.p_hat_std = zip(p_ids, last_fit_raw_.p_hat_std);
    out.p_correlations = zip(p_ids, last_fit_raw_.p_hat_correlations);

    cache.mle_result = out;
    return out;
}

Contour StatisticManager::confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, ContourOptions options) {
    if (!cache.mle_result.fit_ok || !last_fitter_) {
        throw std::runtime_error("Please run compute_MLE before requesting a confidence contour.");
    }

    if (!last_fit_param_index_.contains(p1) || !last_fit_param_index_.contains(p2)) {
        throw std::invalid_argument("Contour requested for parameters that are not in the last fitted parameter set.");
    }

    if (p1 == p2) {
        throw std::invalid_argument("Contour requires two distinct parameters.");
    }

    const std::size_t x_id = last_fit_param_index_.at(p1);
    const std::size_t y_id = last_fit_param_index_.at(p2);

    Contour cl = last_fitter_->contour(
        x_id,
        y_id,
        z,
        bounds,
        options
    );

    return cl;
}

void StatisticManager::print_cache()
{
    for (auto elem : cache.eta_specs_real) {
        std::cout << " eta_specs_real : " << elem.first << " = " << elem.second;
    }
    std::cout << std::endl;

    for (auto elem : cache.SigmaEta) {
        for (auto elem2 : elem.second) {
            std::cout << " SigmaEta : " << elem.first << " | " << elem2.first << " = " << elem2.second;
        }
        std::cout << std::endl;
    }

    for (auto elem : cache.exp_obs) {
        std::cout << " exp_obs : " << elem.first.str() << " = " << elem.second;
    }
    std::cout << std::endl;

    for (auto elem : cache.SigmaObs) {
        for (auto elem2 : elem.second) {
            std::cout << " SigmaObs : " << elem.first.str() << " | " << elem2.first.str() << " = " << elem2.second;
        }
        std::cout << std::endl;
    }

    for (auto elem : cache.p_specs) {
        std::cout << " p_specs : " << elem.first << " = " << elem.second;
    }
    std::cout << std::endl;

    std::cout << "eta size : " << this->cache.eta_specs_real.size() << std::endl;
    std::cout << "etasigma size : " << this->cache.eta_specs_real.size() << " | " << this->cache.SigmaEta.begin()->second.size() << std::endl;
    std::cout << "p_specs size : " << this->cache.p_specs.size() << std::endl;
    std::cout << "exp_obs : " << this->cache.exp_obs.size() << std::endl;
    std::cout << "SigmaObs size : " << this->cache.SigmaObs.size() << " | " << this->cache.SigmaObs.begin()->second.size() << std::endl;
}

void StatisticManager::reload_nuisance_specs() {
    default_nuisance_specs_ = nuisance_reader_->load_default();

    if (current_user_nuisance_file_.has_value()) {
        user_nuisance_specs_ = nuisance_reader_->load_user(*current_user_nuisance_file_);
    } else {
        user_nuisance_specs_ = nuisance_reader_->load_user();
    }

    rebuild_merged_nuisance_specs();
    invalidate_fit_state();
}

void StatisticManager::set_nuisance_user_file(const fs::path& user_yaml_path) {
    current_user_nuisance_file_ = user_yaml_path;
    reload_nuisance_specs();
}

void StatisticManager::clear_nuisance_user_file() {
    current_user_nuisance_file_.reset();
    reload_nuisance_specs();
}

void StatisticManager::rebuild_merged_nuisance_specs() {
    merged_nuisance_specs_ = default_nuisance_specs_;
    for (const auto& [pid, spec] : user_nuisance_specs_) {
        merged_nuisance_specs_[pid] = spec;
    }

}

void StatisticManager::invalidate_fit_state() {
    cache.mle_result = FitResultWithMaps{};

    last_ctx_.reset();
    last_like_.reset();
    last_fitter_.reset();

    last_fit_param_ids_.clear();
    last_nuisance_ids_.clear();
    last_fit_param_index_.clear();
}

const NuisanceSpec* StatisticManager::find_nuisance_spec(const ParamId& pid) const {
    auto it = merged_nuisance_specs_.find(pid);
    if (it == merged_nuisance_specs_.end()) {
        return nullptr;
    }
    return &it->second;
}

MarginalType StatisticManager::resolve_nuisance_marginal_type(const ParamId& pid) const {
    MarginalType mt = MarginalType::GAUSSIAN;

    if (const auto* spec = find_nuisance_spec(pid)) {
        mt = spec->marginal;
    }

    if (config.override_nuisance_marginals.contains(pid)) {
        mt = config.override_nuisance_marginals.at(pid);
    }

    return mt;
}

fit_app::ParameterDefinition StatisticManager::make_nuisance_parameter_definition(const ParamId& pid,
                                                                                  double value,
                                                                                  double sigma_hint) const
{
    fit_app::ParameterDefinition out;
    out.name = param_name(pid);
    out.value = value;

    const double s = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
        ? std::abs(sigma_hint)
        : std::max(1e-3, 0.01 * std::abs(value));

    out.step_hint = s;

    if (const auto* spec = find_nuisance_spec(pid)) {
        out.limits = spec->bounds;
        return out;
    }

    const std::string& nm = out.name;
    if (nm.find("SMINPUTS:3") != std::string::npos) {
        out.limits = std::make_pair(0.05, 0.30);
    } else if (nm.find("MASS:") != std::string::npos ||
               nm.find("FLIFE:") != std::string::npos ||
               nm.find("FCONST:") != std::string::npos ||
               nm.find("FMASS:") != std::string::npos ||
               nm.find("SMINPUTS:5") != std::string::npos ||
               nm.find("SMINPUTS:6") != std::string::npos) {
        out.limits = std::make_pair(std::max(1e-12, value - 5.0 * s), value + 5.0 * s);
    }

    return out;
}

MarginalConfig StatisticManager::make_nuisance_marginal_config(const ParamId& pid,
                                                               MarginalType mt) const
{
    // if (const auto* spec = find_nuisance_spec(pid)) {
    //     return MarginalConfigFactory().create(pid, mt, *spec);
    // }
    return MarginalConfigFactory().create(pid, mt);
}

// MarginalConfig StatisticManager::make_nuisance_marginal_config(const ParamId& pid,
//                                                                MarginalType mt) const
// {
//     if (const auto* spec = find_nuisance_spec(pid)) {
//         NuisanceSpec effective_spec = *spec;
//         effective_spec.marginal = mt; // respecte aussi les overrides éventuels
//         return MarginalConfigFactory().create(pid, mt, effective_spec);
//     }

//     return MarginalConfigFactory().create(pid, mt);
// }

std::map<ParamId, double> StatisticManager::get_all_obss_deps() {
    std::unordered_set<ParamId> eta_infos;

    for (const auto& obsId : obs_int->get_obs_ids()) {
        for (auto paramId : obs_int->get_obs_deps(obsId.s)) {
            eta_infos.insert(paramId);
        }
    }

    std::unordered_set<ParamId> eta_infos_leaf = this->sp->get_all_leaf_sources(eta_infos);

    std::map<ParamId, double> eta_specs_real_leaf;
    std::map<ParamId, double> delta_rel;

    // 1) Préfiltre existant : relevance basée sur l'incertitude relative
    for (const auto& paramId : eta_infos_leaf) {
        // On ne traite jamais les paramètres de fit comme nuisances
        if (cache.p_specs.contains(paramId)) {
            continue;
        }

        const double u   = std::abs(pspp->get_param(paramId)->get_combined_std().real());
        const double val = pspp->get_param(paramId)->get_val().real();

        if (!std::isfinite(u) || fpeq(u, 0.0)) {
            delta_rel[paramId] = 0.0;
            continue;
        }

        if (fpeq(val, 0.0)) {
            // Cas utile : valeur centrale nulle mais sigma non nulle.
            // On le garde au moins jusqu'au filtre de sensibilité.
            eta_specs_real_leaf[paramId] = val;
            delta_rel[paramId] = std::numeric_limits<double>::infinity();
            continue;
        }

        delta_rel[paramId] = std::abs(u / val);
    }

    double delta_rel_max = 0.0;
    for (const auto& [pid, d] : delta_rel) {
        if (std::isfinite(d)) {
            delta_rel_max = std::max(delta_rel_max, d);
        }
    }
    if (!(delta_rel_max > 0.0)) {
        delta_rel_max = 1.0;
    }

    for (const auto& [pid, d] : delta_rel) {
        const double rel_to_max = std::isfinite(d) ? (d / delta_rel_max) : 1.0;
        LOG_INFO("Compared to max relative uncertainty", pid, rel_to_max);
    }

    for (const auto& [pid, d] : delta_rel) {
        const double rel_to_max = std::isfinite(d) ? (d / delta_rel_max) : 1.0;
        if (rel_to_max > config.nuisance_relevance_cutoff || !std::isfinite(d)) {
            eta_specs_real_leaf[pid] = pspp->get_param(pid)->get_val();
        }
    }

    // 2) Nouveau filtre : sensibilité locale du modèle
    //
    // Important :
    // - ce n'est pas un argument "mathématique" de platitude stricte,
    //   car le prior des nuisances reste dans la NLL ;
    // - c'est un argument pratique/numerique : si une nuisance ne bouge
    //   pratiquement pas les observables, elle alourdit le fit pour peu de gain.
    //
    // On le désactive si :
    // - demandé explicitement,
    // - aucun nuisance retenu,
    // - ou pas de paramètres de fit courants disponibles.
    if (config.nuisance_sensitivity_pruning &&
        !eta_specs_real_leaf.empty() &&
        !cache.p_specs.empty())
    {
        // TODO : Makes fit unstable.
        // std::unordered_set<ParamId> shifted {};
        // // Shifting every zero nuisance to avoid false negatives
        // for (auto& [pid, val]: eta_specs_real_leaf) {
        //     if (fpeq(val, 0.0)) {
        //         LOG_INFO("Shifting", pid);
        //         double u_1sigma = std::abs(pspp->get_param(pid)->get_combined_std().real());
        //         eta_specs_real_leaf.at(pid) = u_1sigma / 2;
        //         shifted.emplace(pid);
        //     }
        // }

        const auto exp_obs_map = get_obs_exp();
        const auto unzipped_exp_obs = unzip(exp_obs_map);
        const std::vector<ExperimentObs> obs_ids = unzipped_exp_obs.ids;
        const std::vector<double> exp_obs_vals   = unzipped_exp_obs.vals;

        const auto baseline_pred_map =
            this->obs_int->predict_optimized(cache.p_specs, eta_specs_real_leaf);
        const std::vector<double> baseline_pred =
            ordered_prediction_vector(obs_ids, baseline_pred_map);

        std::map<ParamId, double> screened_eta_specs;

        std::cout << "[FIT] Model-sensitivity pruning on "
                  << eta_specs_real_leaf.size()
                  << " nuisance candidates" << std::endl;

        for (const auto& [pid, nominal] : eta_specs_real_leaf) {
            const double sigma =
                std::abs(pspp->get_param(pid)->get_combined_std().real());

            // Si sigma est nul / non fini, on garde par sécurité.
            if (!std::isfinite(sigma) || fpeq(sigma, 0.0)) {
                screened_eta_specs[pid] = nominal;
                std::cout << "[FIT] sensitivity " << pid
                          << " : sigma not finite or zero -> keep" << std::endl;
                continue;
            }

            const fit_app::ParameterDefinition def =
                make_nuisance_parameter_definition(pid, nominal, sigma);

            double step = config.nuisance_sensitivity_probe_sigmas * sigma;
            if (!std::isfinite(step) || fpeq(step, 0.0)) {
                screened_eta_specs[pid] = nominal;
                std::cout << "[FIT] sensitivity " << pid
                          << " : probe step not finite or zero -> keep" << std::endl;
                continue;
            }

            double eta_plus  = nominal + step;
            double eta_minus = nominal - step;

            if (def.limits.has_value()) {
                const auto [lo, hi] = def.limits.value();
                eta_plus  = std::clamp(eta_plus,  lo, hi);
                eta_minus = std::clamp(eta_minus, lo, hi);
            }

            if (!std::isfinite(eta_plus) ||
                !std::isfinite(eta_minus) ||
                std::abs(eta_plus - eta_minus) < 1e-14)
            {
                screened_eta_specs[pid] = nominal;
                std::cout << "[FIT] sensitivity " << pid
                          << " : effective probe collapsed -> keep" << std::endl;
                continue;
            }

            auto eta_plus_map  = eta_specs_real_leaf;
            auto eta_minus_map = eta_specs_real_leaf;
            eta_plus_map[pid]  = eta_plus;
            eta_minus_map[pid] = eta_minus;

            const auto pred_plus_map =
                this->obs_int->predict_optimized(cache.p_specs, eta_plus_map);
            const auto pred_minus_map =
                this->obs_int->predict_optimized(cache.p_specs, eta_minus_map);

            const std::vector<double> pred_plus =
                ordered_prediction_vector(obs_ids, pred_plus_map);
            const std::vector<double> pred_minus =
                ordered_prediction_vector(obs_ids, pred_minus_map);

            double max_abs_shift = 0.0;
            double max_rel_shift = 0.0;

            for (std::size_t i = 0; i < baseline_pred.size(); ++i) {
                // Approxime le shift 1σ par demi-différence symétrique
                const double one_sigma_shift =
                    0.5 * std::abs(pred_plus[i] - pred_minus[i]);

                const double scale = std::max({
                    std::abs(baseline_pred[i]),
                    std::abs(exp_obs_vals[i]),
                    config.nuisance_sensitivity_scale_floor
                });

                max_abs_shift = std::max(max_abs_shift, one_sigma_shift);
                max_rel_shift = std::max(max_rel_shift, one_sigma_shift / scale);
            }

            const bool keep =
                (max_abs_shift >= config.nuisance_sensitivity_abs_cutoff) ||
                (max_rel_shift >= config.nuisance_sensitivity_rel_cutoff);

            std::cout << "[FIT] sensitivity " << pid
                      << " : max_abs_shift=" << max_abs_shift
                      << ", max_rel_shift=" << max_rel_shift
                      << " -> " << (keep ? "keep" : "drop")
                      << std::endl;

            if (keep) {
                screened_eta_specs[pid] = nominal;
            }
        }

        // Shifting back to zero
        // for (ParamId pid : shifted)
        //     screened_eta_specs.at(pid) = 0.0;

        eta_specs_real_leaf = std::move(screened_eta_specs);
    }

    LOG_INFO("Significant nuisances");
    for (const auto& [pid, val] : eta_specs_real_leaf) {
        LOG_INFO(pid, val);
    }

    return eta_specs_real_leaf;
}

// std::map<ParamId, double> StatisticManager::get_all_obss_deps() {
//     std::unordered_set<ParamId> eta_infos;

//     for (const auto& obsId : obs_int->get_obs_ids()) {
//         for (auto paramId : obs_int->get_obs_deps(obsId.s)) {
//             eta_infos.insert(paramId);
//         }
//     }

//     std::unordered_set<ParamId> eta_infos_leaf = this->sp->get_all_leaf_sources(eta_infos);
//     std::map<ParamId, double> eta_specs_real_leaf;
//     std::map<ParamId, double> delta_rel;

//     for (auto& paramId : eta_infos_leaf) {
//         double u = pspp->get_param(paramId)->get_combined_std().real();
//         double val = pspp->get_param(paramId)->get_val().real();
//         if (fpeq(val, 0.0) && !fpeq(u, 0.0)) {
//             eta_specs_real_leaf[paramId] = val;
//             continue;
//         }
            
//         if (!std::isfinite(u) || fpeq(std::abs(u), 0.0))
//             delta_rel[paramId] = 0;
//         else
//             delta_rel[paramId] = std::abs(u / val);
//     }
        
    
//     using T = std::pair<ParamId, double>;
//     double delta_rel_max = std::max_element(
//         delta_rel.begin(), delta_rel.end(), 
//         [] (const T& p, const T& q) { return p.second < q.second; }
//     )->second;

//     for (auto& paramId : eta_infos_leaf) {
//         LOG_INFO("Compared to max relative uncertainty", paramId, delta_rel[paramId] / delta_rel_max);
//     }

//     for (auto& paramId : eta_infos_leaf) {
//         if (delta_rel[paramId] / delta_rel_max > config.nuisance_relevance_cutoff)
//             eta_specs_real_leaf[paramId] = pspp->get_param(paramId)->get_val();
//     }

//     LOG_INFO("Significant nuisances");
//     for (const auto& [pid, val] : eta_specs_real_leaf) {
//         LOG_INFO(pid, val);
//     }

//     // for (auto it = eta_specs_real_leaf.begin(); it != eta_specs_real_leaf.end(); ) {
//     //     std::ostringstream oss;
//     //     oss << it->first;
//     //     if (oss.str() == "B_Ks:1_5_2") {
//     //         std::cout << "[FIT] DEBUG removing nuisance " << oss.str() << "\n";
//     //         it = eta_specs_real_leaf.erase(it);
//     //     } else {
//     //         ++it;
//     //     }
//     // }
//     return eta_specs_real_leaf;
// }

std::map<ParamId, double> StatisticManager::get_p_specs(const std::vector<ParamId>& p_specs) {
    std::map<ParamId, double> out;
    for (auto elem : p_specs) {            
        out[elem] = pspp->get_param(elem)->get_val();
    }
    return out;
}
std::map<ParamId, std::map<ParamId, double>> StatisticManager::get_all_correlations() {
    std::map<ParamId, std::map<ParamId, double>> res;
    CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);
    res = ct.transform(this->cache.eta_specs_real);
    return res;
}

std::map<ExperimentObs, std::map<ExperimentObs, double>> StatisticManager::get_all_obs_correlations() {
    std::map<ExperimentObs, std::map<ExperimentObs, double>> res;
    CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);

    res = ct.transform(this->cache.exp_obs);
    return res;
}

// std::map<ExperimentObs, double> StatisticManager::get_obs_exp() {
//     std::map<ExperimentObs, double> out;

//     for (const auto& obsId : obs_int->get_obs_ids()) {
//         auto _ = pspp->get_obs_param(obsId);
//         for (auto _2 : _) {
//             out[_2.first] = _2.second->get_val(); 
//         }
//         // out[obsId] = pspp->get_obs_param(obsId)->get_val();
//     }
//     return out;
// }

std::map<ExperimentObs, double> StatisticManager::get_obs_exp() {
    std::map<ExperimentObs, double> out;

    for (const auto& obsId : obs_int->get_obs_ids()) {
        auto exp_params = pspp->get_obs_param(obsId);

        for (const auto& [exp_obs, param] : exp_params) {
            if (selected_experiments_.has_value()
                && !selected_experiments_->contains(exp_obs.experiment)) {
                continue;
            }

            out[exp_obs] = param->get_val();
        }
    }

    if (out.empty()) {
        throw std::runtime_error(
            "StatisticManager::get_obs_exp: experiment selection matched no experimental observable."
        );
    }

    return out;
}

void StatisticManager::prepare_likelihood_for_scan(const std::vector<ParamId>& p_specs) {
    update_cache(p_specs);

    if (cache.p_specs.empty()) {
        throw std::invalid_argument("prepare_likelihood_for_scan called with an empty fit parameter list.");
    }

    auto unzipped_fit_params = unzip(cache.p_specs);
    auto unzipped_nuisances  = unzip(cache.eta_specs_real);
    auto unzipped_exp_obs    = unzip(cache.exp_obs);

    const std::vector<ParamId> p_ids = unzipped_fit_params.ids;
    const std::vector<double> p0     = unzipped_fit_params.vals;

    const std::vector<ParamId> eta_ids = unzipped_nuisances.ids;
    const std::vector<double> eta0     = unzipped_nuisances.vals;

    const std::vector<ExperimentObs> obs_ids = unzipped_exp_obs.ids;
    const std::vector<double> exp_obs_vals   = unzipped_exp_obs.vals;

    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->nuisance_dist = build_nuisance_distribution();
    ctx->exp_obs_dist  = build_exp_data_distribution();
    ctx->exp_obs_values = exp_obs_vals;

    ctx->fp_defs.reserve(p_ids.size());
    for (std::size_t i = 0; i < p_ids.size(); ++i) {
        double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
        ctx->fp_defs.emplace_back(make_fit_param_def(p_ids[i], p0[i], sigma));
    }

    ctx->nuis_defs.reserve(eta_ids.size());
    for (std::size_t i = 0; i < eta_ids.size(); ++i) {
        double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
        ctx->nuis_defs.emplace_back(
            make_nuisance_parameter_definition(eta_ids[i], eta0[i], sigma)
        );
    }

    auto model_fn = [this, obs_ids, p_ids, eta_ids](
        const std::vector<double>& p_vec,
        const std::vector<double>& eta_vec) -> std::vector<double>
    {
        auto pred_map = this->obs_int->predict_optimized(
            zip(p_ids, p_vec),
            zip(eta_ids, eta_vec)
        );

        return ordered_prediction_vector(obs_ids, pred_map);
    };

    last_ctx_ = ctx;
    last_like_ = std::make_shared<BaseLikelihood>(model_fn, ctx, p_ids.size());
    last_fitter_ = std::make_shared<MLFitter>(ctx, model_fn);

    last_fit_param_ids_ = p_ids;
    last_nuisance_ids_ = eta_ids;
    last_fit_param_index_.clear();
    for (std::size_t i = 0; i < p_ids.size(); ++i) {
        last_fit_param_index_[p_ids[i]] = i;
    }

    last_scan_p_ = p0;
    last_scan_eta_ = eta0;
    has_manual_scan_point_ = false;

    std::cout << "[SCAN] Likelihood prepared without MLE.\n";
    std::cout << "[SCAN] n_fit_params = " << p_ids.size()
              << ", n_nuisances = " << eta_ids.size() << "\n";
}

void StatisticManager::set_manual_scan_point(const std::map<ParamId, double>& p_hat,
                                             const std::map<ParamId, double>& eta_hat) {
    if (!last_like_) {
        throw std::runtime_error(
            "Please call prepare_likelihood_for_scan(...) or compute_MLE(...) before set_manual_scan_point(...)."
        );
    }

    last_scan_p_.resize(last_fit_param_ids_.size());
    for (std::size_t i = 0; i < last_fit_param_ids_.size(); ++i) {
        auto it = p_hat.find(last_fit_param_ids_[i]);
        if (it == p_hat.end()) {
            throw std::invalid_argument("Missing manual fit-parameter value for one scanned parameter.");
        }
        last_scan_p_[i] = it->second;
    }

    last_scan_eta_.resize(last_nuisance_ids_.size());
    for (std::size_t i = 0; i < last_nuisance_ids_.size(); ++i) {
        auto it = eta_hat.find(last_nuisance_ids_[i]);
        if (it == eta_hat.end()) {
            throw std::invalid_argument("Missing manual nuisance value for one nuisance parameter.");
        }
        last_scan_eta_[i] = it->second;
    }

    has_manual_scan_point_ = true;
    std::cout << "[SCAN] Manual scan point loaded.\n";
}

LikelihoodScanGrid StatisticManager::scan_likelihood_around_current_point(
    ParamId p1,
    ParamId p2,
    double x_half_width,
    double y_half_width,
    std::size_t nx,
    std::size_t ny
) const {
    if (!last_like_) {
        throw std::runtime_error(
            "Please call prepare_likelihood_for_scan(...) or compute_MLE(...) before requesting a likelihood scan."
        );
    }

    if (!last_fit_param_index_.contains(p1) || !last_fit_param_index_.contains(p2)) {
        throw std::invalid_argument(
            "Likelihood scan requested for parameters that are not in the current prepared parameter set."
        );
    }

    if (p1 == p2) {
        throw std::invalid_argument("Likelihood scan requires two distinct parameters.");
    }

    const std::size_t ix = last_fit_param_index_.at(p1);
    const std::size_t iy = last_fit_param_index_.at(p2);

    std::vector<double> p_ref;
    std::vector<double> eta_ref;

    if (has_manual_scan_point_) {
        p_ref = last_scan_p_;
        eta_ref = last_scan_eta_;
    } else {
        if (last_fit_raw_.p_hat.empty() || last_fit_raw_.eta_hat.empty()) {
            throw std::runtime_error(
                "No reference point available. Use compute_MLE(...) or set_manual_scan_point(...)."
            );
        }
        p_ref = last_fit_raw_.p_hat;
        eta_ref = last_fit_raw_.eta_hat;
    }

    if (ix >= p_ref.size() || iy >= p_ref.size()) {
        throw std::runtime_error("Internal error: parameter index out of range.");
    }

    std::vector<double> theta0 = p_ref;
    theta0.insert(theta0.end(), eta_ref.begin(), eta_ref.end());

    const double nll0 = last_like_->nll(theta0);

    LikelihoodScanGrid out;
    out.x_param = p1;
    out.y_param = p2;
    out.x_center = p_ref[ix];
    out.y_center = p_ref[iy];
    out.nx = nx;
    out.ny = ny;
    out.points.reserve(nx * ny);

    double nll_min = std::numeric_limits<double>::infinity();

    const double x_min = out.x_center - x_half_width;
    const double x_max = out.x_center + x_half_width;
    const double y_min = out.y_center - y_half_width;
    const double y_max = out.y_center + y_half_width;

    for (std::size_t i = 0; i < nx; ++i) {
        const double x = x_min + (x_max - x_min) * static_cast<double>(i) / static_cast<double>(nx - 1);

        for (std::size_t j = 0; j < ny; ++j) {
            const double y = y_min + (y_max - y_min) * static_cast<double>(j) / static_cast<double>(ny - 1);

            std::vector<double> theta = theta0;
            theta[ix] = x;
            theta[iy] = y;

            const double nll = last_like_->nll(theta);

            LikelihoodScanPoint pt;
            pt.x = x;
            pt.y = y;
            pt.nll = nll;

            nll_min = std::min(nll_min, nll);
            out.points.push_back(pt);
        }
    }

    for (auto& pt : out.points) {
        pt.delta_nll = pt.nll - nll_min;
    }

    std::cout << "[SCAN] Built likelihood scan around current point\n";
    std::cout << "[SCAN] center x = " << out.x_center << "\n";
    std::cout << "[SCAN] center y = " << out.y_center << "\n";
    std::cout << "[SCAN] nll at reference point = " << nll0 << "\n";
    std::cout << "[SCAN] min nll on grid        = " << nll_min << "\n";

    return out;
}

void StatisticManager::save_likelihood_scan_csv(const std::string& path,
                                                const LikelihoodScanGrid& grid) const {
    std::ofstream out(path);
    out << "# x=" << grid.x_param << "\n";
    out << "# y=" << grid.y_param << "\n";
    out << "# x_center=" << std::setprecision(17) << grid.x_center << "\n";
    out << "# y_center=" << std::setprecision(17) << grid.y_center << "\n";
    out << "# nx=" << grid.nx << "\n";
    out << "# ny=" << grid.ny << "\n";
    out << "x,y,nll,delta_nll\n";

    out << std::setprecision(17);
    for (const auto& pt : grid.points) {
        out << pt.x << ","
            << pt.y << ","
            << pt.nll << ","
            << pt.delta_nll << "\n";
    }
}

void StatisticManager::select_experiment(const std::string& experiment) {
    select_experiments(std::set<std::string>{experiment});
}

void StatisticManager::select_experiments(const std::vector<std::string>& experiments) {
    select_experiments(std::set<std::string>(
        experiments.begin(),
        experiments.end()
    ));
}

void StatisticManager::select_experiments(const std::set<std::string>& experiments) {
    if (experiments.empty()) {
        throw std::invalid_argument(
            "StatisticManager::select_experiments: empty experiment set."
        );
    }

    selected_experiments_ = experiments;
    invalidate_fit_state();
}

void StatisticManager::select_experiments_all() {
    selected_experiments_.reset();
    invalidate_fit_state();
}

bool StatisticManager::has_experiment_selection() const noexcept {
    return selected_experiments_.has_value();
}

std::set<std::string> StatisticManager::selected_experiments() const {
    if (!selected_experiments_) {
        return {};
    }

    return *selected_experiments_;
}