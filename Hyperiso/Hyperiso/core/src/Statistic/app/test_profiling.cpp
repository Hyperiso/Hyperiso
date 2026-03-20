#include "ILikelihood.h"
#include "BaseLikelihood.h"
#include "IProfilingStrategy.h"
#include "ProfiledLikelihood2D.h"
#include "Profiler.h"
#include "WithGaussianConstraints.h"
#include "GaussianMarginal.h"
#include "GaussianCopula.h"
#include "Math.h"
#include "ContourEngine.h"

int main() {
    unsigned long seed = 1234567890;

    auto model_fn = [] (const std::vector<double>& p, const std::vector<double> eta) {
        double a = p[0];
        double b = p[1];
        double c = p[2];
        double d = eta[0]; 
        double e = eta[1];

        auto f = [a, b, c, d, e] (double x, double y) {
            return a * x * x + b * y * y + c * x * y + d * x + e * y;
        };

        return std::vector<double> ({f(1, 0), f(0, 1), f(1, 1), f(-1, 0), f(0, -1), f(-1, -1), f(1, -1), f(-1, 1)});
    };

    std::vector<std::unique_ptr<IMarginalDistribution>> nuis_marginals;
    nuis_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>(3, 1, seed)));
    nuis_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>(-1, 0.5, seed)));

    RealMatrix R_nuis ({
        {1, 0.1},
        {0.1, 1}
    });

    std::unique_ptr<ICopula> nuis_copula = std::make_unique<GaussianCopula>(seed, R_nuis);
    std::unique_ptr<JointDistribution> nuis_dist = std::make_unique<JointDistribution>(std::move(nuis_marginals), std::move(nuis_copula));

    std::vector<std::unique_ptr<IMarginalDistribution>> obs_marginals;
    obs_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>( 0, 0.1, seed)));
    obs_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>( 0, 0.3, seed)));
    obs_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>( 0, 0.2, seed)));
    obs_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>( 0, 0.1, seed)));
    obs_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>( 0, 0.3, seed)));
    obs_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>( 0, 0.2, seed)));
    obs_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>( 0, 0.5, seed)));
    obs_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>( 0, 0.3, seed)));

    RealMatrix R_obs ({
        {1, 0, 0.1, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0.5, 0},
        {0.1, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0.5, 0, 0, 0, 0, 1, 0.1},
        {0, 0, 0, 0, 0, 0, 0.1, 1}
    });

    std::unique_ptr<ICopula> obs_copula = std::make_unique<GaussianCopula>(seed, R_obs);
    std::unique_ptr<JointDistribution> obs_dist = std::make_unique<JointDistribution>(std::move(obs_marginals), std::move(obs_copula));

    std::shared_ptr<LikelihoodContext> ctx = std::make_shared<LikelihoodContext>();
    ctx->exp_obs_dist = std::move(obs_dist);
    ctx->nuisance_dist = std::move(nuis_dist);
    ctx->exp_obs_values = {5, 0, 3, -1, 2, -1, 9, 1};
    ctx->nuis_defs = {
        fit_app::ParameterDefinition {"d",  3.0, 1.0},
        fit_app::ParameterDefinition {"e", -1.0, 0.5}
    };
    ctx->fp_defs = {
        fit_app::ParameterDefinition {"a"},
        fit_app::ParameterDefinition {"b"},
        fit_app::ParameterDefinition {"c"}
    };

    std::shared_ptr<ILikelihood> base = std::make_shared<BaseLikelihood>(model_fn, ctx, 3);

    std::cout << "Base likelihood:" << std::endl;
    std::cout << "dim = " << base->dim() << std::endl;
    std::cout << "nll(2, 1, -2, 3, -1) = " << base->nll({2, 1, -2, 3, -1}) << std::endl;

    auto of = [base] (std::vector<double> theta) {
        return base->nll(theta);
    };

    MinimizationContext min_ctx;
    MinimizationResult res = minimize_combined(of, {0, 0, 0, 3.0, -1.0}, {1, 1, 1, 1, 1}, min_ctx);

    auto of_profiled = [base, min_ctx] (std::vector<double> p) {
        auto wrapped = [base, p] (std::vector<double> eta) {
            std::vector<double> theta = p;
            theta.insert(theta.end(), eta.begin(), eta.end());
            return base->nll(theta);
        };

        return minimize_combined(wrapped, {3, -1}, {1, 1}, min_ctx).min;
    };

    RealMatrix cov = inverse_hessian(of_profiled, {res.argmin[0], res.argmin[1], res.argmin[2]}, {1, 1, 1});

    RealMatrix corr = cov;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++)
            corr.at(i, j) /= std::sqrt(cov.at(i, i) * cov.at(j, j));
    }

    FitResult fr;
    fr.ell_hat = res.min;
    fr.eta_hat = {res.argmin[3], res.argmin[4]};
    fr.p_hat = {res.argmin[0], res.argmin[1], res.argmin[2]};
    fr.p_hat_std = {std::sqrt(cov.at(0, 0)), std::sqrt(cov.at(1, 1)), std::sqrt(cov.at(2, 2))};
    fr.p_hat_correlations = corr;

    std::cout << "Toy master fit:" << std::endl;
    std::cout << "min nll =" << fr.ell_hat << std::endl;
    std::cout << "p_hat = (" << fr.p_hat[0] << ", " << fr.p_hat[1] << ", " << fr.p_hat[2] << ")" << std::endl;
    std::cout << "eta_hat = (" << fr.eta_hat[0] << ", " << fr.eta_hat[1] << ")" << std::endl;
    std::cout << "p_hat_std = (" << fr.p_hat_std[0] << ", " << fr.p_hat_std[1] << ", " << fr.p_hat_std[2] << ")" << std::endl;
    std::cout << "p_hat_corr:" << std::endl;
    std::cout << fr.p_hat_correlations << std::endl;


    // std::shared_ptr<Profiler> profiler = std::make_shared<Profiler>(fit_app::make_minuit_backend());
    // std::shared_ptr<IProfilingStrategy> slice = std::make_shared<SliceProfilingStrategy>(0, 1, fr);
    // std::shared_ptr<IProfilingStrategy> project = std::make_shared<ProjectionProfilingStrategy>(0, 1, fr);

    // ProfiledLikelihood2D pl_1 (base, profiler, slice);
    // ProfiledLikelihood2D pl_2 (base, profiler, project);

    // std::vector<std::unique_ptr<IMarginalDistribution>> fitted_marginals;
    // fitted_marginals.emplace_back(std::move(std::make_unique<GaussianMarginal>(fr.p_hat[2], fr.p_hat_std[2], seed)));

    // RealMatrix R_fitted = eye(1);

    // std::unique_ptr<ICopula> fitted_copula = std::make_unique<GaussianCopula>(seed, R_fitted);
    // std::unique_ptr<JointDistribution> fitted_dist = std::make_unique<JointDistribution>(std::move(fitted_marginals), std::move(fitted_copula));
    // std::shared_ptr<ILikelihood> constrained_base = std::make_shared<WithGaussianConstraints>(base, std::move(fitted_dist), std::vector<std::size_t> ({2}));

    // ProfiledLikelihood2D pl_3 (constrained_base, profiler, project);

    // double ell_hat_1 = pl_1.profiled_nll(2.1, 0.9);
    // double ell_hat_2 = pl_2.profiled_nll(2.1, 0.9);
    // double ell_hat_3 = pl_3.profiled_nll(2.1, 0.9);

    // std::cout << "Profiled likelihoods:" << std::endl;
    // std::cout << "Method 1 (Slice):" << ell_hat_1 << std::endl;
    // std::cout << "Method 2 (Free projection):" << ell_hat_2 << std::endl;
    // std::cout << "Method 3 (Prior-constrained projection):" << ell_hat_3 << std::endl;

    ContourConfig cc;
    cc.fr = fr;
    cc.x_id = 0;
    cc.y_id = 1;
    cc.primary_contour_method = ContourAlgorithm::AMS;
    // cc.fallback_contour_method = ContourAlgorithm::AMS;
    cc.profiling_method = ProfilingMethod::SLICE;

    ContourEngine ce(base, cc);
    Contour cl = ce.compute_contour(1.0, {1.8, 2.2, 0.8, 1.2}, 200);

    std::ofstream os;
    os.open("contour_proj_constrained.csv");
    os << "x,y\n";
    for(const Path& path: cl.paths) {
        for (const Point& point: path) {
            os << point.first << "," << point.second << "\n";
        }
    }

    os.close();

    return 0;
}