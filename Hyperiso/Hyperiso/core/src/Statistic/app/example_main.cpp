#include <iostream>
#include <random>
#include <gsl/gsl_cdf.h>


#include "BToMuMuToy.h"
#include "LinearAlgebra.h"
#include "Statistics.h"
#include "GaussianApprox.h"
// #include "MonteCarloPredictor.h"
#include "Fit.h"


int main() {
// // Experimental inputs (vector + covariance)
// const Vec Oexp{3.52e-9, 1.3e-10};
// Matrix SigmaO = {
// { (1.2e-10)*(1.2e-10), -0.2*(1.2e-10)*(0.2e-10) },
// { -0.2*(1.2e-10)*(0.2e-10), (0.2e-10)*(0.2e-10) }
// };


// // Nuisances prior mean + covariance (diagonal here, but any SPD is ok)
// const Vec eta_mean{0.194, 0.234, 0.0635, 0.04111, 0.00858};
// Matrix SigmaEta = {
// {0.010*0.010, 0,0,0,0},
// {0, 0.010*0.010, 0,0,0},
// {0,0, 0.014*0.014, 0,0},
// {0,0,0, 0.00077*0.00077, 0},
// {0,0,0,0, 0.00019*0.00019}
// };


// // Build contexts
// auto SigmaO_chol = SPDMatrix::cholesky(SigmaO);
// auto SigmaEta_chol = SPDMatrix::cholesky(SigmaEta);
// LikelihoodContext ctx{Oexp, SigmaO_chol, eta_mean, SigmaEta_chol};


// // Model
// BToMuMuToy model;


// // MC prediction + skewness check
// MCPredictConfig cfg; cfg.draws = 10000; cfg.skew_abs_threshold = 0.2;
// MonteCarloPredictor mc(model, eta_mean, SigmaEta, cfg);
// std::mt19937 rng(12345);


// const Vec p_test{-4.5, 0.0};
// auto summaries = mc.summarize(p_test, rng);
// std::cout << "Skewness[BR_s] = " << summaries[0].skew << " (ok=" << summaries[0].approx_ok << ")\n";
// std::cout << "BR_s_untag = " << summaries[0].mu << " ± " << summaries[0].sigma << "\n";
// std::cout << "BR_d = " << summaries[1].mu << " ± " << summaries[1].sigma << "\n";


// // MLE and profiled likelihood
// MLEstimator est(ctx, [&model](const Vec& p, const Vec& eta){ return model.predict(p, eta); });


// const Vec p0{-4.5, 0.0};
// const Vec eta0 = eta_mean; // good initial guess
// FitResult fr = est.fit(p0, eta0);


// std::cout << "MLE: C10=" << fr.p_hat[0] << ", Cp10=" << fr.p_hat[1]
// << ", ell_hat=" << fr.ell_hat << "\n";


// // Profile T(C10) with Cp10 fixed at 0, confidence interval @95%
// auto T = [&](double c10){ return est.test_statistic(Vec{c10, 0.0}, fr, eta0); };


// const double thr95 = gsl_cdf_chisq_Pinv(0.95, 1);
// double c10_min=-7.0, c10_max=-1.0; int N=200;


// double left=std::numeric_limits<double>::quiet_NaN();
// double right=std::numeric_limits<double>::quiet_NaN();


// double prevp=c10_min, prevT=T(prevp);
// for (int i=1;i<=N;++i) {
// double p = c10_min + (c10_max-c10_min)*i/double(N);
// double t = T(p);
// if (std::isnan(left) && (prevT-thr95)*(t-thr95) <= 0.0) left = prevp;
// if ((prevT-thr95)*(t-thr95) <= 0.0) right = p;
// prevp=p; prevT=t;
// }
// std::cout << "95% CI for C10 (Cp10=0 profiled): [" << left << ", " << right << "]\n";


return 0;
}