#include <iostream>
#include <random>
#include <gsl/gsl_cdf.h>


#include "ObservableInterface.h"
#include "ObservableInterfaceAdapter2.h"
#include "adapters/RvgNuisanceSampler.h"
#include "MonteCarloPredictorGeneric.h"
#include "Fit.h"
#include "BToMuMuToy.h"
#include "LinearAlgebra.h"
#include "DistributionFactory.h"
#include "RandomVectorGenerator.h"
#include "HyperisoMaster.h"
#include "BlockProxy.h"
#include "StatParameterProxy.h"
#include "CovarianceTransformer.h"
#include "StatCorrelationProxy.h"

int main() {

    HyperisoMaster hyp = HyperisoMaster();
    Config config;
    config.model = Model::SM;

    hyp.init("lha/si_input.flha", config);

    std::map<ObservableId, QCDOrder> obs_ids = {
        {ObservableMapper::to_id(Observables::BR_BS_MUMU), QCDOrder::LO},
        {ObservableMapper::to_id(Observables::BR_BS_MUMU_UNTAG), QCDOrder::LO},
        {ObservableMapper::to_id(Observables::BR_BD_MUMU), QCDOrder::LO}
    };


    // ObservableInterface oi;
    // for (auto oid : obs_ids) oi.add_observable(oid, QCDOrder::LO, /*add_dependencies=*/true);

    // for (auto id : oi.get_all_ops_deps(Observables::BR_BS_MUMU)) {
    //     std::cout << id << std::endl;
    // }

    //Build experimental vector/covariance
    Vec Oexp{3.52e-9, 3.52e-9, 1.3e-10}; 
    Matrix SigmaO = {
        { (1.2e-10)*(1.2e-10), -0.2*(1.2e-10)*(1.2e-10), -0.2*(1.2e-10)*(0.2e-10) },
        { -0.2*(1.2e-10)*(1.2e-10), (1.2e-10)*(1.2e-10), 0.0 },
        { -0.2*(1.2e-10)*(0.2e-10), 0.0, (0.2e-10)*(0.2e-10) }
    };


    // Nuisances 
    Vec eta_mean{0.194, 0.234, 0.0635, 0.04111, 0.00858};
    Matrix SigmaEta = {
        {0.010*0.010, 0,0,0,0},
        {0, 0.010*0.010, 0,0,0},
        {0,0, 0.014*0.014, 0,0},
        {0,0,0, 0.00077*0.00077, 0},
        {0,0,0,0, 0.00019*0.00019}
    };

    std::vector<ParamId> p_specs = {
        {ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD", LhaID(3051313,4137,0,2)}, // C10
        {ParameterType::WILSON, "BPrimeCoefficients_B_SCALE_STANDARD", LhaID(3051313,4234,0,2)} // C10'
    };
    std::vector<ParamId> eta_specs = {
        {ParameterType::FLAVOR, "FCONST", /*code_fBd*/ 521},
        {ParameterType::FLAVOR, "FCONST", /*code_fBs*/ 531},
        {ParameterType::DECAY, "B_ll", /*code_y_s*/ 1},
        {ParameterType::SM, "VCKM", /*code_Vts*/ LhaID(2,1)},
        {ParameterType::SM, "VCKM", /*code_Vtd*/ LhaID(2,0)}
    };

    StatParameterProxy spp = StatParameterProxy();
    std::vector<ParamId> eta_specs_real;
    Vec eta_mean_real;
    for (auto elem : obs_ids) {
        for (auto _ : oi.get_all_ops_deps(elem))
        if (!(std::find(eta_specs_real.begin(), eta_specs_real.end(), _) != eta_specs_real.end())) {
            eta_specs_real.push_back(_);
        }
    }

    
    std::shared_ptr<IStatCorrelationProxy> pscp = std::make_shared<StatCorrelationProxy>();
    std::shared_ptr<IStatParameterProxy> pspp = std::make_shared<StatParameterProxy>();

    CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);
    std::vector<ParamId> eta_specs_real_with_corr = ct.check_if_corr(eta_specs_real);

    for (auto elem : eta_specs_real_with_corr) {
        // std::cout << elem << std::endl;
        eta_mean_real.push_back(spp(elem, DataType::VALUE));
    }

    Matrix SigmaEtaReal = ct.transform(eta_specs_real_with_corr);

    for (int i = 0; i < eta_specs_real_with_corr.size(); i++) {
        for (int j = 0; j < eta_specs_real_with_corr.size(); j++) {
            std::cout << "[" << eta_specs_real_with_corr[i] << ", " << eta_specs_real_with_corr[j] << "] = " << SigmaEtaReal[i][j] << " | ";
            // std::cout << "[" << i << ", " << j << "] = " << SigmaEtaReal[i][j] << " | ";
        }
        std::cout << std::endl;
    }

    ObservableInterfaceAdapterObs model(obs_ids, p_specs, eta_specs_real_with_corr);

    std::cout << "creating RandomVectorGenerator" << std::endl;

    unsigned int seed = std::random_device{}();
    auto dist = DistributionFactory::create("gaussian", seed);
    auto decomp = std::make_unique<CholeskyDecomposition>();

    RandomVectorGenerator rvg(std::move(dist), std::move(decomp));

    std::cout << "RandomVectorGenerator created" << std::endl;

    RvgNuisanceSampler sampler(rvg);

    std::cout << "Creating MonteCarloPredictor" << std::endl;

    // MC prediction with pluggable sampler
    MonteCarloPredictor2 mc(model, sampler, eta_mean_real, SigmaEtaReal, {10, 0.2});
    std::mt19937 rng(1234);

    std::cout << "MonteCarloPredictor created" << std::endl;

    Vec p_test{-4.5, 0.0};
    auto sums = mc.summarize(p_test, rng);

    std::cout << "summarize ented" << std::endl; 

    std::cout << "Skewness[0]=" << sums[0].skew << " ok=" << sums[0].approx_ok << std::endl;

    for (auto sum : sums) {
        std::cout << "value = " << sum.mu << " +- " << sum.sigma << std::endl;
    }

    std::cout << "Now doing likelihood : " << std::endl; 
    // Likelihood/MLE/intervals
    SPDMatrix SO = SPDMatrix::cholesky(SigmaO);
    SPDMatrix SE = SPDMatrix::cholesky(SigmaEtaReal);
    LikelihoodContext ctx{Oexp, SO, eta_mean_real, SE};
    MLEstimator est(ctx, [&model](const Vec& p, const Vec& eta){ return model.predict(p, eta); });

    std::cout << "Now doing MLE : " << std::endl; 
    Vec p0{-4.5, 0.0}; Vec eta0 = eta_mean_real;
    auto fr = est.fit(p0, eta0, 100);

    std::cout << "MLE fit done: " << std::endl; 

    std::cout << "MLE: C10=" << fr.p_hat[0] << ", Cp10=" << fr.p_hat[1]
    << ", ell_hat=" << fr.ell_hat << std::endl;

    std::cout << "finding thr95 " << std::endl; 
    const double thr95 = gsl_cdf_chisq_Pinv(0.95, 1);

    std::cout << "thr95 : " << thr95 << std::endl; 

    auto T = [&](double c10){ return est.test_statistic(Vec{c10, 0.0}, fr, eta0); };
    double a=-7,b=-1; int N=10; double left=std::nan(""), right=std::nan("");
    double prev=a, prevT=T(prev);
    for(int i=1;i<=N;++i){ double x=a+(b-a)*i/double(N); double t=T(x);
    if (std::isnan(left) && (prevT-thr95)*(t-thr95)<=0) left=prev; if ((prevT-thr95)*(t-thr95)<=0) right=x; prev=x; prevT=t; }
    std::cout << "95% CI C10|Cp10=0: ["<<left<<","<<right<<"]";


    return 0;
}