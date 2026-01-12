#include "JointDistribution.h"
#include "MarginalFactory.h"
#include "CopulaFactory.h"
#include "Matrix.h"
#include "FlatMarginal.h"
#include "GaussianMarginal.h"
#include "SplitGaussianMarginal.h"
#include "GaussianCopula.h"

int main() {

    unsigned int seed = 123456789;

    RealMatrix R ({
        {1.0, 0.8}, 
        {0.8, 1.0}}
    );

    StudentTCopulaConfig c_cfg;
    c_cfg.R = R;
    c_cfg.nu = 4;

    GaussianMarginalCfg m_cfg_1 {1.0, 0.2};
    GaussianMarginalCfg m_cfg_2 {3.0, 0.5}; 

    auto m_1 = DistributionFactory::create(MarginalType::GAUSSIAN, m_cfg_1, seed);
    auto m_2 = DistributionFactory::create(MarginalType::GAUSSIAN, m_cfg_2, seed);
    auto cop = CopulaFactory::create(CopulaType::STUDENT_T, c_cfg, seed);

    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
    marginals.push_back(std::move(m_1));
    marginals.push_back(std::move(m_2));

    JointDistribution rvg(std::move(marginals), std::move(cop));

    LOG_INFO("JointDistribution initialized");

    std::vector<Vector> smpl = rvg.sample(10000);
    
    std::ofstream os;
    os.open("sample.csv");

    for (const Vector& z : smpl) {
        os << z[0] << "," << z[1] << "\n";
    }

    os.close();

    return 0;
}