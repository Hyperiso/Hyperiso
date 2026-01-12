#include "RandomVectorGenerator.h"
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
        {1.0, 0.5}, 
        {0.5, 1.0}}
    );

    GaussianCopulaConfig c_cfg {R};

    GaussianMarginalCfg m_cfg_1 {1.0, 2.0};
    FlatMarginalCfg m_cfg_2 {3.0, 4.0}; 

    auto m_1 = DistributionFactory::create(MarginalType::GAUSSIAN, m_cfg_1, seed);
    auto m_2 = DistributionFactory::create(MarginalType::FLAT, m_cfg_2, seed);
    auto cop = CopulaFactory::create(CopulaType::GAUSSIAN, c_cfg, seed);

    JointDistribution rvg {{std::move(m_1), std::move(m_2)}, std::move(cop)};

    std::vector<Vector> smpl = rvg.sample(1000);
    
    std::ofstream os;
    os.open("sample.csv");

    for (const Vector& z : smpl) {
        os << z[0] << "," << z[1] << "\n";
    }

    os.close();

    return 0;
}