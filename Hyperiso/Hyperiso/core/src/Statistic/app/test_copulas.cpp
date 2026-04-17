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
        {1.0, 0.7}, 
        {0.7, 1.0}}
    );

    GaussianCopulaConfig c_cfg;
    c_cfg.R = R;
    // c_cfg.nu = 4;

    FlatMarginalCfg m_cfg_1 {1.0, 2.0};
    FlatMarginalCfg m_cfg_2 {2.0, 4.0}; 

    auto m_1 = MarginalFactory::create(MarginalType::FLAT, m_cfg_1, seed);
    auto m_2 = MarginalFactory::create(MarginalType::FLAT, m_cfg_2, seed);
    auto cop = CopulaFactory::create(CopulaType::GAUSSIAN, c_cfg, seed);

    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
    marginals.push_back(std::move(m_1));
    marginals.push_back(std::move(m_2));

    JointDistribution rvg(std::move(marginals), std::move(cop));

    LOG_INFO("JointDistribution initialized");

    // Sampling joint distribution
    std::vector<Vector> smpl = rvg.sample(10000);
    
    std::ofstream os;
    os.open("sample.csv");

    for (const Vector& z : smpl) {
        os << z[0] << "," << z[1] << "\n";
    }

    os.close();

    // Scanning pdf
    double x_1 {0.0}, x_2 {1.0};
    double dx = 1e-2;
    double x1_max {2.0};
    double x2_max {5.0};
    std::size_t n1 = x1_max / dx;
    std::size_t n2 = x2_max / dx;

    printf("(%li,%li)\n", n1, n2);

    os.open("logpdf.csv");

    for (size_t i = 0; i < n1; i++) {
        x_2 = 0.0;
        for (size_t j = 0; j < n2; j++) {
            os << x_1 << "," << x_2 << "," << rvg.logpdf({x_1, x_2}) << "\n"; 
            x_2 += dx;            
        }
        x_1 += dx;
    }

    os.close();
    
    return 0;
}