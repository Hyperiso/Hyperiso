
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "BaseLikelihood.h"
#include "GradientHelper.h"
#include "GaussianCopula.h"
#include "GaussianMarginal.h"

using Vec = std::vector<double>;
using Mat = std::vector<Vec>;

static Vec matvec_std(const Mat& A, const Vec& x) {
    if (A.empty()) return {};
    if (A[0].size() != x.size()) {
        throw std::runtime_error("matvec_std: dimension mismatch");
    }

    Vec out(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i) {
        for (std::size_t j = 0; j < x.size(); ++j) {
            out[i] += A[i][j] * x[j];
        }
    }
    return out;
}

static Vec add_std(const Vec& a, const Vec& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("add_std: dimension mismatch");
    }

    Vec out(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        out[i] = a[i] + b[i];
    }
    return out;
}

static Vec sub_std(const Vec& a, const Vec& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("sub_std: dimension mismatch");
    }

    Vec out(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        out[i] = a[i] - b[i];
    }
    return out;
}

static double dot_std(const Vec& a, const Vec& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("dot_std: dimension mismatch");
    }

    double out = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        out += a[i] * b[i];
    }
    return out;
}

static Mat transpose_std(const Mat& A) {
    if (A.empty()) return {};

    Mat T(A[0].size(), Vec(A.size(), 0.0));

    for (std::size_t i = 0; i < A.size(); ++i) {
        for (std::size_t j = 0; j < A[0].size(); ++j) {
            T[j][i] = A[i][j];
        }
    }

    return T;
}

static Mat matmul_std(const Mat& A, const Mat& B) {
    if (A.empty() || B.empty()) return {};
    if (A[0].size() != B.size()) {
        throw std::runtime_error("matmul_std: dimension mismatch");
    }

    Mat C(A.size(), Vec(B[0].size(), 0.0));

    for (std::size_t i = 0; i < A.size(); ++i) {
        for (std::size_t k = 0; k < B.size(); ++k) {
            for (std::size_t j = 0; j < B[0].size(); ++j) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}

static Mat matadd_std(const Mat& A, const Mat& B) {
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        throw std::runtime_error("matadd_std: dimension mismatch");
    }

    Mat C = A;

    for (std::size_t i = 0; i < A.size(); ++i) {
        for (std::size_t j = 0; j < A[0].size(); ++j) {
            C[i][j] += B[i][j];
        }
    }

    return C;
}

static Vec solve_linear_std(Mat A, Vec b) {
    const std::size_t n = b.size();

    if (A.size() != n || A[0].size() != n) {
        throw std::runtime_error("solve_linear_std: matrix must be square");
    }

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        double best = std::abs(A[col][col]);

        for (std::size_t row = col + 1; row < n; ++row) {
            const double candidate = std::abs(A[row][col]);
            if (candidate > best) {
                best = candidate;
                pivot = row;
            }
        }

        if (best < 1e-14) {
            throw std::runtime_error("solve_linear_std: singular matrix");
        }

        if (pivot != col) {
            std::swap(A[pivot], A[col]);
            std::swap(b[pivot], b[col]);
        }

        const double diag = A[col][col];

        for (std::size_t j = col; j < n; ++j) {
            A[col][j] /= diag;
        }
        b[col] /= diag;

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) continue;

            const double factor = A[row][col];

            for (std::size_t j = col; j < n; ++j) {
                A[row][j] -= factor * A[col][j];
            }
            b[row] -= factor * b[col];
        }
    }

    return b;
}

static RealMatrix identity_real_matrix(std::size_t n) {
    RealMatrix R(n, n);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            R.at(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }

    return R;
}

static std::unique_ptr<JointDistribution> make_independent_gaussian_joint(
    const Vec& means,
    const Vec& sigmas
) {
    if (means.size() != sigmas.size()) {
        throw std::runtime_error("make_independent_gaussian_joint: dimension mismatch");
    }

    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
    marginals.reserve(means.size());

    for (std::size_t i = 0; i < means.size(); ++i) {
        marginals.emplace_back(
            std::make_unique<GaussianMarginal>(means[i], sigmas[i])
        );
    }

    RealMatrix R = identity_real_matrix(means.size());

    std::unique_ptr<ICopula> copula =
        std::make_unique<GaussianCopula>(std::random_device{}(), R);

    return std::make_unique<JointDistribution>(
        std::move(marginals),
        std::move(copula)
    );
}

static fit_app::ParameterDefinition make_param_def(
    const std::string& name,
    double value,
    double step_hint
) {
    fit_app::ParameterDefinition def;
    def.name = name;
    def.value = value;
    def.step_hint = step_hint;
    def.fixed = false;
    return def;
}

struct KnownLinearGaussianModel {
    Mat A = {
        {1.0,  2.0},
        {0.5, -1.0},
        {1.5,  0.0},
        {0.0,  1.0}
    };

    Mat B = {
        {1.0,  0.0,  0.5},
        {0.0,  1.0, -0.2},
        {0.3, -0.4,  1.0},
        {1.0,  0.2,  0.0}
    };

    Vec y = {1.0, -0.5, 0.7, 0.2};

    Vec p0 = {0.0, 0.0};
    Vec eta0 = {0.1, -0.2, 0.05};

    Vec obs_sigmas = {0.5, 1.2, 0.8, 1.5};
    Vec eta_sigmas = {1.0, 0.7, 1.3};

    Vec predict(const Vec& p, const Vec& eta) const {
        return add_std(matvec_std(A, p), matvec_std(B, eta));
    }

    Mat W_obs_std() const {
        Mat W(obs_sigmas.size(), Vec(obs_sigmas.size(), 0.0));
        for (std::size_t i = 0; i < obs_sigmas.size(); ++i) {
            W[i][i] = 1.0 / (obs_sigmas[i] * obs_sigmas[i]);
        }
        return W;
    }

    Mat W_eta_std() const {
        Mat W(eta_sigmas.size(), Vec(eta_sigmas.size(), 0.0));
        for (std::size_t i = 0; i < eta_sigmas.size(); ++i) {
            W[i][i] = 1.0 / (eta_sigmas[i] * eta_sigmas[i]);
        }
        return W;
    }

    Mat exact_H_eta() const {
        return matadd_std(
            matmul_std(matmul_std(transpose_std(B), W_obs_std()), B),
            W_eta_std()
        );
    }

    Vec exact_g_eta(const Vec& p) const {
        const Vec r0 = sub_std(predict(p, eta0), y);
        return matvec_std(transpose_std(B), matvec_std(W_obs_std(), r0));
    }

    Vec exact_eta_hat(const Vec& p) const {
        const Mat H = exact_H_eta();
        const Vec g = exact_g_eta(p);
        const Vec Hinv_g = solve_linear_std(H, g);

        Vec out = eta0;
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] -= Hinv_g[i];
        }
        return out;
    }

    double exact_profiled_nll_using_base_constant(
        const IProfileableLikelihood& like,
        const Vec& p
    ) const {
        const double nll0_from_base = like.nll_from_split(p, eta0);
        const Mat H = exact_H_eta();
        const Vec g = exact_g_eta(p);
        const Vec Hinv_g = solve_linear_std(H, g);
        return nll0_from_base - 0.5 * dot_std(g, Hinv_g);
    }
};

static std::shared_ptr<BaseLikelihood> build_test_likelihood(
    const KnownLinearGaussianModel& kgm
) {
    auto ctx = std::make_shared<LikelihoodContext>();

    ctx->exp_obs_values = kgm.y;

    // exp_obs_dist est la distribution des résidus r = f(p, eta) - y,
    // donc elle est centrée en 0.
    ctx->exp_obs_dist = make_independent_gaussian_joint(
        Vec(kgm.y.size(), 0.0),
        kgm.obs_sigmas
    );

    // nuisance_dist est centrée en eta0.
    ctx->nuisance_dist = make_independent_gaussian_joint(
        kgm.eta0,
        kgm.eta_sigmas
    );

    ctx->fp_defs.push_back(make_param_def("p0", kgm.p0[0], 1e-2));
    ctx->fp_defs.push_back(make_param_def("p1", kgm.p0[1], 1e-2));

    ctx->nuis_defs.push_back(make_param_def("eta0", kgm.eta0[0], kgm.eta_sigmas[0]));
    ctx->nuis_defs.push_back(make_param_def("eta1", kgm.eta0[1], kgm.eta_sigmas[1]));
    ctx->nuis_defs.push_back(make_param_def("eta2", kgm.eta0[2], kgm.eta_sigmas[2]));

    ModelFn model_fn = [kgm](const Vec& p, const Vec& eta) {
        return kgm.predict(p, eta);
    };

    return std::make_shared<BaseLikelihood>(model_fn, ctx, 2);
}

static double max_abs_diff(const Vec& a, const Vec& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("max_abs_diff: dimension mismatch");
    }

    double out = 0.0;

    for (std::size_t i = 0; i < a.size(); ++i) {
        out = std::max(out, std::abs(a[i] - b[i]));
    }

    return out;
}

int main(int argc, char** argv) {
    std::string csv_name = "classes_laplace_contour.csv";
    if (argc >= 2) {
        csv_name = argv[1];
    }

    KnownLinearGaussianModel kgm;
    std::shared_ptr<BaseLikelihood> like = build_test_likelihood(kgm);

    const int n0 = 161;
    const int n1 = 161;

    const double p0_min = -1.5;
    const double p0_max =  1.5;
    const double p1_min = -1.5;
    const double p1_max =  1.5;

    std::ofstream csv(csv_name);
    if (!csv) {
        throw std::runtime_error("Could not open output CSV");
    }

    csv << std::setprecision(17);

    csv << "p0,p1,"
        << "nll0_from_base,"
        << "profiled_laplace_classes,profiled_exact,diff_profiled,"
        << "eta0_laplace,eta1_laplace,eta2_laplace,"
        << "eta0_exact,eta1_exact,eta2_exact,"
        << "eta_max_abs_diff,"
        << "direct_nll_at_eta_laplace,direct_nll_at_eta_exact\n";

    double min_laplace = std::numeric_limits<double>::infinity();
    double min_exact = std::numeric_limits<double>::infinity();
    double max_profile_diff = 0.0;
    double max_eta_diff = 0.0;

    for (int i = 0; i < n0; ++i) {
        const double p0 =
            p0_min + (p0_max - p0_min) * static_cast<double>(i) / static_cast<double>(n0 - 1);

        for (int j = 0; j < n1; ++j) {
            const double p1 =
                p1_min + (p1_max - p1_min) * static_cast<double>(j) / static_cast<double>(n1 - 1);

            const Vec p = {p0, p1};

            // Vraie brique de ton pipeline :
            // BaseLikelihood + JointDistribution::curvature + GradientHelper::laplace_profile_eta.
            const LaplaceProfileComputation lap = laplace_profile_eta(*like, p);

            // Référence analytique connue pour ce modèle linéaire-gaussien.
            const double exact_profiled =
                kgm.exact_profiled_nll_using_base_constant(*like, p);

            const Vec exact_eta = kgm.exact_eta_hat(p);

            const double nll0_from_base = like->nll_from_split(p, kgm.eta0);
            const double diff_profiled = lap.nll_hat - exact_profiled;
            const double eta_diff = max_abs_diff(lap.eta_hat, exact_eta);

            const double direct_lap = like->nll_from_split(p, lap.eta_hat);
            const double direct_exact = like->nll_from_split(p, exact_eta);

            min_laplace = std::min(min_laplace, lap.nll_hat);
            min_exact = std::min(min_exact, exact_profiled);
            max_profile_diff = std::max(max_profile_diff, std::abs(diff_profiled));
            max_eta_diff = std::max(max_eta_diff, eta_diff);

            csv << p0 << "," << p1 << ","
                << nll0_from_base << ","
                << lap.nll_hat << "," << exact_profiled << "," << diff_profiled << ","
                << lap.eta_hat[0] << "," << lap.eta_hat[1] << "," << lap.eta_hat[2] << ","
                << exact_eta[0] << "," << exact_eta[1] << "," << exact_eta[2] << ","
                << eta_diff << ","
                << direct_lap << "," << direct_exact << "\n";
        }
    }

    std::cout << "Wrote: " << csv_name << "\n";
    std::cout << "min_laplace=" << std::setprecision(17) << min_laplace << "\n";
    std::cout << "min_exact=" << std::setprecision(17) << min_exact << "\n";
    std::cout << "max_abs_diff_profiled=" << max_profile_diff << "\n";
    std::cout << "max_abs_diff_eta=" << max_eta_diff << "\n";

    if (max_profile_diff > 1e-6 || max_eta_diff > 1e-6) {
        std::cerr << "[FAIL] Laplace profiling through project classes does not match exact reference.\n";
        return 1;
    }

    std::cout << "[OK] BaseLikelihood + JointDistribution + GradientHelper match the exact reference.\n";
    return 0;
}
