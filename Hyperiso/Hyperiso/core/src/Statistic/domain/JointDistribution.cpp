#include "JointDistribution.h"

namespace {

constexpr double kUClip = 1e-13;

double finite_or_throw(double v, const std::string& label, std::size_t i) {
    if (!std::isfinite(v)) {
        std::ostringstream oss;
        oss << "JointDistribution::curvature non-finite " << label
            << " at index " << i;
        throw std::runtime_error(oss.str());
    }
    return v;
}

} // namespace

JointDistribution::JointDistribution(
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals,
    std::unique_ptr<ICopula> copula) :
    marginals_(std::move(marginals)), copula_(std::move(copula))
{}

std::vector<std::vector<double>> JointDistribution::sample(std::size_t n) const {
    std::vector<std::vector<double>> u = this->copula_->sample_u(n);
    std::vector<std::vector<double>> x (n, std::vector<double>(u[0].size(), 0.0));

    for (size_t i = 0; i < marginals_.size(); i++) {
        for (size_t j = 0; j < n; j++) {
            x[j][i] = marginals_.at(i)->ppf(u[j][i]);
        }
    }

    return x;
}

std::vector<double> JointDistribution::sample() const {
    std::vector<double> u = this->copula_->sample_u();
    std::vector<double> x (u.size(), 0.0);

    for (size_t i = 0; i < marginals_.size(); i++) {
        x[i] = marginals_.at(i)->ppf(u[i]);
    }

    return x;
}

double JointDistribution::logpdf(std::vector<double> x) const {
    if (x.size() != marginals_.size()) 
        throw std::invalid_argument("Wrong size of random vector.");

    std::vector<double> u = std::vector<double>(x.size(), 0.0);
    double log_marg {0.0};

    for (size_t i = 0; i < marginals_.size(); i++) {
        u[i] = std::clamp(marginals_[i]->cdf(x[i]), kUClip, 1.0 - kUClip);
        log_marg += marginals_[i]->logpdf(x[i]);
    }

    return log_marg + copula_->log_density(u);    
}

RealMatrix JointDistribution::curvature(std::vector<double> x) const {
    if (x.size() != marginals_.size()) {
        throw std::invalid_argument("Wrong size of random vector.");
    }

    const std::size_t d = x.size();

    std::vector<double> u(d, 0.0);
    RealMatrix W(d, d);

    Vector f_i(d, 0.0);
    Vector df_i(d, 0.0);
    Vector d2logf_i(d, 0.0);

    for (std::size_t i = 0; i < d; ++i) {
        const double raw_u = marginals_[i]->cdf(x[i]);
        u[i] = std::clamp(raw_u, kUClip, 1.0 - kUClip);

        // Important numerical detail:
        // If raw_u is far outside the clipped range, evaluating f, f', f''
        // at the original x can underflow to zero while the copula derivatives
        // are evaluated at the clipped u. That creates terms like 0 * inf or 0/0.
        //
        // Therefore the curvature uses the same effective point as the clipped u.
        // For normal marginals this is equivalent to saturating the local curvature
        // in extreme tails, and prevents non-finite W entries.
        const double x_eff = marginals_[i]->ppf(u[i]);
        PDFDiff fdf = marginals_[i]->f_df_ddf(x_eff);

        if (!std::isfinite(fdf.f) || !(fdf.f > 0.0)) {
            std::ostringstream oss;
            oss << "JointDistribution::curvature invalid marginal density"
                << " at index " << i
                << " x=" << x[i]
                << " x_eff=" << x_eff
                << " u=" << u[i]
                << " f=" << fdf.f;
            throw std::runtime_error(oss.str());
        }

        finite_or_throw(fdf.df, "marginal df", i);
        finite_or_throw(fdf.ddf, "marginal ddf", i);

        f_i[i] = fdf.f;
        df_i[i] = fdf.df;

        const double dlogf = fdf.df / fdf.f;
        d2logf_i[i] = fdf.ddf / fdf.f - dlogf * dlogf;

        finite_or_throw(f_i[i], "marginal f", i);
        finite_or_throw(df_i[i], "marginal df", i);
        finite_or_throw(d2logf_i[i], "marginal d2logf", i);
    }

    LogDensityDiff cdc = copula_->log_c_dc_ddc(u);

    for (std::size_t i = 0; i < d; ++i) {
        finite_or_throw(cdc.dlog_c.at(i, 0), "copula dlog_c", i);

        for (std::size_t j = 0; j < d; ++j) {
            const double ddlogc = cdc.ddlog_c.at(i, j);
            if (!std::isfinite(ddlogc)) {
                std::ostringstream oss;
                oss << "JointDistribution::curvature non-finite copula ddlog_c"
                    << " at (" << i << "," << j << ")";
                throw std::runtime_error(oss.str());
            }

            double wij = -ddlogc * f_i[i] * f_i[j];

            if (i == j) {
                wij += -cdc.dlog_c.at(i, 0) * df_i[i] - d2logf_i[i];
            }

            if (!std::isfinite(wij)) {
                std::ostringstream oss;
                oss << "JointDistribution::curvature produced non-finite W"
                    << " at (" << i << "," << j << ")"
                    << " u_i=" << u[i]
                    << " u_j=" << u[j]
                    << " f_i=" << f_i[i]
                    << " f_j=" << f_i[j]
                    << " df_i=" << df_i[i]
                    << " d2logf_i=" << d2logf_i[i]
                    << " dlogc_i=" << cdc.dlog_c.at(i, 0)
                    << " ddlogc_ij=" << ddlogc;
                throw std::runtime_error(oss.str());
            }

            W.at(i, j) = wij;
        }
    }

    // Force exact symmetry. Curvature is a Hessian, so asymmetries here are numerical.
    for (std::size_t i = 0; i < d; ++i) {
        for (std::size_t j = i + 1; j < d; ++j) {
            const double v = 0.5 * (W.at(i, j) + W.at(j, i));
            W.at(i, j) = v;
            W.at(j, i) = v;
        }
    }

    return W;  
}

std::size_t JointDistribution::dim() {
    return marginals_.size();
}

std::vector<double> JointDistribution::get_stds() {
    std::vector<double> stds;
    for (auto& m : this->marginals_)
        stds.emplace_back(m->std());

    return stds;
}
