#include "MnContourExtractor.h"

Contour MnContourExtractor::extract(const ScalarField2D &field, const ContourRequest &cr) {
    auto minuit_backend = fit_app::make_minuit_backend();

    fit_app::ContourOptionsBackEnd co;
    co.npoints = 4 * cr.resolution;
    co.up = cr.level;

    auto f = fit_app::LambdaObjectiveFunction(
        [field] (std::vector<double>theta) { return field(theta[0], theta[1]); },
        cr.level
    );

    fit_app::FitOptions fopt;
    fopt.run_hesse = false;
    fopt.verbose = false;
    std::vector<fit_app::ParameterDefinition> p_defs = {cr.p_defs[0], cr.p_defs[1]};
    fit_app::BackendFitResult bfr = minuit_backend->minimize(f, p_defs, fopt);
    fit_app::BackendContourResult cres = minuit_backend->contour(f, bfr, 0, 1, co);

    Contour cont;
    cont.level = cr.level;
    cont.paths = {cres.points};
    cont.success = cres.success;

    return cont;
}