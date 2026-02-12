#include "ScanRunner.h"

ScanVar to_scan_var(const YamlScanParam& p) {
    ScanVar v;
    v.block = p.block_name;
    v.pdg = p.pdg_code;
    v.name = p.block_name + ":" + std::to_string(p.pdg_code);
    v.min = p.min_val;
    v.max = p.max_val;
    v.step = p.step_val;
    return v;
}


DataSet ScanRunner::run(const OutputSpec& spec) {
    DataSet ds;

    ds.vars.reserve(scan_params_.size());
    for (auto& p : scan_params_) ds.vars.push_back(to_scan_var(p));

    if (extractor_) ds.outputs_schema = extractor_->schema(spec);

    std::vector<double> x(ds.vars.size(), 0.0);
    std::vector<double> cache(ds.vars.size(), 0.0);

    for (size_t i = 0; i < scan_params_.size(); ++i) {
        cache[i] = upp_.get_value(scan_params_[i].block_name, scan_params_[i].pdg_code).value_or(0.0);
    }

    recurse_dim(0, ds, x, spec);

    for (size_t i = 0; i < scan_params_.size(); ++i) {
        upp_.set_value(scan_params_[i].block_name, scan_params_[i].pdg_code, cache[i]);
    }

    return ds;
}


void ScanRunner::recurse_dim(size_t dim, DataSet& ds, std::vector<double>& x, const OutputSpec& spec) {
    if (dim >= scan_params_.size()) {
        DataPoint p;
        p.x = x;
        if (extractor_) extractor_->extract(p.y, spec);
        ds.points.push_back(std::move(p));
        return;
    }

    auto& sp = scan_params_[dim];
    for (double val = sp.min_val; val < sp.max_val; val += sp.step_val) {
        upp_.set_value(sp.block_name, sp.pdg_code, val);
        x[dim] = val;
        recurse_dim(dim + 1, ds, x, spec);
    }
}