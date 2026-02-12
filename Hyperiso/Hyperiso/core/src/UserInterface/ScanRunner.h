#ifndef SCAN_RUNNER_H
#define SCAN_RUNNER_H

#include "YamlInputReader.h"
#include "UserParameterProxy.h"
#include "IOutputExtractor.h"
#include "ObjectsOutputs.h"

inline ScanVar to_scan_var(const YamlScanParam& p);

class ScanRunner {
public:
    ScanRunner(UserParameterProxy& upp,
               std::vector<YamlScanParam> scan_params,
               std::shared_ptr<IOutputExtractor> extractor)
    : upp_(upp), scan_params_(std::move(scan_params)), extractor_(std::move(extractor))
    {}

    DataSet run(const OutputSpec& spec);

private:
    UserParameterProxy& upp_;
    std::vector<YamlScanParam> scan_params_;
    std::shared_ptr<IOutputExtractor> extractor_;

    void recurse_dim(size_t dim, DataSet& ds, std::vector<double>& x, const OutputSpec& spec);
};

#endif