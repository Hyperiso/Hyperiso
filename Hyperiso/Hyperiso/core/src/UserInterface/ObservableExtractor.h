#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <optional>
#include <memory>

#include "IOutputExtractor.h"
#include "ObservableInterface.h"
#include "Include.h"

std::string obs_key(const ObservableId& obs, QCDOrder ord);
std::string obs_key_bin(const ObservableId& obs, QCDOrder ord, double bmin, double bmax);

class ObservableExtractor : public IOutputExtractor {
public:
    ObservableExtractor(ObservableInterface& oi, QCDOrder ord)
        : oi_(oi), ord_(ord) {}

    std::vector<std::string> schema(const OutputSpec& spec) const override;

    void extract(std::unordered_map<std::string, Value>& outY,
                 const OutputSpec& spec) const override;

private:
    ObservableInterface& oi_;
    QCDOrder ord_;
};