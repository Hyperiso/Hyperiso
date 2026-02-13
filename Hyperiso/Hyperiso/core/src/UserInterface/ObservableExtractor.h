#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <optional>
#include <memory>

#include "IOutputExtractor.h"
#include "ObservableInterface.h"
#include "Include.h" // pour QCDOrder, OrderMapper, Value, etc.

// clé output
inline std::string obs_key(const std::string& obsName,
                           QCDOrder ord,
                           const std::optional<std::pair<double,double>>& bin)
{
    std::string k = "OBS:";
    k += obsName;
    k += ":";
    k += OrderMapper::str(ord);

    if (bin.has_value()) {
        k += ":BIN[";
        k += std::to_string(bin->first);
        k += ",";
        k += std::to_string(bin->second);
        k += "]";
    }
    return k;
}

class ObservableExtractor : public IOutputExtractor {
public:
    // ord = ordre demandé (parser.getValue("qcd_order"))
    ObservableExtractor(ObservableInterface& oi, QCDOrder ord)
    : oi_(oi), ord_(ord) {}

    std::vector<std::string> schema(const OutputSpec& spec) const override;
    void extract(std::unordered_map<std::string, Value>& outY,
                 const OutputSpec& spec) const override;

private:
    ObservableInterface& oi_;
    QCDOrder ord_;
};
