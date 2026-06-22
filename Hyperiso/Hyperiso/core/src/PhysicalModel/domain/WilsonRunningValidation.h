#ifndef WILSON_RUNNING_VALIDATION_H
#define WILSON_RUNNING_VALIDATION_H

#include <algorithm>
#include <functional>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "Include.h"
#include "Wilson.h"
#include "SourcesView.h"

namespace WilsonRunningValidation {

inline std::string context(const std::string& group_name, WilsonBasis basis, QCDOrder order) {
    std::ostringstream os;
    os << "group=" << group_name
       << ", basis=" << WilsonBasisMapper::str(basis)
       << ", order=" << OrderMapper::str(order);
    return os.str();
}

inline void require_running_function(
    const std::map<QCDOrder,
        std::function<std::unordered_map<WCoefId, scalar_t>(
            const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>&,
            const BlockSrc&)>>& funcs,
    QCDOrder order,
    const std::string& group_name,
    WilsonBasis basis)
{
    const auto it = funcs.find(order);
    if (it == funcs.end() || !it->second) {
        throw std::invalid_argument(
            "Missing Wilson running function for " + context(group_name, basis, order));
    }
}

inline void require_matching_input_complete(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
    const std::vector<WCoefId>& members,
    QCDOrder order,
    const std::string& group_name,
    WilsonBasis basis)
{
    const auto it = matching.find(order);
    if (it == matching.end()) {
        throw std::invalid_argument(
            "Missing Wilson matching input order for " + context(group_name, basis, order));
    }

    for (const auto& coef_id : members) {
        if (!it->second.contains(coef_id)) {
            throw std::invalid_argument(
                "Missing Wilson matching input coefficient '" + coef_id.str() +
                "' for " + context(group_name, basis, order));
        }
    }
}

inline void require_running_result_known_members(
    const std::unordered_map<WCoefId, scalar_t>& result,
    const std::vector<WCoefId>& members,
    QCDOrder order,
    const std::string& group_name,
    WilsonBasis basis)
{
    for (const auto& [coef_id, _] : result) {
        if (std::find(members.begin(), members.end(), coef_id) == members.end()) {
            throw std::invalid_argument(
                "Wilson running function returned coefficient '" + coef_id.str() +
                "' outside the group member list for " + context(group_name, basis, order));
        }
    }
}

} // namespace WilsonRunningValidation

#endif // WILSON_RUNNING_VALIDATION_H
