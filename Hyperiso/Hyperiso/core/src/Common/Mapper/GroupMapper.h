#ifndef GROUP_MAPPER_H
#define GROUP_MAPPER_H

#include <sstream>

#include "wgroup_ids.hpp"
#include "scaletype_ids.hpp"
#include "wilsonbasis_ids.hpp"

/**
 * @file GroupMapperExt.h
 * @brief Extended helper for composing Wilson group labels.
 *
 * GroupMapperExt builds composite names that combine:
 *   - a Wilson group name (e.g. "BCoefficients"),
 *   - a scale label (e.g. "EW_SCALE", "B_SCALE"),
 *   - optionally a Wilson basis label at hadronic scale
 *     (e.g. "STANDARD", "TRADITIONAL").
 *
 * The resulting string typically has the form:
 *   "<group>_<scale>[_<basis>]"
 * with the basis part present only for hadronic scales.
 */

/**
 * @class GroupMapperExt
 * @brief Extension of GroupMapper for building composite group strings.
 */
class GroupMapperExt : public GroupMapper {
public:
    /**
     * @brief Returns a composite label "GROUP_SCALE[_BASIS]" for a Wilson group.
     *
     * Examples:
     *   - str(WGroup::B, ScaleType::MATCHING)        -> "BCoefficients_EW_SCALE"
     *   - str(WGroup::B, ScaleType::HADRONIC)        -> "BCoefficients_B_SCALE_STANDARD"
     *   - str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_TRADITIONAL)
     *                                            -> "BCoefficients_B_SCALE_TRADITIONAL"
     *
     * @param g      Wilson operator group.
     * @param scale  Scale type (matching / hadronic).
     * @param basis  Wilson basis (used only if scale == HADRONIC).
     * @return Composite string label.
     */
    static std::string str(WGroup g, ScaleType scale, WilsonBasis basis = WilsonBasis::B_STANDARD) {
        std::stringstream ss;
        ss << GroupMapper::str(g)
           << "_" << ScaleTypeMapper::str(scale)
           << (scale == ScaleType::HADRONIC ? "_" + WilsonBasisMapper::str(basis) : "");
        return ss.str();
    }
};

#endif