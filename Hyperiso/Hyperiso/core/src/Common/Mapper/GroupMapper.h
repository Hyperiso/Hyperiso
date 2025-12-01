#pragma once
#include "wgroup_ids.hpp"
#include "scaletype_ids.hpp"
#include "wilsonbasis_ids.hpp"
#include <sstream>

class GroupMapperExt : public GroupMapper {
public:
    static std::string str(WGroup g, ScaleType scale, WilsonBasis basis = WilsonBasis::B_STANDARD) {
        std::stringstream ss;
        ss << GroupMapper::str(g)
           << "_" << ScaleTypeMapper::str(scale)
           << (scale == ScaleType::HADRONIC ? "_" + WilsonBasisMapper::str(basis) : "");
        return ss.str();
    }
};

