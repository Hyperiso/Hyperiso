#ifndef GROUP_MAPPER_H
#define GROUP_MAPPER_H

#include "EnumMapper.h"
#include "GeneralEnum.h"
#include "ScaleTypeMapper.h"
#include "BWilsonBasisMapper.h"

class GroupMapper : public EnumMapperBase<WGroup, GroupMapper> {
public:
    static const std::map<WGroup, std::string>& mapping() {
        static const std::map<WGroup, std::string> m = {
            {WGroup::B, "BCoefficients"},
            {WGroup::BPrime, "BPrimeCoefficients"},
            {WGroup::BScalar, "BScalarCoefficients"},
            {WGroup::Blnu, "BlnuCoefficients"},
            {WGroup::BCLNU, "BclnuCoefficients"},
        };
        return m;
    }

    static const std::map<std::string, WGroup>& inverse_mapping() {
        static const std::map<std::string, WGroup> inv = invert_map(mapping());
        return inv;
    }

    static std::string str(WGroup group) {
        return EnumMapperBase::str(group);
    }

    static std::string str(WGroup group, ScaleType scale, bool full=false, BWilsonBasis basis=BWilsonBasis::STANDARD) {
        std::stringstream ss;
        ss << mapping().at(group)
           << "_" << ScaleTypeMapper::str(scale)
           << (full ? "_FULL" : "")
           << (scale == ScaleType::HADRONIC ? "_" + BWilsonBasisMapper::str(basis) : "");

        return ss.str();
    }
};

#endif