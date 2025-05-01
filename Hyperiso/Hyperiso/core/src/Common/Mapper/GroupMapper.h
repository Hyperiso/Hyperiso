#ifndef GROUP_MAPPER_H
#define GROUP_MAPPER_H


#include "EnumMapper.h"
#include "GeneralEnum.h"

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
};

#endif