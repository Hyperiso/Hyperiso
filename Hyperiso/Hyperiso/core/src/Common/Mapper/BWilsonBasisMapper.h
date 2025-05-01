#ifndef BWILSON_BASIS_MAPPER_H
#define BWILSON_BASIS_MAPPER_H

#include "EnumMapper.h"
#include "GeneralEnum.h"

class BWilsonBasisMapper : public EnumMapperBase<BWilsonBasis, BWilsonBasisMapper> {
public:
    static const std::map<BWilsonBasis, std::string>& mapping();
    static const std::map<std::string, BWilsonBasis>& inverse_mapping();
};

#endif