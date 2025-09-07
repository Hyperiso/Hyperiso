#pragma once
#include "generic_mapper.hpp"
#include "Map.h"

struct ScaleTypeTag {};
using ScaleTypeId = IdOf<ScaleTypeTag>;

class ScaleTypeMapper
: public GenericMapperNoExt<ScaleTypeTag, ScaleType, scaletype_mapping>
{
public:
    using Base = GenericMapperNoExt<ScaleTypeTag, ScaleType, scaletype_mapping>;
    using Base::str;
    static std::string block(ScaleType t){ return Base::str(t); }
};