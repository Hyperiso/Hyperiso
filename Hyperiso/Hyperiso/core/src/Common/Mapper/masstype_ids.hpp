
#pragma once
#include "generic_mapper.hpp"
#include "Map.h"

struct MassTypeTag {};
using MassTypeId = IdOf<MassTypeTag>;

class MassTypeMapper
: public GenericMapperNoExt<MassTypeTag, MassType, masstype_mapping>
{ public: using Base=GenericMapperNoExt<MassTypeTag,MassType,masstype_mapping>; using Base::str; };
