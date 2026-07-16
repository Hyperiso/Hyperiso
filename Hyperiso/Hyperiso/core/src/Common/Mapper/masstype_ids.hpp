
#ifndef MASS_TYPE_IDS_H
#define MASS_TYPE_IDS_H

#include "generic_mapper.h"
#include "Map.h"

/**
 * @file MassTypeMapper.h
 * @brief Mapper for MassType <-> string identifiers.
 */
struct MassTypeTag {};
using MassTypeId = IdOf<MassTypeTag>;

class MassTypeMapper
: public GenericMapperNoExt<MassTypeTag, MassType, masstype_mapping>
{ public: using Base=GenericMapperNoExt<MassTypeTag,MassType,masstype_mapping>; using Base::str; };

#endif