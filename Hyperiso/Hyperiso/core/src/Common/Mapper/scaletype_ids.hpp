#ifndef SCALE_TYPE_IDS_H
#define SCALE_TYPE_IDS_H

#include "generic_mapper.h"
#include "Map.h"

/**
 * @file ScaleTypeMapper.h
 * @brief Mapper for ScaleType <-> string identifiers.
 *
 * Provides utilities for converting ScaleType to text (e.g. "EW_SCALE",
 * "B_SCALE") and performing case-insensitive name lookup.
 */
struct ScaleTypeTag {};

/** @brief Strongly typed identifier for ScaleType elements. */
using ScaleTypeId = IdOf<ScaleTypeTag>;

class ScaleTypeMapper
: public GenericMapperNoExt<ScaleTypeTag, ScaleType, scaletype_mapping>
{
public:
    using Base = GenericMapperNoExt<ScaleTypeTag, ScaleType, scaletype_mapping>;
    using Base::str;

    /// Returns the canonical string block name associated with a ScaleType.
    static std::string block(ScaleType t){ return Base::str(t); }
};

#endif