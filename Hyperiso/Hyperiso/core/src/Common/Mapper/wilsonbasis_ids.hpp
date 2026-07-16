#ifndef WILSON_BASIS_IDS_H
#define WILSON_BASIS_IDS_H

#include "generic_mapper.h"
#include "Map.h"

/**
 * @file WilsonBasisMapper.h
 * @brief Mapper for WilsonBasis <-> string identifiers.
 *
 * This is a thin wrapper around GenericMapperNoExt providing a dedicated
 * type-safe identifier (WilsonBasisId) and exposing the standard mapping
 * utilities (string lookup, enumeration listing, etc.).
 */
struct WilsonBasisTag {};

/** @brief Strongly typed identifier for WilsonBasis elements. */
using WilsonBasisId = IdOf<WilsonBasisTag>;

/**
 * @class WilsonBasisMapper
 * @brief Mapping interface for WilsonBasis.
 *
 * Inherits all functionality from GenericMapperNoExt:
 *   - enum -> string conversion,
 *   - case-insensitive string lookup -> SymbolId,
 *   - listing builtin and custom identifiers.
 */
class WilsonBasisMapper
: public GenericMapperNoExt<WilsonBasisTag, WilsonBasis, wilsonbasis_mapping>
{ public: using Base=GenericMapperNoExt<WilsonBasisTag,WilsonBasis,wilsonbasis_mapping>; using Base::str; };

#endif