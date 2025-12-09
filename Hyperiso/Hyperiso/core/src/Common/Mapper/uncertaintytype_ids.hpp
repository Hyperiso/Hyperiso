#ifndef UNCERTAINTY_TYPE_IDS_H
#define UNCERTAINTY_TYPE_IDS_H

#include "generic_mapper.h"
#include "Map.h"
#include "GeneralEnum.h"

/**
 * @file uncertaintytype_ids.hpp
 * @brief Mapping and helpers for uncertainty types.
 *
 * This header defines:
 *   - UncTag / UncertaintyTypeId: strongly typed identifiers for uncertainty types,
 *   - uncertaintytype_mapping(): builtin enum ↔ string mapping,
 *   - UncertaintyTypeMapper: mapper on top of GenericMapperNoExt with an
 *     additional mapping to DataType.
 */

/** @brief Tag type for uncertainty-type identifiers. */
struct UncTag {};

/** @brief Strongly typed identifier for UncertaintyType values. */
using UncertaintyTypeId = IdOf<UncTag>;

/**
 * @brief Builtin mapping UncertaintyType -> string name.
 *
 * The names are intended for use in user-facing I/O or configuration:
 *   - STAT     -> "STATISTICAL"
 *   - SYST     -> "SYSTEMATICS"
 *   - COMBINED -> "COMBINED"
 */
inline const std::map<UncertaintyType, std::string>& uncertaintytype_mapping() {
    static const std::map<UncertaintyType, std::string> m = {
        {UncertaintyType::STAT,     "STATISTICAL"},
        {UncertaintyType::SYST,     "SYSTEMATICS"},
        {UncertaintyType::COMBINED, "COMBINED"},
    };
    return m;
}

/**
 * @class UncertaintyTypeMapper
 * @brief Mapper for UncertaintyType <-> UncertaintyTypeId <-> string.
 *
 * This is a thin wrapper around GenericMapperNoExt with:
 *   - Tag   = UncTag
 *   - EnumT = UncertaintyType
 *   - MapFn = uncertaintytype_mapping
 *
 * It also exposes a separate mapping from UncertaintyType to DataType
 * (e.g. STD_STAT, STD_SYST, ...).
 */
class UncertaintyTypeMapper
: public GenericMapperNoExt<UncTag, UncertaintyType, uncertaintytype_mapping>
{
public:
    using Base = GenericMapperNoExt<UncTag, UncertaintyType, uncertaintytype_mapping>;
    using Base::str;

    /**
     * @brief Returns the builtin mapping UncertaintyType -> DataType.
     *
     * This is typically used when interfacing to the internal data
     * representation of uncertainties in datasets.
     */
    static const std::map<UncertaintyType, DataType>& data_type_mapping(){
        static const std::map<UncertaintyType, DataType> m = {
            {UncertaintyType::STAT,     DataType::STD_STAT},
            {UncertaintyType::SYST,     DataType::STD_SYST},
            {UncertaintyType::COMBINED, DataType::STD_COMBINED},
        };
        return m;
    }

    /**
     * @brief Returns the DataType corresponding to a given uncertainty type.
     *
     * @param t UncertaintyType value.
     * @return Matching DataType.
     *
     * @throws std::out_of_range if @p t is not in the mapping.
     */
    static DataType d_type(UncertaintyType t){ return data_type_mapping().at(t); }
};

#endif