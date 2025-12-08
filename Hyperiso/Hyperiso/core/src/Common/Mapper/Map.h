#ifndef MAP_H
#define MAP_H

#include <map>
#include <string>
#include <vector>

#include "GeneralEnum.h"
#include "LhaID.h"

/**
 * @file Map.h
 * @brief Static mappings between enum classes and string (or FLHA) identifiers.
 *
 * This header declares a set of helper functions returning references
 * to static maps. These maps provide a bidirectional bridge between
 * internal enum types (e.g. Observables, ParameterType, Model, ...)
 * and their string / FLHA / auxiliary representations used in input
 * or output formats.
 *
 * All returned maps are initialized once on first use via function-local
 * static storage and then reused.
 */

/**
 * @brief Returns the mapping between Observables and their string names.
 */
const std::map<Observables, std::string>& observable_mapping();

/**
 * @brief Returns the mapping between Observables and their FLHA identifiers.
 */
const std::map<Observables, LhaID>&       observable_flha_mapping();

/**
 * @brief Returns the mapping between QCDOrder and its string representation.
 */
const std::map<QCDOrder,        std::string>& order_mapping();

/**
 * @brief Returns the mapping between WGroup and its string representation.
 */
const std::map<WGroup,          std::string>& group_mapping();

/**
 * @brief Returns the mapping between WCoef and its string name.
 */
const std::map<WCoef,           std::string>& wcoef_mapping();

/**
 * @brief Returns the mapping between WCoef and its FLHA (block, index) pair.
 */
const std::map<WCoef, std::pair<int,int>>&    wcoef_flha_mapping();

/**
 * @brief Returns the mapping between ParameterType and its string name.
 */
const std::map<ParameterType,   std::string>& parametertype_mapping();

/**
 * @brief Returns the mapping between Model and its string name.
 */
const std::map<Model,           std::string>& model_mapping();

/**
 * @brief Returns the mapping between WilsonBasis and its string name.
 */
const std::map<WilsonBasis,     std::string>& wilsonbasis_mapping();

/**
 * @brief Returns the mapping between ContributionType and its string name.
 */
const std::map<ContributionType,std::string>& contributiontype_mapping();

/**
 * @brief Returns the mapping between MassType and its string name.
 */
const std::map<MassType,        std::string>& masstype_mapping();

/**
 * @brief Returns the mapping between ScaleType and its string name.
 */
const std::map<ScaleType,       std::string>& scaletype_mapping();

/**
 * @brief Returns the mapping between Decays and their string names.
 */
const std::map<Decays,          std::string>& decays_mapping();

/**
 * @brief Returns the mapping between Decays and the list of Observables
 *        associated with each decay channel.
 */
const std::map<Decays, std::vector<Observables>>& decay_observable_mapping();

#endif // MAP_H
