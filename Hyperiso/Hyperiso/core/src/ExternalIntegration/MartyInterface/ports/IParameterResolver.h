#ifndef IPARAMETER_RESOLVER_H
#define IPARAMETER_RESOLVER_H

#include <unordered_map>
#include <string>
#include <vector>

#include "Extractor.h"
#include "LhaID.h"

/**
 * @file IParameterResolver.h
 * @brief Declares the interface for resolving parameter names to blocks and LHA IDs.
 *
 * The ::IParameterResolver interface is responsible for translating
 * MARTY-exposed parameter names into Hyperiso-specific identifiers
 * (::ResolvedParam) based on mapping databases.
 */

/**
 * @struct ResolvedParam
 * @ingroup MappingModule
 * @brief Represents a resolved parameter in Hyperiso format.
 *
 * This structure holds:
 *  - @ref block : the SLHA block name,
 *  - @ref code : the ::LhaID (indices within the block),
 *  - @ref is_complex : whether the parameter is complex,
 *  - @ref is_bsm : whether the parameter belongs to a BSM model.
 */
struct ResolvedParam {
    std::string block;  ///< SLHA block name.
    LhaID code;         ///< LHA ID within the block.
    bool is_complex;    ///< True if the parameter is complex-valued.
    bool is_bsm;        ///< True if the parameter is BSM (not SM).
};

/**
 * @class IParameterResolver
 * @ingroup MappingModule
 * @brief Interface for resolving external parameter names into ResolvedParam.
 *
 * Implementations typically combine multiple mapping databases (SM + model)
 * and apply model-dependent rules to decide whether a parameter is SM or BSM.
 */
class IParameterResolver {
public:
    virtual ~IParameterResolver() = default;

    /**
     * @brief Resolves a list of extracted parameters.
     *
     * For each entry in @p params, this method determines:
     *  - its block name,
     *  - LHA index,
     *  - complexity (real/complex),
     *  - SM vs BSM assignment.
     *
     * @param params    List of parameters as extracted from MARTY code.
     * @param modelIsSM True if the current calculation is SM-only.
     * @return Map from parameter name to ::ResolvedParam.
     */
    virtual std::unordered_map<std::string, ResolvedParam>
    resolve(const std::vector<Extractor::Parameter>& params,
            bool modelIsSM) const = 0;
    
    /**
     * @brief Clones the resolver.
     *
     * Used by higher-level components (::Interpreter) to maintain
     * independent copies without knowing the concrete type.
     *
     * @return A new heap-allocated copy of this resolver.
     */
    virtual std::unique_ptr<IParameterResolver> clone() const = 0;
};

#endif  // IPARAMETER_RESOLVER_H
