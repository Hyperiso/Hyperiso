#ifndef SOURCES_VIEW_H
#define SOURCES_VIEW_H

#include <unordered_map>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <stdexcept>
#include <vector>
#include <utility>

#include "Block.h"
#include "Parameters.h"


/**
 * @brief Joins a list of strings with ", ".
 *
 * Utility function mainly used to format error messages with
 * "available keys" lists.
 *
 * @param v Vector of strings to join.
 * @return A single string of the form "a, b, c".
 */
std::string join(const std::vector<std::string>& v);

/**
 * @brief Convenience stringification for LhaID.
 *
 * Delegates to LhaID::to_string().
 *
 * @param id LhaID to stringify.
 * @return String representation of the ID.
 */
std::string to_string(const LhaID& id);

/**
 * @brief Convenience stringification for ParamId.
 *
 * Produces a human-readable identifier of the form:
 *   <ParameterType>::<block>::<code>
 * using ParameterTypeMapper::str for the type.
 *
 * @param id ParamId to stringify.
 * @return String representation of the ParamId.
 */
std::string to_string(const ParamId& id);

/**
 * @class BlockSrc
 * @brief Lightweight view over a set of source blocks.
 *
 * This class is a non-owning wrapper around a
 * `std::unordered_map<std::string, std::shared_ptr<Block>>`.
 *
 * It provides:
 *   - checked access to a block by name, with good error messages
 *   - checked access to parameters in those blocks
 *   - convenient `get_val` helpers from various ID forms
 *
 * It is typically used inside dependent-block update lambdas:
 *
 * @code
 * void recalc(const BlockSrc& src, std::shared_ptr<DependentBlock> self) {
 *     double alpha = src.get_val("SMINPUTS", 1);
 *     ...
 * }
 * @endcode
 */
class BlockSrc {
public:
    /**
     * @brief Constructs a BlockSrc view from an existing map.
     *
     * @param m   Map of block-name to Block pointers (not owned).
     * @param ctx Optional textual context added to error messages
     *            (e.g. name of the DependentBlock being updated).
     */
    explicit BlockSrc(const std::unordered_map<std::string, std::shared_ptr<Block>>& m,
                      std::string ctx = {});
    
    /**
     * @brief Checks whether a block with a given name is present.
     *
     * @param name Name of the block.
     * @return True if the underlying map contains this block.
     */
    bool has_block(std::string_view name) const;
    
    /**
     * @brief Returns a reference to a block by name.
     *
     * If the block is missing, throws std::invalid_argument with a message
     * mentioning the context (if provided) and the list of available blocks.
     *
     * @param name Name of the block (case-sensitive as stored in the map).
     * @return Const reference to shared_ptr<Block>.
     *
     * @throws std::invalid_argument if the block is not found.
     */
    const std::shared_ptr<Block>& block(std::string_view name) const;
    
    /**
     * @brief Retrieves a parameter from a given block by LhaID.
     *
     * If the block or parameter is missing, throws std::invalid_argument with
     * a detailed message including available IDs.
     *
     * @param blk  Name of the block.
     * @param code LhaID of the parameter.
     * @return Shared pointer to the Parameter.
     *
     * @throws std::invalid_argument if the block or parameter is not found.
     */
    std::shared_ptr<Parameter> get_param(std::string_view blk, const LhaID& code) const;
    
    /**
     * @brief Helper: get parameter using an integer code.
     *
     * Constructs LhaID(code) and forwards to the main overload.
     *
     * @param blk  Name of the block.
     * @param code Integer PDG-like code.
     * @return Shared pointer to the Parameter.
     */
    std::shared_ptr<Parameter> get_param(std::string_view blk, int code) const;

    /**
     * @brief Helper: get parameter using a pair of integers.
     *
     * Constructs an LhaID(code.first, code.second) and forwards.
     *
     * @param blk  Name of the block.
     * @param code Pair of indices.
     * @return Shared pointer to the Parameter.
     */
    std::shared_ptr<Parameter> get_param(std::string_view blk, std::pair<int,int> code) const;
    
    /**
     * @brief Retrieves the value of a parameter given as an initializer_list.
     *
     * Example:
     * @code
     * auto v = src.get_val("MASS", {25});        // H mass
     * auto v2 = src.get_val("YU", {3,3});       // Yukawa entry
     * @endcode
     *
     * @param blk  Name of the block.
     * @param code List of integer components used to build an LhaID.
     * @return The current parameter value (scalar_t).
     *
     * @throws std::invalid_argument if the block or parameter is missing.
     */
    scalar_t get_val(std::string_view blk, std::initializer_list<long> code) const;
    
    /**
     * @brief Retrieves the value of a parameter given as an initializer_list.
     *
     * Example:
     * @code
     * auto v = src.get_val("MASS", {25});        // H mass
     * auto v2 = src.get_val("YU", {3,3});       // Yukawa entry
     * @endcode
     *
     * @param blk  Name of the block.
     * @param code List of positive integer components used to build an LhaID.
     * @return The current parameter value (scalar_t).
     *
     * @throws std::invalid_argument if the block or parameter is missing.
     */
    scalar_t get_val(std::string_view blk,
                           std::initializer_list<std::size_t> code) const;
    
    /**
     * @brief Retrieves the value of a parameter given as an initializer_list.
     *
     * Example:
     * @code
     * auto v = src.get_val("MASS", {25});        // H mass
     * auto v2 = src.get_val("YU", {3,3});       // Yukawa entry
     * @endcode
     *
     * @param blk  Name of the block.
     * @param code List of positive integer components used to build an LhaID.
     * @return The current parameter value (scalar_t).
     *
     * @throws std::invalid_argument if the block or parameter is missing.
     */
    scalar_t get_val(std::string_view blk, std::initializer_list<int> code) const;

    /**
     * @brief Retrieves the value of a parameter identified by an LhaID.
     *
     * @param blk  Name of the block.
     * @param code LhaID of the parameter.
     * @return The current parameter value as double.
     */
    double get_val(std::string_view blk, const LhaID& code) const;

    /**
     * @brief Helper: get value from an integer code.
     *
     * @param blk  Name of the block.
     * @param code Integer PDG-like code.
     * @return The current parameter value as double.
     */
    double get_val(std::string_view blk, int code) const;

    /**
     * @brief Helper: get value from a pair of integers.
     *
     * @param blk  Name of the block.
     * @param code Pair of indices.
     * @return The current parameter value as double.
     */
    double get_val(std::string_view blk, std::pair<int,int> code) const;
    
    /**
     * @brief Returns the underlying map reference.
     *
     * This is useful if the caller needs to iterate or inspect the raw
     * collection of blocks.
     *
     * @return Const reference to the underlying map.
     */
    const std::unordered_map<std::string, std::shared_ptr<Block>>& raw() const;

private:
    const std::unordered_map<std::string, std::shared_ptr<Block>>* m_;  ///< Non-owning pointer to block map.
    std::string ctx_;                                                   ///< Optional context for error messages.
};

/**
 * @class ParamSrc
 * @brief Lightweight view over a set of source parameters keyed by ParamId.
 *
 * This class wraps a map
 * `std::unordered_map<ParamId, std::shared_ptr<Parameter>>` without owning it,
 * and provides:
 *
 *   - presence checks,
 *   - safe retrieval with detailed error messages,
 *   - convenient `get_val` helpers accepting (type, block, code...).
 *
 * It is particularly useful for dependency graphs where derived parameters
 * are functions of other parameters identified by full ParamId.
 */
class ParamSrc {
public:
    /**
     * @brief Constructs a ParamSrc view from an existing map.
     *
     * @param m   Map of ParamId to Parameter pointers (not owned).
     * @param ctx Optional context string used in error messages.
     */
    explicit ParamSrc(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& m,
                      std::string ctx = {});
    
    /// @brief Default constructor: creates an empty, unbound view.
    ParamSrc() = default;

    /**
     * @brief Checks whether a parameter with the given ID exists.
     *
     * @param id ParamId to test.
     * @return True if the map contains this ID.
     */
    bool has(const ParamId& id) const;
    
    /**
     * @brief Returns a parameter reference or throws if missing.
     *
     * On failure, throws std::invalid_argument with a message that includes
     * the context (if any) and the list of available ParamIds.
     *
     * @param id ParamId to retrieve.
     * @return Const reference to the shared_ptr<Parameter>.
     *
     * @throws std::invalid_argument if the parameter is not found.
     */
    const std::shared_ptr<Parameter>& require(const ParamId& id) const;
    
    /**
     * @brief Convenience wrapper for require().
     *
     * @param id ParamId to retrieve.
     * @return Shared pointer to the Parameter.
     *
     * @throws std::invalid_argument if the parameter is not found.
     */
    std::shared_ptr<Parameter> get_param(const ParamId& id) const;

    /**
     * @brief Retrieves the current value of a parameter.
     *
     * @param id ParamId to retrieve.
     * @return The parameter value (scalar_t).
     */
    scalar_t get_val(const ParamId& id) const;
    
    /**
     * @brief Helper: get value from typed (ParameterType, block, code-list).
     *
     * Constructs a ParamId and forwards to the main get_val().
     *
     * @param t     ParameterType (owner / category).
     * @param block Block name.
     * @param code  List of integer components used to build LhaID.
     * @return The parameter value.
     */
    scalar_t get_val(ParameterType t, std::string_view block, std::initializer_list<long> code) const;

    /**
     * @brief Helper: get value from typed (ParameterType, block, LhaID).
     */
    scalar_t get_val(ParameterType t, std::string_view block, const LhaID& code) const;

    /**
     * @brief Helper: get value from typed (ParameterType, block, int).
     */
    scalar_t get_val(ParameterType t, std::string_view block, int code) const;

    /**
     * @brief Helper: get value from typed (ParameterType, block, pair<int,int>).
     */
    scalar_t get_val(ParameterType t, std::string_view block, std::pair<int,int> code) const;
    
    /**
     * @brief Access to the underlying parameter map.
     *
     * @return Const reference to the internal map.
     */
    const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& raw() const;

private:
    const std::unordered_map<ParamId, std::shared_ptr<Parameter>>* m_;  ///< Non-owning pointer to parameter map.
    std::string ctx_;                                                   ///< Optional context for error messages.
};

/**
 * @brief Helper function to build a BlockSrc view from a block map.
 *
 * @param m   Map of block name to Block pointer.
 * @param ctx Optional context string for error reporting.
 * @return A BlockSrc instance referencing @p m.
 */
BlockSrc as_block_src(const std::unordered_map<std::string, std::shared_ptr<Block>>& m, std::string ctx = {});

/**
 * @brief Helper function to build a ParamSrc view from a parameter map.
 *
 * @param m   Map of ParamId to Parameter pointer.
 * @param ctx Optional context string for error reporting.
 * @return A ParamSrc instance referencing @p m.
 */
ParamSrc as_param_src(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& m, std::string ctx = {});

#endif