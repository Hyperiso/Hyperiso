#ifndef HYPERISO_LHA_ELEMENTS_H
#define HYPERISO_LHA_ELEMENTS_H

#include <string>
#include <optional>
#include <vector>
#include <memory>
#include <map>
#include <sstream>
#include <iostream>
#include <limits>

#include "Logger.h"
#include "LhaID.h"
#include "LhaBlockPrototype.h"
#include "DBNode.h"

/**
 * @enum RenormalizationScheme
 * @brief Enumeration of possible renormalization schemes for LHA elements.
 *
 * The numeric values are interpreted directly from the corresponding
 * integer column in the LHA/FLHA block (when present).
 */
enum class RenormalizationScheme {
    POLE, MSBAR, DRBAR, ONE_S, KIN, INVARIANT, MOM, SMOM, NONE
};

/**
 * @class AbstractElement
 * @brief Abstract base class for a single element of an LHA/FLHA block.
 *
 * An AbstractElement represents one logical "entry" in a block:
 *   - it has a unique identifier (LhaID) built from the line indices,
 *   - it can expose an optional scale and binning,
 *   - it can be serialized to a string or to a DBNode.
 *
 * Concrete implementations (LhaElement<T>) provide the actual value type.
 */
class AbstractElement {
protected:
    const LhaID id;  /**< Unique identifier for the element. */

public:
    /**
     * @brief Constructs an AbstractElement with a given identifier.
     *
     * @param id Unique LhaID associated with this element.
     */
    inline explicit AbstractElement(const LhaID& id_) : id(id_) {}

    /**
     * @brief Returns the identifier of the element.
     *
     * The identifier encodes the integer indices of the line
     * (e.g. PDG codes, matrix indices), excluding value/scale/scheme/bin.
     *
     * @return LhaID of this element.
     */
    inline LhaID getId() const { return this->id; }

    /**
     * @brief Returns the binning associated to this element.
     *
     * The convention is defined by derived classes. For LhaElement<T>,
     * if no bin is defined, a dummy pair (-1.0, -1.0) is returned.
     *
     * @return Pair (low, high) bin edges.
     */
    virtual std::pair<double, double> getBinning() const = 0;

    /**
     * @brief Returns the renormalization scale of this element.
     *
     * For blocks without a scale, derived classes typically return 0.
     *
     * @return Scale value Q.
     */
    virtual double getScale() const = 0;

    /**
     * @brief Converts the element to a line-like string representation.
     *
     * The exact format is implementation-dependent; for LhaElement<T> it
     * typically contains:
     *   - the ID (indices),
     *   - the value,
     *   - optionally Q and the renormalization scheme.
     *
     * @return String representation of the element.
     */
    virtual std::string toString() const = 0;

    /**
     * @brief Converts the element into a DBNode for database/storage use.
     *
     * The node typically stores:
     *   - "central_value"           : the value,
     *   - "scale"                   : Q (if available),
     *   - "renormalization_scheme"  : scheme (if available),
     *   - "bin_low", "bin_high"     : bin edges (if available).
     *
     * @return Shared pointer to a newly created Node.
     */
    virtual std::shared_ptr<Node> toDBNode() const = 0;

    /**
     * @brief Virtual destructor.
     */
    virtual ~AbstractElement() = default;

};

/**
 * @class LhaElement
 * @brief Template concrete implementation of an LHA/FLHA block element.
 *
 * LhaElement<T> represents an element whose numerical content is of type T
 * (e.g. double or std::string). It is constructed from:
 *   - a Prototype describing the block structure (column indices),
 *   - a parsed line (vector of strings).
 *
 * From these, it:
 *   - determines the identifier (LhaID) by inspecting columns other than
 *     value/scale/scheme/bin,
 *   - extracts the value of type T,
 *   - optionally extracts:
 *       * the renormalization scale Q,
 *       * the renormalization scheme,
 *       * a binning (low, high).
 *
 * @tparam T Value type stored by the element (e.g. double or std::string).
 */
template <typename T>
class LhaElement : public AbstractElement {
private:
    T value;                                        /**< Value of the element. */
    std::optional<RenormalizationScheme> rScheme;   /**< Optional renormalization scheme. */
    std::optional<double> Q;                        /**< Optional scale value. */
    std::optional<std::pair<double, double>> bin;   /**< Optional energy binning.  */

    /**
     * @brief Encodes a unique LhaID for this line based on the prototype and raw tokens.
     *
     * All tokens that are not:
     *   - the value column,
     *   - the scale column,
     *   - the scheme column,
     *   - the bin columns (low, high),
     *   - the global scale (if globalScale is true and index 0),
     * and which can be interpreted as integers, are collected into a vector
     * of longs that forms the underlying LhaID.
     *
     * @param prototype Block prototype describing column roles.
     * @param line      Parsed line (string tokens).
     * @return LhaID built from the remaining integer tokens.
     */
    LhaID encodeId(const Prototype& prototype, const std::vector<std::string>& line);

public:
    /**
     * @brief Constructs an LhaElement from a prototype and a line of data.
     *
     * The constructor:
     *   - normalizes the indices using the Prototype (value, scale, scheme, bin),
     *   - builds the element ID via encodeId(),
     *   - parses the scale Q (either global header scale or per-line),
     *   - parses the renormalization scheme (if present),
     *   - parses the binning (if present),
     *   - converts the value token to type T.
     *
     * @param prototype Description of the block structure.
     * @param line      Vector of tokens representing one line of the block.
     *
     * @throws std::runtime_error if indices are inconsistent or conversion fails.
     */
    LhaElement(const Prototype& prototype, const std::vector<std::string>& line);

    /**
     * @brief Returns the stored value of the element.
     *
     * @return Value of type T.
     */
    inline T getValue() const { return this->value; }

    /**
     * @brief Returns the renormalization scheme of the element.
     *
     * If no scheme is defined, RenormalizationScheme::NONE is returned.
     *
     * @return RenormalizationScheme of this element.
     */
    inline RenormalizationScheme getScheme() const { return this->rScheme.has_value() ? this->rScheme.value() : RenormalizationScheme::NONE; }
    
    /**
     * @brief Returns the renormalization scale Q of the element.
     *
     * If no scale is defined, 0.0 is returned by convention.
     *
     * @return Scale value Q.
     */
    inline double getScale() const override { return this->Q.has_value() ? this->Q.value() : 0; }

    /**
     * @brief Returns the binning associated with the element.
     *
     * If no binning is defined, the default pair (-1.0, -1.0) is returned.
     *
     * @return Pair (low, high) bin edges.
     */
    inline std::pair<double, double> getBinning() const override { return this->bin.value_or(std::pair(-1.0, -1.0)); }
    
    /**
     * @brief Returns a line-like string representation of the element.
     *
     * The exact format is:
     *   "<id>\t<value>[\t<scale>][\t<scheme>]\n"
     *
     * where the optional parts are printed only if available.
     *
     * @return String representation of this element.
     */
    std::string toString() const override;

    /**
     * @brief Converts the element to a DBNode representation.
     *
     * The returned Node contains at least:
     *   - "central_value" : T
     * and optionally:
     *   - "scale"                 : double,
     *   - "renormalization_scheme": int,
     *   - "bin_low", "bin_high"   : double.
     *
     * @return Shared pointer to a new Node containing the element data.
     */
    std::shared_ptr<Node> toDBNode() const override;
    
};

/**
 * @class LhaElementFactory
 * @brief Factory for creating LHA/FLHA elements from a prototype and a line.
 *
 * The factory is responsible for choosing the appropriate template
 * instantiation of LhaElement<T> based on the block prototype. For example:
 *   - "FCINFO", "FMODSEL" and "SPINFO" are treated as string-valued blocks,
 *   - most other blocks are treated as double-valued.
 */
class LhaElementFactory {
    public:
        /**
         * @brief Creates an element of the appropriate type for a given line.
         *
         * The decision on the value type (e.g. double vs std::string) is based
         * on the prototype's blockName. The created element is always returned
         * via a shared pointer to AbstractElement.
         *
         * @param prototype Description of the LHA/FLHA block.
         * @param line      Parsed line tokens for a single entry.
         * @return Shared pointer to the created AbstractElement.
         *
         * @throws std::runtime_error If the line is inconsistent with the prototype.
         */
        static std::shared_ptr<AbstractElement> createElement(const Prototype& prototype, const std::vector<std::string>& line);
    };

#endif // HYPERISO_LHA_ELEMENTS_H
