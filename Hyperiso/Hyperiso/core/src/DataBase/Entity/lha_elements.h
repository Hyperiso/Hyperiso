#ifndef HYPERISO_LHA_ELEMENTS_H
#define HYPERISO_LHA_ELEMENTS_H

#include <string>
#include <optional>
#include <vector>
#include <memory>
#include <map>
#include "Logger.h"

/**
 * @struct Prototype
 * @brief Represents the structure of an LHA block, specifying columns for values, scales, and renormalization groups.
 */
struct Prototype {
    std::string blockName;      /**< Block name, case insensitive. */
    int itemCount {2};          /**< Number of columns in the block. */
    int valueIdx {1};           /**< Column index for value. */
    int scaleIdx {-1};          /**< Column index for scale, -1 if scale-independent. */
    int rgIdx {-1};             /**< Column index for renormalization group, -1 if irrelevant. */
    bool globalScale {false};   /**< Indicates if the block uses a global scale (Q= in header). */
};

/**
 * @struct LhaID
 * @brief Represents an identifier of a LHA element, possibly containing several sub-ids 
 */
struct LhaID {
    std::vector<long> parts;     /**< Collection of sub-ids. */

    /**
     * @brief Constructs a LhaID with specified sub-ids
     * @param parts List of sub-ids of the element
     */
    LhaID(const std::vector<long>& parts) : parts(std::move(parts)) {}

    /**
     * @brief Constructs a LhaID with a single identifier
     * @param id Identifier of the element
     */
    LhaID(long id) : parts({id}) {}

    
    /**
     * @brief Allows for implicit conversion of a trivial LhaID to an integer 
     */
    operator long() const {
        if (this->parts.size() > 1) {
            LOG_WARN("Casting nontrivial LhaID to int discards information.");
        }
        return this->parts.at(0);
    };

    inline friend bool operator==(const LhaID& lhs, const LhaID& rhs) { return lhs.parts == rhs.parts; };
    inline friend bool operator!=(const LhaID& lhs, const LhaID& rhs) { return !(lhs == rhs); };

    friend std::ostream& operator<<(std::ostream&, const LhaID&);
};

/**
 * @enum RenormalizationScheme
 * @brief Enumeration of different renormalization schemes for LHA elements.
 */
enum class RenormalizationScheme {
    POLE, MSBAR, DRBAR, ONE_S, KIN, INVARIANT, MOM, SMOM, NONE
};

/**
 * @class AbstractElement
 * @brief Abstract base class for elements within an LHA block.
 * 
 * Provides a common interface for all elements, including an identifier and
 * a method to convert the element to a string.
 */
class AbstractElement {
protected:
    const LhaID id;  /**< Unique identifier for the element. */

public:
    /**
     * @brief Constructs an AbstractElement with a specified block and identifier.
     * @param id Unique identifier for the element. 
     */
    inline explicit AbstractElement(const LhaID& id) : id(id) {}

    /**
     * @brief Retrieves the identifier of the element.
     * @return String containing the element's ID.
     */
    inline LhaID getId() const { return this->id; }

    /**
     * @brief Converts the element to a string representation.
     * @return String representing the element.
     */
    virtual std::string toString() const = 0;

    /**
     * @brief Virtual destructor for AbstractElement.
     */
    virtual ~AbstractElement() = default;

};

/**
 * @class LhaElement
 * @brief Template class for elements within an LHA block.
 * 
 * Represents an element with a specific type and optional renormalization scheme and scale.
 * 
 * @tparam T Type of the element's value (e.g., `double` or `std::string`).
 */
template <typename T>
class LhaElement : public AbstractElement {
private:
    T value;                                      /**< Value of the element. */
    std::optional<RenormalizationScheme> rScheme; /**< Optional renormalization scheme. */
    std::optional<double> Q;                      /**< Optional scale value. */

    /**
     * @brief Encodes a unique ID for the element based on its data.
     * @param block Pointer to the LhaBlock containing the element.
     * @param line Vector of strings representing the line data for the element.
     * @return Encoded ID as an initialized LhaID.
     */
    LhaID encodeId(const Prototype& prototype, const std::vector<std::string>& line);

public:
    /**
     * @brief Constructs an LhaElement with specified block and data line.
     * @param prototype Prototype of the LhaBlock containing the element.
     * @param line Vector of strings representing the line data.
     */
    LhaElement(const Prototype& prototype, const std::vector<std::string>& line);

    /**
     * @brief Retrieves the value of the element.
     * @return The value of the element.
     */
    inline T getValue() const { return this->value; }

    /**
     * @brief Retrieves the renormalization scheme of the element.
     * @return The renormalization scheme or `RenormalizationScheme::NONE` if not set.
     */
    inline RenormalizationScheme getScheme() const { return this->rScheme.has_value() ? this->rScheme.value() : RenormalizationScheme::NONE; }
    
    /**
     * @brief Retrieves the scale associated with the element.
     * @return The scale value or `0` if not set.
     */
    inline double getScale() const { return this->Q.has_value() ? this->Q.value() : 0; }

    /**
     * @brief Converts the element to a string representation.
     * @return String representing the element.
     */
    std::string toString() const override;
    
};

#endif // HYPERISO_LHA_ELEMENTS_H
