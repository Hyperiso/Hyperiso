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
     * @brief Retrieves the energy scale of the element.
     * @return The energy binning of the element, if defined.
     */
    virtual std::pair<double, double> getBinning() const = 0;

    /**
     * @brief Retrieves the energy binning of the element.
     * @return The energy scale of the element, if defined.
     */
    virtual double getScale() const = 0;

    /**
     * @brief Converts the element to a string representation.
     * @return String representing the element.
     */
    virtual std::string toString() const = 0;

    /**
     * @brief Converts the element to a DBNode representation.
     * @return DBNode representing the element.
     */
    virtual std::shared_ptr<Node> toDBNode() const = 0;

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
    std::optional<std::pair<double, double>> bin;         /**< Optional energy binning.  */

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
    inline double getScale() const override { return this->Q.has_value() ? this->Q.value() : 0; }

    /**
     * @brief Retrieves the scale associated with the element.
     * @return The scale value or `0` if not set.
     */
    inline std::pair<double, double> getBinning() const override { return this->bin.value_or(std::pair(-1.0, -1.0)); }
    /**
     * @brief Converts the element to a string representation.
     * @return String representing the element.
     */
    std::string toString() const override;

    /**
     * @brief Converts the element to a DBNode representation.
     * @return DBNode representing the element.
     */
    std::shared_ptr<Node> toDBNode() const override;
    
};

/**
 * @class LhaElementFactory
 * @brief Factory class for creating LHA elements.
 * 
 * Provides a static method for creating elements of the appropriate type based on the block and data line.
 */
class LhaElementFactory {
    public:
        /**
         * @brief Creates an LHA element based on the block and line data.
         * 
         * Determines the element type (e.g., `double` or `std::string`) based on the block name
         * and creates an appropriate instance.
         * 
         * @param block Pointer to the LhaBlock that will contain the new element.
         * @param line Vector of strings representing the line data.
         * @return Unique pointer to the created `AbstractElement`.
         */
        static std::shared_ptr<AbstractElement> createElement(const Prototype& prototype, const std::vector<std::string>& line);
    };

#endif // HYPERISO_LHA_ELEMENTS_H
