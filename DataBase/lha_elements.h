#ifndef HYPERISO_LHA_ELEMENTS_H
#define HYPERISO_LHA_ELEMENTS_H

#include <string>
#include <optional>
#include <vector>
#include <memory>
#include <map>
#include "Logger.h"

class LhaBlock;

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
    const std::string id;  /**< Unique identifier for the element. */
    const LhaBlock* block; /**< Pointer to the LhaBlock that contains this element. */

public:
    /**
     * @brief Constructs an AbstractElement with a specified block and identifier.
     * @param block Pointer to the LhaBlock containing the element.
     * @param id Unique identifier for the element.
     */
    inline explicit AbstractElement(LhaBlock* block, const std::string& id) : id(id),block(block)  {}

    /**
     * @brief Retrieves the identifier of the element.
     * @return String containing the element's ID.
     */
    inline std::string getId() const { return this->id; }

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
     * @return Encoded ID as a string.
     */
    std::string encodeId(LhaBlock* block, const std::vector<std::string>& line);

public:
    /**
     * @brief Constructs an LhaElement with specified block and data line.
     * @param block Pointer to the LhaBlock containing the element.
     * @param line Vector of strings representing the line data.
     */
    LhaElement(LhaBlock* block, const std::vector<std::string>& line);

    /**
     * @brief Retrieves the LhaBlock containing this element.
     * @return Pointer to the containing LhaBlock.
     */
    inline LhaBlock* getBlock() const { return this->block; }

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
    static std::unique_ptr<AbstractElement> createElement(LhaBlock* block, const std::vector<std::string>& line);
};

#endif // HYPERISO_LHA_ELEMENTS_H
