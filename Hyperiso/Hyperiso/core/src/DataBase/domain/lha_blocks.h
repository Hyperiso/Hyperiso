#ifndef HYPERISO_LHA_BLOCKS_H
#define HYPERISO_LHA_BLOCKS_H

#include <string>
#include <memory>
#include <vector>
#include <sstream>
#include <algorithm>

#include "lha_elements.h"
#include "LhaBlockPrototype.h"

/**
 * @class LhaBlock
 * @brief Represents an LHA block with multiple elements, based on a given prototype.
 */
class LhaBlock {
private:
    Prototype prototype; /**< The prototype defining the structure of the block. */
    std::vector<std::shared_ptr<AbstractElement>> entries; /**< List of elements contained in the block. */

public:
    /**
     * @brief Constructs an LhaBlock using a specified prototype.
     * @param prototype The prototype defining block structure.
     */
    inline explicit LhaBlock(const Prototype& prototype) : prototype(prototype) {}

    /**
     * @brief Reads data and populates the block with elements from given lines.
     * @param lines A vector of strings, each representing a line to be added as an element.
     */
    void readData(const std::vector<std::vector<std::string>>& lines);

    /**
     * @brief Checks whether an element is present in the block.
     * @param id Identifier of the element.
     * @return `true` if the given element is present in the block, `false` otherwise.
     */
    bool hasElement(const LhaID& id) const;

    /**
     * @brief Checks whether an element is present in the block.
     * @param id Parts of the identifier of the element.
     * @return `true` if the given element is present in the block, `false` otherwise.
     */
    bool hasElement(const std::vector<long>& id) const { return hasElement(LhaID(id)); };

    /**
     * @brief Retrieves an element by its identifier.
     * @param id Identifier of the element.
     * @return Pointer to the `AbstractElement` if found; otherwise, `nullptr`.
     */
    AbstractElement* get(const LhaID& id) const;

    /**
     * @brief Retrieves an element by its identifier.
     * @param id Parts of the identifier of the element.
     * @return Pointer to the `AbstractElement` if found; otherwise, `nullptr`.
     */
    AbstractElement* get(const std::vector<long>& id) const { return get(LhaID(id)); }

    /**
     * @brief Retrieves all entries in the block.
     * @return Pointer to a vector of unique pointers to `AbstractElement`.
     */
    const std::vector<std::shared_ptr<AbstractElement>>* getEntries() const;

    /**
     * @brief Adds a new element to the block from a line of data.
     * @param line The line data used to create the new element.
     */
    void addElement(const std::vector<std::string>& line);

    /**
     * @brief Retrieves the prototype associated with this block.
     * @return The `Prototype` instance.
     */
    inline Prototype getPrototype() { return this->prototype; };

    /**
     * @brief Converts the block and its entries to a string representation.
     * @return A string containing the block name and its entries.
     */
    std::string toString() const;

    /**
     * @brief Converts the block and its entries to a string representation.
     * @return A string containing the block name and its entries.
     */
    std::shared_ptr<Node> toDBNode() const;

    /**
     * @brief Destructor for LhaBlock, clears all entries.
     */
    inline ~LhaBlock() {entries.clear();}
};

#endif // HYPERISO_LHA_BLOCKS_H