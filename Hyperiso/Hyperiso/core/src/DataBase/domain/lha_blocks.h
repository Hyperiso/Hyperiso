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
 * @file lha_blocks.h
 * @brief Representation of a single LHA/FLHA block and its entries.
 *
 * An LhaBlock groups together all elements belonging to a given block
 * (e.g. "SMINPUTS", "FWCOEF", "FMASS"), using a Prototype to interpret
 * the structure of each line. Elements are stored as AbstractElement
 * instances created via LhaElementFactory.
 */

/**
 * @class LhaBlock
 * @brief Represents an LHA/FLHA block with multiple typed elements.
 *
 * An LhaBlock:
 *   - is defined by a Prototype that describes the block layout
 *     (value column, scale column, etc.),
 *   - owns a list of entries (AbstractElement) created from parsed lines,
 *   - provides lookup by LhaID,
 *   - can be serialized to a human-readable string or to a DBNode.
 */
class LhaBlock {
private:
    Prototype prototype; /**< The prototype defining the structure of the block. */
    std::vector<std::shared_ptr<AbstractElement>> entries; /**< List of elements contained in the block. */

public:
    /**
     * @brief Constructs an LhaBlock using a specified prototype.
     *
     * The block starts empty; elements are later added via readData()
     * or addElement().
     *
     * @param prototype Prototype defining the block structure.
     */
    inline explicit LhaBlock(const Prototype& prototype_) : prototype(prototype_) {}

    /**
     * @brief Reads data and populates the block with elements from raw lines.
     *
     * Each non-empty line is passed to LhaElementFactory::createElement()
     * together with this block's prototype to create the appropriate
     * LhaElement<T>, which is then stored in the block.
     *
     * @param lines Vector of tokenized lines; each inner vector represents
     *              one block entry (split into string fields).
     */
    void readData(const std::vector<std::vector<std::string>>& lines);

    /**
     * @brief Checks whether an element is present in the block.
     *
     * This performs a linear search over the stored entries and compares
     * their LhaID with @p id.
     *
     * @param id Identifier of the element.
     * @return true if an element with this id exists, false otherwise.
     */
    bool hasElement(const LhaID& id) const;

    /**
     * @brief Convenience overload: checks presence from raw ID parts.
     *
     * Constructs an LhaID from @p id and forwards to hasElement(const LhaID&).
     *
     * @param id Components of the identifier.
     * @return true if an element with this id exists, false otherwise.
     */
    bool hasElement(const std::vector<long>& id) const { return hasElement(LhaID(id)); };

    /**
     * @brief Retrieves an element by its identifier.
     *
     * Performs a linear search and returns a non-owning pointer to the
     * stored element if found.
     *
     * @param id Identifier of the element.
     * @return Pointer to AbstractElement if found, nullptr otherwise.
     */
    AbstractElement* get(const LhaID& id) const;

    /**
     * @brief Convenience overload: retrieves an element from raw ID parts.
     *
     * Constructs an LhaID from @p id and forwards to get(const LhaID&).
     *
     * @param id Components of the identifier.
     * @return Pointer to AbstractElement if found, nullptr otherwise.
     */
    AbstractElement* get(const std::vector<long>& id) const { return get(LhaID(id)); }

    /**
     * @brief Returns a read-only view of all entries in the block.
     *
     * The lifetime of the returned pointer is tied to the LhaBlock instance.
     * It must not be deleted by the caller.
     *
     * @return Pointer to the internal vector of shared_ptr<AbstractElement>.
     */
    const std::vector<std::shared_ptr<AbstractElement>>* getEntries() const;

    /**
     * @brief Adds a new element to the block from a line of data.
     *
     * The line is interpreted according to this block's Prototype and
     * passed to LhaElementFactory::createElement().
     *
     * @param line Tokens representing one block entry.
     *
     * @throws std::runtime_error if the line is incompatible with the prototype.
     */
    void addElement(const std::vector<std::string>& line);

    /**
     * @brief Returns the Prototype associated with this block.
     *
     * The prototype is returned by value (copy).
     *
     * @return Prototype describing this block.
     */
    inline Prototype getPrototype() { return this->prototype; };

    /**
     * @brief Serializes the block and all its entries to a string.
     *
     * The format is:
     *   "Block <BLOCKNAME>:\n"
     *   "<entry1>\n"
     *   "<entry2>\n"
     *   ...
     * where each entry line is produced by AbstractElement::toString().
     *
     * @return Human-readable string representation of the block.
     */
    std::string toString() const;

    /**
     * @brief Converts the block and its entries to a DBNode representation.
     *
     * The resulting node has the following structure:
     *   - If the prototype has globalScale == true, the node stores:
     *       "scale" : <double>    (taken from the first entry)
     *   - Under a child block named prototype.blockName, each entry is
     *     stored as:
     *       "<id>" : <Node>       (Node produced by entry->toDBNode())
     *
     * @return Shared pointer to the root Node representing this block.
     */
    std::shared_ptr<Node> toDBNode() const;

    /**
     * @brief Destructor for LhaBlock.
     *
     * Clears the entries vector; shared_ptr semantics ensure any remaining
     * elements are properly destroyed when no longer referenced.
     */
    inline ~LhaBlock() {entries.clear();}
};

#endif // HYPERISO_LHA_BLOCKS_H