#ifndef HYPERISO_LHA_BLOCKS_H
#define HYPERISO_LHA_BLOCKS_H

#include <string>
#include <memory>
#include <vector>
#include <sstream>
#include <algorithm>

#include "lha_elements.h"

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

// SLHA Block prototypes
const Prototype MODSEL = Prototype{"MODSEL"};
const Prototype SMINPUTS = Prototype{"SMINPUTS"};
const Prototype VCKMIN = Prototype{"VCKMIN"};   /**< SLHA2 */
const Prototype UPMNSIN = Prototype{"UPMNSIN"}; /**< SLHA2 */
const Prototype MINPAR = Prototype{"MINPAR"};
const Prototype EXTPAR = Prototype{"EXTPAR"};
const Prototype MASS = Prototype{"MASS"};
const Prototype NMIX = Prototype{"NMIX", 3, 2};
const Prototype UMIX = Prototype{"UMIX", 3, 2};
const Prototype VMIX = Prototype{"VMIX", 3, 2};
const Prototype STOPMIX = Prototype{"STOPMIX", 3, 2};
const Prototype SBOTMIX = Prototype{"SBOTMIX", 3, 2};
const Prototype STAUMIX = Prototype{"STAUMIX", 3, 2};
const Prototype ALPHA = Prototype{"ALPHA", 1, 0};
const Prototype HMIX = Prototype{"HMIX", 3, 2, -1, -1, true};
const Prototype GAUGE = Prototype{"GAUGE", 3, 2, -1, -1, true};
const Prototype MSOFT = Prototype{"MSOFT", 3, 2, -1, -1, true};
const Prototype AU = Prototype{"AU", 4, 3, -1, -1, true};
const Prototype AD = Prototype{"AD", 4, 3, -1, -1, true};
const Prototype AE = Prototype{"AE", 4, 3, -1, -1, true};
const Prototype YU = Prototype{"YU", 4, 3, -1, -1, true};
const Prototype YD = Prototype{"YD", 4, 3, -1, -1, true};
const Prototype YE = Prototype{"YE", 4, 3, -1, -1, true};
const Prototype SPINFO = Prototype{"SPINFO"};

const std::vector<Prototype> SLHA_BLOCKS = {MODSEL, SMINPUTS, VCKMIN, UPMNSIN, MINPAR, EXTPAR, MASS, NMIX, UMIX, VMIX, STOPMIX, SBOTMIX, STAUMIX, ALPHA, HMIX, GAUGE, MSOFT, AU, AD, AE, YU, YD, YE, SPINFO};

// FLHA Block prototypes
const Prototype FCINFO = Prototype{"FCINFO"};
const Prototype FMODSEL = Prototype{"FMODSEL"};
const Prototype FMASS = Prototype{"FMASS", 4, 1, 2, 3};
const Prototype FLIFE = Prototype{"FLIFE"};
const Prototype FCONST = Prototype{"FCONST", 5, 2, 3, 4};
const Prototype FCONSTRATIO = Prototype{"FCONSTRATIO", 7, 4, 5, 6};
const Prototype FBAG = Prototype{"FBAG", 5, 2, 3, 4};
const Prototype FWCOEF = Prototype{"FWCOEF", 6, 5, -1, -1, true};
const Prototype IMFWCOEF = Prototype{"IMFWCOEF", 6, 5, -1, -1, true};
const Prototype FOBS = Prototype{"FOBS", 9, 2, 3};
const Prototype FOBSERR = Prototype{"FOBSERR", 9, 2, 3};
const Prototype FOBSSM = Prototype{"FOBSSM", 9, 2, 3};
const Prototype FDIPOLE = Prototype{"FDIPOLE", 4, 3};
const Prototype FPARAM = Prototype{"FPARAM", 9, 2, 3};

const std::vector<Prototype> FLHA_BLOCKS = {FCINFO, FMODSEL, FMASS, FLIFE, FCONST, FCONSTRATIO, FBAG, FWCOEF, IMFWCOEF, FOBS, FOBSERR, FOBSSM, FDIPOLE, FPARAM};

/**
 * @class LhaBlock
 * @brief Represents an LHA block with multiple elements, based on a given prototype.
 */
class LhaBlock {
private:
    Prototype prototype; /**< The prototype defining the structure of the block. */
    std::vector<std::unique_ptr<AbstractElement>> entries; /**< List of elements contained in the block. */

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
     * @brief Retrieves an element by its identifier.
     * @param id Identifier of the element.
     * @return Pointer to the `AbstractElement` if found; otherwise, `nullptr`.
     */
    AbstractElement* get(const std::string& id) const;

    /**
     * @brief Retrieves all entries in the block.
     * @return Pointer to a vector of unique pointers to `AbstractElement`.
     */
    const std::vector<std::unique_ptr<AbstractElement>>* getEntries() const;

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
     * @brief Destructor for LhaBlock, clears all entries.
     */
    inline ~LhaBlock() {entries.clear();}
};

#endif // HYPERISO_LHA_BLOCKS_H