#ifndef LHABLOCKPROTOTYPE_H
#define LHABLOCKPROTOTYPE_H

#include <string>
#include <unordered_set>
#include "BlockName.h"

/**
 * @struct Prototype
 * @brief Represents the structure of an LHA block, specifying columns for values, scales, and renormalization groups.
 */
struct Prototype {
    BlockName blockName;      /**< Block name, case insensitive. */
    size_t itemCount {2};       /**< Number of columns in the block. */
    size_t valueIdx {1};        /**< Column index for value. */
    int scaleIdx {-1};          /**< Column index for scale, -1 if scale-independent. */
    int rgIdx {-1};             /**< Column index for renormalization group, -1 if irrelevant. */
    int binIdx {-1};          /**< Column index for bin lower bound, -1 if unbinned. */
    bool globalScale {false};   /**< Indicates if the block uses a global scale (Q= in header). */

    bool operator==(const Prototype& other) const {
        return blockName == other.blockName &&
               itemCount == other.itemCount &&
               valueIdx == other.valueIdx &&
               scaleIdx == other.scaleIdx &&
               binIdx == other.binIdx &&
               rgIdx == other.rgIdx &&
               globalScale == other.globalScale;
    }
};

namespace std {
    template <>
    struct hash<Prototype> {
        std::size_t operator()(const Prototype& p) const noexcept {
            std::size_t h1 = std::hash<BlockName>{}(p.blockName);
            std::size_t h2 = std::hash<size_t>{}(p.itemCount);
            std::size_t h3 = std::hash<size_t>{}(p.valueIdx);
            std::size_t h4 = std::hash<int>{}(p.scaleIdx);
            std::size_t h5 = std::hash<int>{}(p.rgIdx);
            std::size_t h6 = std::hash<bool>{}(p.globalScale);

            return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3) ^ (h5 << 4) ^ (h6 << 5);
        }
    };
}

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
const Prototype NMAMIX = Prototype{"NMAMIX", 3, 2};
const Prototype NMHMIX = Prototype{"NMHMIX", 3, 2};
const Prototype STOPMIX = Prototype{"STOPMIX", 3, 2};
const Prototype SBOTMIX = Prototype{"SBOTMIX", 3, 2};
const Prototype STAUMIX = Prototype{"STAUMIX", 3, 2};
const Prototype ALPHA = Prototype{"ALPHA", 1, 0};
const Prototype HMIX = Prototype{"HMIX", 3, 2, -1, -1, -1, true};
const Prototype GAUGE = Prototype{"GAUGE", 3, 2, -1, -1, -1, true};
const Prototype MSOFT = Prototype{"MSOFT", 3, 2, -1, -1, -1, true};
const Prototype AU = Prototype{"AU", 4, 3, -1, -1, -1, true};
const Prototype AD = Prototype{"AD", 4, 3, -1, -1, -1, true};
const Prototype AE = Prototype{"AE", 4, 3, -1, -1, -1, true};
const Prototype YU = Prototype{{"YU", "UCOUPL"}, 4, 3, -1, -1, -1, true};
const Prototype YD = Prototype{{"YD", "DCOUPL"}, 4, 3, -1, -1, -1, true};
const Prototype YE = Prototype{{"YE", "LCOUPL", "YL"}, 4, 3, -1, -1, -1, true};

const std::unordered_set<Prototype> LHA_BLOCKS = {MODSEL, SMINPUTS, VCKMIN, UPMNSIN, MASS, GAUGE};
const std::unordered_set<Prototype> SLHA_BLOCKS = {MINPAR, EXTPAR, NMIX, UMIX, VMIX, STOPMIX, SBOTMIX, STAUMIX, ALPHA, HMIX, MSOFT, AU, AD, AE, YU, YD, YE, NMHMIX, NMAMIX}; // MAJ : SPINFO/FCINFO pb (values are string)

// FLHA Block prototypes
const Prototype FCINFO = Prototype{"FCINFO"};
const Prototype FMODSEL = Prototype{"FMODSEL"};
const Prototype FMASS = Prototype{"FMASS", 4, 1, 2, 3};
const Prototype FLIFE = Prototype{"FLIFE"};
const Prototype FCONST = Prototype{"FCONST", 5, 2, 3, 4};
const Prototype FCONSTRATIO = Prototype{"FCONSTRATIO", 7, 4, 5, 6};
const Prototype FBAG = Prototype{"FBAG", 5, 2, 3, 4};
// const Prototype FWCOEF = Prototype{"FWCOEF", 6, 5, -1, -1, true};
// const Prototype IMFWCOEF = Prototype{"IMFWCOEF", 6, 5, -1, -1, true};
const Prototype FWCOEF   = Prototype{"FWCOEF",   7, 6, -1, -1, -1, true};
const Prototype IMFWCOEF = Prototype{"IMFWCOEF", 7, 6, -1, -1, -1, true};
const Prototype FOBS = Prototype{"FOBS", 11, 2, 3, -1, 4};
const Prototype FOBSERR = Prototype{"FOBSERR", 11, 2, 3, -1, 4};
const Prototype FOBSSM = Prototype{"FOBSSM", 11, 2, 3, -1, 4};
const Prototype FDIPOLE = Prototype{"FDIPOLE", 4, 3};
const Prototype FPARAM = Prototype{"FPARAM", 9, 2, 3};

const std::unordered_set<Prototype> FLHA_BLOCKS = {FCINFO, FMODSEL, FMASS, FLIFE, FCONST, FCONSTRATIO, FBAG, FWCOEF, IMFWCOEF, FOBS, FOBSERR, FOBSSM, FDIPOLE, FPARAM};


#endif // LHABLOCKPROTOTYPE_H
