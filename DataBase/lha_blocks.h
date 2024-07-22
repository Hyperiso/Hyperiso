#ifndef HYPERISO_LHA_BLOCKS_H
#define HYPERISO_LHA_BLOCKS_H

#include <string>
#include <memory>
#include <vector>

#include "lha_elements.h"

struct Prototype {
    std::string blockName;      // Case insensitive
    int itemCount {2};          // Number of columns
    int valueIdx {1};           // Value column 
    int scaleIdx {-1};          // Scale column, -1 if scale independent
    int rgIdx {-1};             // Renormalization group column, -1 if irrelevant
    bool globalScale {false};   // Blocks with Q= in header
};

// SLHA
const Prototype MODSEL = Prototype{"MODSEL"};
const Prototype SMINPUTS = Prototype{"SMINPUTS"};
const Prototype VCKMIN = Prototype{"VCKMIN"};   // SLHA2
const Prototype UPMNSIN = Prototype{"UPMNSIN"}; // SLHA2
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

// FLHA Exclusive
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


class LhaBlock {
private:
    Prototype prototype;
    std::vector<std::unique_ptr<AbstractElement>> entries;

public:
    inline explicit LhaBlock(const Prototype& prototype) : prototype(prototype) {}
    void readData(const std::vector<std::vector<std::string>>& lines);
    AbstractElement* get(const std::string& id) const;
    const std::vector<std::unique_ptr<AbstractElement>>* getEntries() const;
    void addElement(const std::vector<std::string>& line);
    inline Prototype getPrototype() { return this->prototype; };
    std::string toString() const;

    inline ~LhaBlock() {entries.clear();}
};

#endif // HYPERISO_LHA_BLOCKS_H