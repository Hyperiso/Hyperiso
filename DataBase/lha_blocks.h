#if !defined(HYPERISO_LHA_BLOCKS_H)
#define HYPERISO_LHA_BLOCKS_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <algorithm>

#include "lha_elements.h"

// TODO : Add SLHA block names 
enum class BlockId {  FCINFO, FMODSEL, SMINPUTS, VCKMIN, UPMNSIN, VCKM, IMVCKM,
                        UPMNS, IMUPMNS, FMASS, MASS, FLIFE, FCONST, FCONSTRATIO, FBAG,
                        FWCOEF, IMFWCOEF, FOBS, FOBSERR, FOBSSM, FDIPOLE, FPARAM};

const std::vector<std::string> blockNames = {
    "FCINFO", "FMODSEL", "SMINPUTS", "VCKMIN", "UPMNSIN",
    "VCKM", "IMVCKM", "UPMNS", "IMUPMNS", "FMASS", "MASS",
    "FLIFE", "FCONST", "FCONSTRATIO", "FBAG", "FWCOEF",
    "IMFWCOEF", "FOBS", "FOBSERR", "FOBSSM", "FDIPOLE", "FPARAM"
};

class BlockIdHelper {
public:
    static std::string getBlockName(BlockId value) {
        return blockNames[static_cast<int>(value)];
    }

    static BlockId getBlockId(const std::string& value) {
        int idx = std::find(blockNames.begin(), blockNames.end(), value) - blockNames.begin();
        return static_cast<BlockId>(idx);
    }
};

class LhaBlock {
protected:
    std::vector<std::unique_ptr<AbstractElement>> entries;
    BlockId id;

public:
    LhaBlock(BlockId id) : id(id) {}
    void readData(std::ifstream& file);
    static std::vector<std::string> parseLine(const std::string& line);
    virtual std::unique_ptr<AbstractElement> createElement(const std::vector<std::string>& words) = 0;
    AbstractElement* get(std::string_view id) const;
    std::string toString() const;
};

class MassBlock : public LhaBlock {
public:
    MassBlock(BlockId id) : LhaBlock(id) {}
    std::unique_ptr<AbstractElement> createElement(const std::vector<std::string>& words) override;
};

class InfoBlock : public LhaBlock {
public:
    InfoBlock(BlockId id) : LhaBlock(id) {}
    std::unique_ptr<AbstractElement> createElement(const std::vector<std::string>& words) override;
}; 

class LhaBlockFactory {
public:
    static std::unique_ptr<LhaBlock> createBlock(BlockId id, bool isFLHA);
};

#endif // HYPERISO_LHA_BLOCKS_H