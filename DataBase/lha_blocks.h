#ifndef HYPERISO_LHA_BLOCKS_H
#define HYPERISO_LHA_BLOCKS_H

#include <string>
#include <memory>
#include <vector>

#include "lha_elements.h"

class LhaBlock {
private:
    std::vector<std::unique_ptr<AbstractElement>> entries; 
    std::string name;

public:
    inline explicit LhaBlock(const std::string& name) : name(name) {}
    void readData(const std::vector<std::vector<std::string>>& lines);
    AbstractElement* get(const std::string& id) const;
    void addElement(const std::vector<std::string>& line);
    std::string toString() const;
};

#endif // HYPERISO_LHA_BLOCKS_H