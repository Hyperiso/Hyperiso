#include <string>
#include <sstream>
#include <memory>
#include <algorithm>

#include "lha_blocks.h"
#include "lha_elements.h"

AbstractElement* LhaBlock::get(const std::string& id) const
{
    auto p = [id](const std::unique_ptr<AbstractElement>& e) { 
        return e->getId() == id; 
    };
    return (*std::find_if(entries.begin(), entries.end(), p)).get();
}

std::string LhaBlock::toString() const {
    std::stringstream stream;
    stream << "Block " << this->name << ":\n";
    for (const auto& entry: entries) {
        stream << entry->toString();
    }
    return stream.str();
}

void LhaBlock::addElement(const std::vector<std::string> &line)
{
    auto elt = LhaElementFactory::createElement(this->name, line);
    this->entries.emplace_back(std::move(elt));
}

void LhaBlock::readData(const std::vector<std::vector<std::string>> &lines)
{
    for (auto line : lines) {
        if (!line.empty()) {
            this->addElement(line);
        }
    }
}
