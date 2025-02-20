#include "lha_blocks.h"

AbstractElement* LhaBlock::get(const LhaID& id) const {
    auto p = [id](const std::shared_ptr<AbstractElement>& e) { 
        return e->getId() == id; 
    };
    auto element = std::find_if(entries.begin(), entries.end(), p);
    return element != entries.end() ? (*element).get() : nullptr;
}

const std::vector<std::shared_ptr<AbstractElement>>* LhaBlock::getEntries() const {
    return &(this->entries);
}

std::string LhaBlock::toString() const {
    std::stringstream stream;
    stream << "Block " << this->prototype.blockName << ":\n";
    for (const auto& entry: entries) {
        stream << entry->toString();
    }

    return stream.str();
}

void LhaBlock::addElement(const std::vector<std::string> &line) {
    auto elt = LhaElementFactory::createElement(this, line);
    this->entries.emplace_back(std::move(elt));
}

void LhaBlock::readData(const std::vector<std::vector<std::string>> &lines) {
    for (auto line : lines) {
        if (!line.empty()) { 
            this->addElement(line);
        }
    }
}

bool LhaBlock::hasElement(const LhaID& id) const {
    for (const auto& e : this->entries) {
        if (e->getId() == id) return true;
    }
    return false;
}
