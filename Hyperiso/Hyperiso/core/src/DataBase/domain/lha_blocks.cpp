#include "lha_blocks.h"

AbstractElement* LhaBlock::get(const LhaID& id) const {
    LOG_DEBUG("Trying to retrieve element", id, "within block", this->prototype.blockName);
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

std::shared_ptr<Node> LhaBlock::toDBNode() const {
    Node node;
    std::map<BlockName, Node::Value> elts_as_nodes;
    for (const auto& elt : this->entries) {
        elts_as_nodes.emplace(elt->getId().to_string(), elt->toDBNode());
    }
    if (!this->entries.empty() && this->prototype.globalScale) {
        double scale = entries[0]->getScale();
        node.set(scale, "scale");
    }
    node.setGroup({this->prototype.blockName}, elts_as_nodes);
    return std::make_shared<Node>(node);
}

void LhaBlock::addElement(const std::vector<std::string> &line) {
    auto elt = LhaElementFactory::createElement(this->prototype, line);
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
