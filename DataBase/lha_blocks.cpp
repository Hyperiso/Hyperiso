#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <memory>
#include <algorithm>

#include "lha_blocks.h"
#include "lha_elements.h"

std::vector<std::string> LhaBlock::parseLine(const std::string &line)
{
    std::istringstream stream(line);
    std::vector<std::string> words {};
    std::string word;
    while (stream >> word) {
        if (word.at(0) == '#') break;
        words.emplace_back(word);
    }
    return words;
}

AbstractElement* LhaBlock::get(std::string_view id) const
{
    auto p = [id](const std::unique_ptr<AbstractElement>& e) { 
        return e->getId() == id; 
    };
    return (*std::find_if(entries.begin(), entries.end(), p)).get();
}

std::string LhaBlock::toString() const {
    std::stringstream stream;
    stream << "Block " << BlockIdHelper::getBlockName(this->id) << ":\n";
    for (auto& entry: entries) {
        stream << entry->toString();
    }
    return stream.str();
}

void LhaBlock::readData(std::ifstream& file) {
    std::string line;
    while (std::getline(file, line)) {
        if (line.size() != 0 && line.at(0) != '#') {
            auto words = LhaBlock::parseLine(line);
            auto entry = createElement(words);
            this->entries.emplace_back(std::move(entry));
        }
        if (tolower(file.peek()) == 'b' || file.peek() == EOF) break;
    }
}

std::unique_ptr<AbstractElement> MassBlock::createElement(const std::vector<std::string>& words) {
    if (words.size() != 4) {
        std::cout << "Error: invalid mass entry\n";
        // Logger::error("Invalid particle mass entry");
    }
    double value = std::stod(words.at(1));
    RenormalizationScheme scheme = static_cast<RenormalizationScheme>(std::stoi(words.at(2)));
    double scale = std::stod(words.at(3));
    if (value < 0) {
        // Logger::error("Invalid mass value: negative");
    }
    if (scale < 0) {
        // Logger::warn("Invalid renormalization scale: negative. Set to 0.");
    }
    return std::make_unique<GeneralElement<double>>(words.at(0), value, scheme, scale);
}

std::unique_ptr<AbstractElement> InfoBlock::createElement(const std::vector<std::string>& words) {
    if (words.size() != 2) {
        // Logger::error("Invalid information entry");
    }
    return std::make_unique<TypedElement<std::string>>(words.at(0), words.at(1));
}

std::unique_ptr<LhaBlock> LhaBlockFactory::createBlock(BlockId id, bool isFLHA) {
    if (id == BlockId::MASS && !isFLHA || id == BlockId::FMASS && isFLHA) {
        return std::make_unique<MassBlock>(id);
    } else if (id == BlockId::FCINFO && isFLHA) {
        return std::make_unique<InfoBlock>(id);
    }
    // TODO : Parse all block names (switch ?)
    return nullptr;
}