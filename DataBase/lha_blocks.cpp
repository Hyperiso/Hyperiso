#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <memory>
#include <algorithm>
#include <span>
#include <complex>
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
            if (!words.empty()) {
                auto entry = createElement(words);
                this->entries.emplace_back(std::move(entry));
            }  
        }
        if (tolower(file.peek()) == 'b' || file.peek() == EOF) break;
    }
}

int WilsonBlock::findImOffset(std::ifstream& file) {
    int currentPos = file.tellg();
    int offset = 0;
    file.seekg(0, std::ios_base::beg);
    std::string line;
    
    while (std::getline(file, line)) {
        if (tolower(line[0]) != 'b') continue;
        std::string _, name;
        std::istringstream stream(line);
        stream >> _ >> name;
        if (name == "IMWCOEFF") {
            offset = file.tellg() - currentPos;
        }
    }
    file.seekg(currentPos, std::ios_base::beg);
    return offset;
}

// TODO: read scale Q
void WilsonBlock::readData(std::ifstream &file) {
    std::string line;
    int imOffset = this->findImOffset(file);
    while (std::getline(file, line)) {
        if (line.size() != 0 && line.at(0) != '#') {
            auto words = LhaBlock::parseLine(line);
            file.seekg(imOffset, std::ios_base::cur);
            std::getline(file, line);
            auto wordsIm = LhaBlock::parseLine(line);
            file.seekg(-imOffset, std::ios_base::cur);
            words.emplace_back(wordsIm.at(4));
            if (!words.empty()) {
                auto entry = createElement(words);
                this->entries.emplace_back(std::move(entry));
            }  
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
    return std::make_unique<GeneralElement<double>>(words.at(0), value, scheme, scale);
}

std::unique_ptr<AbstractElement> InfoBlock::createElement(const std::vector<std::string>& words) {
    if (words.size() != 2) {
        // Logger::error("Invalid information entry");
        std::cout << "Error: invalid info entry\n";
    }
    return std::make_unique<TypedElement<std::string>>(words.at(0), words.at(1));
}

std::unique_ptr<AbstractElement> SingleValueBlock::createElement(const std::vector<std::string>& words) {
    if (words.size() != 2) {
        // Logger::error("Invalid information entry");
        std::cout << "Error: invalid value entry\n";
    }
    double value = std::stod(words.at(1));
    return std::make_unique<TypedElement<double>>(words.at(0), value);
}

std::unique_ptr<AbstractElement> DecayConstBlock::createElement(const std::vector<std::string>& words) {
    if (words.size() != 5) {
        // Logger::error("Invalid information entry");
        std::cout << "Error: invalid decay constant entry\n";
    }
    double value = std::stod(words.at(2));
    RenormalizationScheme scheme = static_cast<RenormalizationScheme>(std::stoi(words.at(3)) + 5);
    double scale = std::stod(words.at(4));
    std::ostringstream id;
    id << std::setw(9) << std::setfill('0') << words.at(0) << std::setw(2) << words.at(1);
    return std::make_unique<GeneralElement<double>>(id.str(), value, scheme, scale);
}

std::unique_ptr<AbstractElement> DecayConstRatioBlock::createElement(const std::vector<std::string>& words) {
    if (words.size() != 7) {
        // Logger::error("Invalid information entry");
        std::cout << "Error: invalid decay constant ratio entry\n";
    }
    double value = std::stod(words.at(4));
    RenormalizationScheme scheme = static_cast<RenormalizationScheme>(std::stoi(words.at(5)) + 5);
    double scale = std::stod(words.at(6));
    std::ostringstream id;
    id << std::setfill('0') 
        << std::setw(9) << words.at(0) 
        << std::setw(9) << words.at(1) 
        << std::setw(2) << words.at(2) 
        << std::setw(2) << words.at(3);
    return std::make_unique<GeneralElement<double>>(id.str(), value, scheme, scale);
}

std::unique_ptr<AbstractElement> BagBlock::createElement(const std::vector<std::string> &words)
{
    if (words.size() != 5) {
        // Logger::error("Invalid information entry");
        std::cout << "Error: invalid bag parameter entry\n";
    }
    double value = std::stod(words.at(2));
    RenormalizationScheme scheme = static_cast<RenormalizationScheme>(std::stoi(words.at(3)) + 5);
    double scale = std::stod(words.at(4));
    std::ostringstream id;
    id << std::setw(9) << std::setfill('0') << words.at(0) << std::setw(2) << words.at(1);
    return std::make_unique<GeneralElement<double>>(id.str(), value, scheme, scale);
}

std::unique_ptr<AbstractElement> ObsBlock::createElement(const std::vector<std::string> &words)
{
    if (words.size() < 6 || words.size() != (5 + std::stoi(words.at(4)))) {
        // Logger::error("Invalid information entry");
        std::cout << "Error: invalid observable entry\n";
    }
    double value = std::stod(words.at(2));
    double scale = std::stod(words.at(3));
    std::ostringstream id;
    id << std::setfill('0') << std::setw(9) << words.at(0) << std::setw(2) << words.at(1) << words.at(4); 
    for (const auto& w : std::span(words.begin() + 5, words.end())) {
        std::cout << w << std::endl;
        id << std::setw(9) << w;
    }
    return std::make_unique<ScaleDependentElement<double>>(id.str(), value, scale);
}

std::unique_ptr<AbstractElement> WilsonBlock::createElement(const std::vector<std::string> &words)
{
    if (words.size() != 7) {
        // Logger::error("Invalid information entry");
        std::cout << "Error: invalid Wilson coefficient entry\n";
    }
    complex_t value = complex_t(std::stod(words.at(4)), std::stod(words.at(5)));
    double scale = std::stod(words.at(6));
    std::ostringstream id;
    id << std::setfill('0') 
        << std::setw(8) << words.at(0) 
        << std::setw(4) << words.at(1)
        << std::setw(2) << words.at(2) 
        << words.at(3); 
    return std::make_unique<ScaleDependentElement<complex_t>>(id.str(), value, scale);
}


bool checkFLHA(BlockId id, bool isFLHA) {
    if (!isFLHA) {
        // Logger::error(BlockIdHelper::getBlockName(id) + " block not allowed in SLHA file");
    }
    return isFLHA;
}

std::unique_ptr<LhaBlock> LhaBlockFactory::createBlock(BlockId id, bool isFLHA) {
    switch (id) {
        case BlockId::MASS:
            return std::make_unique<MassBlock>(id);
        case BlockId::FMASS:
            return checkFLHA(id, isFLHA) ? std::make_unique<MassBlock>(id) : nullptr;
        case BlockId::FCINFO:
            return checkFLHA(id, isFLHA) ? std::make_unique<InfoBlock>(id) : nullptr;
        case BlockId::SMINPUTS:
            return std::make_unique<SMInputsBlock>(id);
        case BlockId::VCKMIN:
            return std::make_unique<CKMParamBlock>(id);
        case BlockId::UPMNSIN:
            return std::make_unique<PMNSParamBlock>(id);
        case BlockId::FLIFE:
            return checkFLHA(id, isFLHA) ? std::make_unique<LifetimeBlock>(id) : nullptr;
        case BlockId::FCONST:
            return checkFLHA(id, isFLHA) ? std::make_unique<DecayConstBlock>(id) : nullptr;
        case BlockId::FCONSTRATIO:
            return checkFLHA(id, isFLHA) ? std::make_unique<DecayConstRatioBlock>(id) : nullptr;
        case BlockId::FBAG:
            return checkFLHA(id, isFLHA) ? std::make_unique<BagBlock>(id) : nullptr;
        case BlockId::FOBS:
            return checkFLHA(id, isFLHA) ? std::make_unique<ObsBlock>(id) : nullptr;
        case BlockId::FOBSERR:
            return checkFLHA(id, isFLHA) ? std::make_unique<ObsBlock>(id) : nullptr;
        case BlockId::FOBSSM:
            return checkFLHA(id, isFLHA) ? std::make_unique<ObsBlock>(id) : nullptr;
        default:
            return nullptr;
    }
}
