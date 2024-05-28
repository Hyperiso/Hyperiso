#include "lha_elements.h"
#include "lha_blocks.h"
#include <sstream>
#include <iostream>

template <typename U>
struct StringConverter {
    static U convert(const std::string& str) {
        static_assert(std::is_same_v<U, void>, "Unsupported conversion");
    }
};

template <>
struct StringConverter<double> {
    static double convert(const std::string& str) {
        return std::stod(str);
    }
};

template <>
struct StringConverter<std::string> {
    static std::string convert(const std::string& str) {
        return str;
    }
};

template<typename T>
LhaElement<T>::LhaElement(LhaBlock* block, const std::vector<std::string>& line) 
        : AbstractElement(block, encodeId(block, line)) {
    if (block->getPrototype().scaleIdx != -1) {
        this->Q.emplace(std::stod(line.at(block->getPrototype().scaleIdx)));
        if (block->getPrototype().rgIdx != -1) {
            // WRONG !!! To be corrected.
            this->rScheme.emplace(static_cast<RenormalizationScheme>(stoi(line.at(block->getPrototype().rgIdx))));
        }
    } else if (block->getPrototype().globalScale) {
        this->Q.emplace(std::stod(line.at(0)));
    }
    this->value = StringConverter<T>::convert(line.at(block->getPrototype().valueIdx));
}

template <typename T>
std::string LhaElement<T>::encodeId(LhaBlock* block, const std::vector<std::string>& line) {
    std::stringstream stream;

    Prototype p = block->getPrototype();
    for (int i=0; i!=line.size(); ++i) {
        if (i != p.valueIdx && i != p.scaleIdx && i != p.rgIdx) {
            if (p.globalScale && i == 0) continue;
            stream << line.at(i) << "|";
        }
    }

    auto id = stream.str();
    if (id.length() > 0) {
        id.erase(id.length() - 1);
    }
    return id;
}

template <typename T>
std::string LhaElement<T>::toString() const {
    std::stringstream stream;
    stream << this->getId() << '\t' << this->getValue();
    if (Q.has_value()) {
        stream << '\t' << this->getScale();
    }   
    if (rScheme.has_value()) {
        stream << '\t' << static_cast<int>(this->getScheme());
    }
    stream << "\n";
    return stream.str();
}

std::unique_ptr<AbstractElement> LhaElementFactory::createElement(LhaBlock* block, const std::vector<std::string>& line) {
    if (block->getPrototype().blockName == "FCINFO" || block->getPrototype().blockName == "FMODSEL" || block->getPrototype().blockName == "SPINFO") {
        return std::make_unique<LhaElement<std::string>>(block, line);
    } else {
        Logger::getInstance()->info(std::to_string(line.size())); 
        return std::make_unique<LhaElement<double>>(block, line);
    }
}