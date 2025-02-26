#include "lha_elements.h"
#include "lha_blocks.h"
#include <sstream>
#include <iostream>
#include <limits>

template <typename U>
struct StringConverter {
    static U convert(const std::string& str) {
        static_assert(std::is_same_v<U, void>, "Unsupported conversion");
    }
};

template <>
struct StringConverter<double> {
    static double convert(const std::string& str) {
        try {
            return std::stod(str);
        } catch(...) {
            return std::numeric_limits<double>::max(); //small patch
        }
    }
};

template <>
struct StringConverter<std::string> {
    static std::string convert(const std::string& str) {
        return str;
    }
};

template<typename T>
LhaElement<T>::LhaElement(const Prototype& prototype, const std::vector<std::string>& line) 
        : AbstractElement(encodeId(prototype, line)) {
    if (prototype.scaleIdx != -1) {
        this->Q.emplace(std::stod(line.at(prototype.scaleIdx)));
        if (prototype.rgIdx != -1) {
            // WRONG !!! To be corrected.
            this->rScheme.emplace(static_cast<RenormalizationScheme>(stoi(line.at(prototype.rgIdx))));
        }
    } else if (prototype.globalScale) {
        this->Q.emplace(std::stod(line.at(0)));
    }

    this->value = StringConverter<T>::convert(line.at(prototype.valueIdx));
}

template <typename T>
LhaID LhaElement<T>::encodeId(const Prototype& prototype, const std::vector<std::string>& line) {
    std::vector<long> sub_ids;
    for (size_t i=0; i!=line.size(); ++i) {
        if (i != prototype.valueIdx && i != prototype.scaleIdx && i != prototype.rgIdx) {
            if (prototype.globalScale && i == 0) continue;
            sub_ids.emplace_back(std::stol(line.at(i)));
        }
    }
    return LhaID(sub_ids);
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

std::shared_ptr<AbstractElement> LhaElementFactory::createElement(LhaBlock* block, const std::vector<std::string>& line) {
    if (block->getPrototype().blockName == "FCINFO" || block->getPrototype().blockName == "FMODSEL" || block->getPrototype().blockName == "SPINFO") {
        return std::make_shared<LhaElement<std::string>>(block->getPrototype(), line);
    } else {
        return std::make_shared<LhaElement<double>>(block->getPrototype(), line);
    }
}
