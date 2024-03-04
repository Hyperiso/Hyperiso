#include "lha_elements.h"
#include <sstream>

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
LhaElement<T>::LhaElement(const std::string& block, const std::vector<std::string>& line) : block(block), AbstractElement(encodeId(block, line)) {
    if (SCALE_BLOCKS.contains(block)) {
        this->Q.emplace(std::stod(line.at(SCALE_BLOCKS.at(block))));
        if (SCHEME_BLOCKS.contains(block)) {
            // WRONG !!! To be corrected.
            this->rScheme.emplace(static_cast<RenormalizationScheme>(stoi(line.at(SCHEME_BLOCKS.at(block)))));
        }
    }
    int value_index = VALUE_POS.contains(block) ? VALUE_POS.at(block) : 1;
    this->value = StringConverter<T>::convert(line.at(value_index));
}

template <typename T>
std::string LhaElement<T>::encodeId(const std::string& block, const std::vector<std::string>& line) {
    std::stringstream stream;

    if (block == "FCINFO" || block == "FMODSEL" || block == "SMINPUTS" || block == "VCKMIN" || block == "UPMNSIN" || block == "FMASS" || block == "FLIFE") {
        stream << line.at(0);
    } else if (block == "FCONST" || block == "FBAG") {
        stream << line.at(0) << "|" << line.at(1);
    } else if (block == "FCONSTRATIO") {
        stream << line.at(0) << "|" << line.at(1) << "|" << line.at(2) << "|" << line.at(3);
    } else if (block.ends_with("FWCOEF")) {
        stream << line.at(1) << "|" << line.at(2) << "|" << line.at(3) << "|" << line.at(4);
    } else if (block.starts_with("FOBS")) {
        stream << line.at(0) << "|" << line.at(1) << "|" << line.at(4);
        for (int i=5; i!=5 + std::stoi(line.at(4)); ++i) {
            stream << "|" << line.at(i);
        }
    } else if (block == "FDIPOLE") {
        stream << line.at(0) << "|" << line.at(1) << "|" << line.at(2);
    } else if (block == "FPARAM") {
        stream << line.at(0);
        Logger::getInstance()->warn("Treatment of block FPARAM should be user-defined. Please check results.");
    } else {
        Logger::getInstance()->error("Unknown block. Please define appropriate behaviour.");
    }

    return stream.str();
}

template <typename T>
std::string LhaElement<T>::toString() const {
    std::stringstream stream;
    stream << this->getId() << '\t' << this->getValue() << '\t' << static_cast<int>(this->getScheme()) << '\t' << this->getScale() << "\n";
    return stream.str();
}

std::unique_ptr<AbstractElement> LhaElementFactory::createElement(const std::string& blockName, const std::vector<std::string>& line) {
    if (blockName == "FCINFO" || blockName == "FMODSEL") {
        return std::make_unique<LhaElement<std::string>>(blockName, line);
    } else {
        return std::make_unique<LhaElement<double>>(blockName, line);
    }
}