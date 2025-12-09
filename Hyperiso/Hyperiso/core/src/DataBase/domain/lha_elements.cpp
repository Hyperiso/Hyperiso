#include "lha_elements.h"

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
        } catch (const std::invalid_argument&) {
            throw std::runtime_error("Invalid double: '" + str + "'");
        } catch (const std::out_of_range&) {
            throw std::runtime_error("Out-of-range double: '" + str + "'");
        }
    }
};


template <>
struct StringConverter<std::string> {
    static std::string convert(const std::string& str) {
        return str;
    }
};

static inline void normalize_indices(const Prototype& p,
                                     const std::vector<std::string>& line,
                                     size_t& vIdx, int& sIdx, int& rIdx, int& bIdx)
{
    vIdx = p.valueIdx < line.size() ? p.valueIdx : (line.empty() ? 0 : line.size() - 1);

    sIdx = (p.scaleIdx >= 0 && static_cast<size_t>(p.scaleIdx) < line.size()) ? p.scaleIdx : -1;
    rIdx = (p.rgIdx    >= 0 && static_cast<size_t>(p.rgIdx)    < line.size()) ? p.rgIdx    : -1;
    bIdx = (p.binIdx   >= 0 && static_cast<size_t>(p.binIdx)   < line.size()) ? p.binIdx   : -1;
}

template<typename T>
LhaElement<T>::LhaElement(const Prototype& prototype, const std::vector<std::string>& line)
    : AbstractElement(encodeId(prototype, line))
{
    size_t vIdx; int sIdx, rIdx, bIdx;
    normalize_indices(prototype, line, vIdx, sIdx, rIdx, bIdx);

    if (prototype.globalScale) {
        if (line.empty()) throw std::runtime_error("Global-scale block: empty line in " + prototype.blockName);
        this->Q.emplace(std::stod(line.at(0))); // Q injecté en [0]
    } else if (sIdx != -1) {
        this->Q.emplace(std::stod(line.at(sIdx)));
    }
    if (rIdx != -1) {
        this->rScheme.emplace(static_cast<RenormalizationScheme>(std::stoi(line.at(rIdx))));
    }

    if (bIdx != -1) {
        this->bin.emplace(std::pair(std::stod(line.at(bIdx)), std::stod(line.at(bIdx + 1))));
    }

    if (vIdx >= line.size())
        throw std::runtime_error("valueIdx out of range in " + prototype.blockName);
    this->value = StringConverter<T>::convert(line.at(vIdx));
}

template <typename T>
LhaID LhaElement<T>::encodeId(const Prototype& prototype, const std::vector<std::string>& line) {
    size_t vIdx; int sIdx, rIdx, bIdx;
    normalize_indices(prototype, line, vIdx, sIdx, rIdx, bIdx);

    std::vector<long> sub_ids;
    for (size_t i = 0; i < line.size(); ++i) {
        if (i == vIdx || static_cast<int>(i) == sIdx || static_cast<int>(i) == rIdx) continue;
        if (bIdx != -1 && (static_cast<int>(i) == bIdx  || static_cast<int>(i) == bIdx + 1)) continue;
        if (prototype.globalScale && i == 0) continue;

        const auto& s = line[i];
        if (s.find_first_of(".eEdD") != std::string::npos) continue;

        try {
            // std::cout << "truc: " <<  s << std::endl;
            sub_ids.emplace_back(std::stol(s));
        } catch (...) {
           LOG_WARN("Non-integer ID token in ", prototype.blockName, ": '", s, "'");
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

template <typename T>
std::shared_ptr<Node> LhaElement<T>::toDBNode() const {
    Node node;
    node.set(this->getValue(), "central_value");
    if (Q.has_value()) {
        node.set(this->getScale(), "scale");
    }   
    if (rScheme.has_value()) {
        node.set(static_cast<int>(this->getScheme()), "renormalization_scheme");
    }
    if (bin.has_value()) {
        auto bin = this->getBinning();
        node.set(bin.first, "bin_low");
        node.set(bin.second, "bin_high");
    }  
    return std::make_shared<Node>(node);
}

std::shared_ptr<AbstractElement> LhaElementFactory::createElement(const Prototype& prototype, const std::vector<std::string>& line) {
    if (prototype.blockName == "FCINFO" || prototype.blockName == "FMODSEL" || prototype.blockName == "SPINFO") {
        return std::make_shared<LhaElement<std::string>>(prototype, line);
    } else {
        return std::make_shared<LhaElement<double>>(prototype, line);
    }
}
