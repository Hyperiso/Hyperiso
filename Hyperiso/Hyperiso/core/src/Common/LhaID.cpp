#include "LhaID.h"

std::ostream &operator<<(std::ostream &os, const LhaID &id) {
    os << id.to_string();
    return os;
}

LhaID::LhaID(const std::string &str_id) {
    std::cout << "herhe : " << str_id << " : " << str_id.size() << std::endl;
    if (str_id.size() == 0) {
        parts = {};
    } else {
        for (const auto &num : split(str_id, '_')) {
            parts.emplace_back(std::stol(num));
        }
    }
}

std::string LhaID::to_string() const {
    std::stringstream ss;
    if (!this->parts.empty()) {
        ss << this->parts.at(0);
        for (size_t i = 1; i < this->parts.size(); i++) {
            ss << '_' << this->parts.at(i);
        }
    }
    return ss.str();
}

LhaID::operator long() const {
    if (this->parts.size() > 1) {
        LOG_WARN("Casting nontrivial LhaID to int discards information.");
        for (auto& part : this->parts){
            std::cout << part<< std::endl;
        }
    }
    
    return this->parts.at(0);
};