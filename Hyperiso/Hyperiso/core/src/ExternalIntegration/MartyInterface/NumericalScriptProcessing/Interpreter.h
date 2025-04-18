#ifndef INTERPRETER_H
#define INTERPRETER_H

#include "Extractor.h"
#include <string>
#include <unordered_map>
#include "MappingDataBase.h"
#include "ModelAPI.h"

class Interpreter {
public:
    struct InterpretedParam {
        std::string block;
        LhaID code;
        bool is_bsm;
        bool is_complex;

        bool operator==(const InterpretedParam& other) const {
            return block == other.block &&
                   code == other.code &&
                   is_bsm == other.is_bsm &&
                   is_complex == other.is_complex;
        }
    };

    Interpreter(const std::string& model = "SM");
    std::unordered_map<std::string, InterpretedParam> interpret(std::vector<Extractor::Parameter>& params);

private:
    std::shared_ptr<MappingDatabase> modelDatabase;
    std::shared_ptr<MappingDatabase> defaultDatabase;
};

namespace std {
    template <>
    struct hash<Interpreter::InterpretedParam> {
        std::size_t operator()(const Interpreter::InterpretedParam& p) const {
            std::size_t h1 = std::hash<std::string>{}(p.block);
            std::size_t h2 = std::hash<LhaID>{}(p.code);
            std::size_t h3 = std::hash<bool>{}(p.is_bsm);

            std::size_t seed = h1;
            seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);

            return seed;
        }
    };
}

#endif // INTERPRETER_H
