#pragma once
#include <functional>
#include <memory>
#include <unordered_map>
#include "GroupDefinition.h"
#include "Wilson.h" // WilsonCoefficient

using CoefPtr   = std::shared_ptr<WilsonCoefficient>;
using CoefMaker = std::function<CoefPtr(const BuildContext&, WCoef)>;

class CoefficientRegistry {
public:
    void register_creator(WCoef c, Model m, Backend b, CoefMaker mk) {
        table_[key(c,m,b)] = std::move(mk);
    }
    CoefPtr create(const BuildContext& ctx, WCoef c) const {
        if (auto it = table_.find(key(c, ctx.model, ctx.backend)); it != table_.end())
            return it->second(ctx, c);

        // fallback Marty -> Builtin
        if (ctx.backend == Backend::Marty) {
            if (auto it2 = table_.find(key(c, ctx.model, Backend::Builtin)); it2 != table_.end())
                return it2->second(ctx, c);
        }
        // fallback modèle -> SM
        if (ctx.model != Model::SM) {
            if (auto it3 = table_.find(key(c, Model::SM, ctx.backend)); it3 != table_.end())
                return it3->second(ctx, c);
        }
        throw std::runtime_error("No factory registered for coefficient");
    }
private:

    struct Key {
        int c;
        int m;
        int b;
        bool operator==(const Key& o) const noexcept {
            return c==o.c && m==o.m && b==o.b;
        }
    };
    struct KeyHash {
        std::size_t operator()(const Key& k) const noexcept {
            // hash (boost-like)
            std::size_t h = std::hash<int>{}(k.c);
            h ^= std::hash<int>{}(k.m) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
            h ^= std::hash<int>{}(k.b) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
            return h;
        }
    };
    static Key key(WCoef c, Model m, Backend b) {
        return Key{(int)c, (int)m, (int)b};
    }
    std::unordered_map<Key, CoefMaker, KeyHash> table_;
};

void register_B(CoefficientRegistry&);
void register_BPrime(CoefficientRegistry&);
void register_BScalar(CoefficientRegistry&);
void register_CC_bc(CoefficientRegistry&);
void register_CC_bu(CoefficientRegistry&);
void register_CC_cs(CoefficientRegistry&);
void register_CC_cd(CoefficientRegistry&);
void register_CC_su(CoefficientRegistry&);
void register_MesonMixing(CoefficientRegistry&);
void register_K(CoefficientRegistry& reg);
