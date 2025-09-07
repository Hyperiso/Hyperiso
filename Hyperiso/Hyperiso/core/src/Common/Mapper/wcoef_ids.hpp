#pragma once
#include "generic_mapper.hpp"
#include "Map.h"
#include "hash_utils.hpp"
#include "GeneralEnum.h"
#include "General.h"

struct WCoefTag {};
using WCoefId = SymbolId<WCoefTag>;

class WCoefMapper
: public GenericMapperWithExt<
      WCoefTag, WCoef,
      std::pair<int,int>, PairHash,
      wcoef_mapping, wcoef_flha_mapping
  >
{
public:
    using Base = GenericMapperWithExt<WCoefTag,WCoef,std::pair<int,int>,PairHash,wcoef_mapping,wcoef_flha_mapping>;
    using Base::str; using Base::enum_elt; using Base::from_external;
    using Base::external_of; using Base::to_id; using Base::enum_of;

    static WCoef from_flha(int a,int b){
        static const auto inv = []{
            std::map<std::pair<int,int>, WCoef> m;
            for (auto& [e,p] : wcoef_flha_mapping()) m.emplace(p,e);
            return m;
        }();
        auto it = inv.find({a,b});
        if (it==inv.end()) throw std::out_of_range("Unsupported FLHA key");
        return it->second;
    }
    static std::optional<WCoefId> from_flha_key(int a,int b){ return Base::from_external({a,b}); }

    static LhaID flha_full(WCoef e, QCDOrder q, ContributionType c){
        auto [x,y] = wcoef_flha_mapping().at(e);
        return LhaID{x,y, int(q)-1, int(c)};
    }
    static LhaID flha_full(const WCoefId& id, QCDOrder q, ContributionType c){
        auto k = external_of(id);
        if (!k) throw std::runtime_error("No FLHA key for "+id.str());
        return LhaID{k->first, k->second, int(q)-1, int(c)};
    }

    static size_t n_wilsons(){ return wcoef_mapping().size(); }

        static const std::map<std::pair<int,int>, WCoef>& inverse_flha_mapping() {
        static const std::map<std::pair<int,int>, WCoef> inv = []{
            std::map<std::pair<int,int>, WCoef> m;
            for (const auto& [e, ab] : wcoef_flha_mapping()) m.emplace(ab, e);
            return m;
        }();
        return inv;
    }

    static const std::vector<WCoef>& B_group(){
        static const std::vector<WCoef> g = {
            WCoef::C1,WCoef::C2,WCoef::C3,WCoef::C4,WCoef::C5,WCoef::C6,
            WCoef::C7,WCoef::C8,WCoef::C9,WCoef::C10
        }; return g;
    }
    static const std::vector<WCoef>& B_prime_group(){
        static const std::vector<WCoef> g = {
            WCoef::CP1,WCoef::CP2,WCoef::CP3,WCoef::CP4,WCoef::CP5,WCoef::CP6,
            WCoef::CP7,WCoef::CP8,WCoef::CP9,WCoef::CP10,WCoef::CPQ1,WCoef::CPQ2
        }; return g;
    }
    static const std::vector<WCoef>& B_scalar_group(){
        static const std::vector<WCoef> g = { WCoef::CQ1, WCoef::CQ2 }; return g;
    }
    static const std::vector<WCoef>& b_clnu_group(){
        static const std::vector<WCoef> g = { WCoef::C_V1, WCoef::C_V2, WCoef::C_S1, WCoef::C_S2, WCoef::C_T }; return g;
    }
    static const std::vector<WCoef>& meson_mixing_group(){
        static const std::vector<WCoef> g = {
            WCoef::C_BD_1, WCoef::C_BD_2, WCoef::C_BD_3, WCoef::C_BD_4, WCoef::C_BD_5,
            WCoef::CT_BD_1,WCoef::CT_BD_2,WCoef::CT_BD_3,
            WCoef::C_BS_1, WCoef::C_BS_2, WCoef::C_BS_3, WCoef::C_BS_4, WCoef::C_BS_5,
            WCoef::CT_BS_1,WCoef::CT_BS_2,WCoef::CT_BS_3,
            WCoef::C_SD_1, WCoef::C_SD_2, WCoef::C_SD_3, WCoef::C_SD_4, WCoef::C_SD_5,
            WCoef::CT_SD_1,WCoef::CT_SD_2,WCoef::CT_SD_3,
            WCoef::C_CU_1, WCoef::C_CU_2, WCoef::C_CU_3, WCoef::C_CU_4, WCoef::C_CU_5,
            WCoef::CT_CU_1,WCoef::CT_CU_2,WCoef::CT_CU_3
        }; return g;
    }

    static std::vector<WCoef> get_group(WGroup g){
        switch (g){
            case WGroup::B: return {B_group().begin(), B_group().end()};
            case WGroup::BPrime: return {B_prime_group().begin(), B_prime_group().end()};
            case WGroup::BScalar: return {B_scalar_group().begin(), B_scalar_group().end()};
            case WGroup::BCC: return {b_clnu_group().begin(), b_clnu_group().end()};
            case WGroup::MESON_MIXING: return {meson_mixing_group().begin(), meson_mixing_group().end()};
            default: LOG_ERROR("Invalid WGroup","get_group couldn't find your group"); return {};
        }
    }
};
