#ifndef WCOEF_IDS_H
#define WCOEF_IDS_H

#include "generic_mapper.h"
#include "Map.h"
#include "hash_utils.hpp"
#include "GeneralEnum.h"
#include "LhaID.h"

/**
 * @file wcoef_ids.h
 * @brief Mapping and helpers for Wilson coefficients.
 *
 * This header defines:
 *   - WCoefTag / WCoefId: strongly typed identifiers for Wilson coefficients,
 *   - WCoefMapper: mapper between
 *       * WCoef enum,
 *       * string identifiers,
 *       * FLHA base indices (std::pair<int,int>),
 *     plus a set of helpers for:
 *       * inverse FLHA lookup,
 *       * building full FLHA LhaID keys including QCD order and contribution type,
 *       * grouping coefficients by WGroup.
 */

/** @brief Tag type for Wilson-coefficient identifiers. */
struct WCoefTag {};

/** @brief Strongly typed identifier for Wilson coefficients (string-based). */
using WCoefId = SymbolId<WCoefTag>;

/**
 * @class WCoefMapper
 * @brief High-level mapper for WCoef <-> text <-> FLHA base indices.
 *
 * Specialization of GenericMapperWithExt with:
 *   - Tag         = WCoefTag
 *   - EnumT       = WCoef
 *   - ExternalKey = std::pair<int,int> (base FLHA indices)
 *   - Hash        = PairHash
 *   - MapFn       = wcoef_mapping
 *   - ExtFn       = wcoef_flha_mapping
 *
 * In addition to the generic mapping interface, it provides:
 *   - utilities to convert to/from FLHA base keys (int,int),
 *   - utilities to build full FLHA LhaID (including QCD order & contribution),
 *   - convenient accessors for specific Wilson groups (B, BPrime, CC_bc, ...),
 *   - an inverse FLHA mapping cache.
 */
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

    /**
     * @brief Returns the WCoef associated with a given FLHA base key (a,b).
     *
     * Uses a static inverse map built once from wcoef_flha_mapping().
     *
     * @param a First FLHA index.
     * @param b Second FLHA index.
     * @return Corresponding WCoef enum.
     *
     * @throws std::out_of_range if the pair (a,b) is not supported.
     */
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

    /**
     * @brief Registry-based lookup of WCoefId from a FLHA base key (a,b).
     *
     * Uses the external-key mapping in the GenericMapperWithExt registry.
     *
     * @param a First FLHA index.
     * @param b Second FLHA index.
     * @return Optional WCoefId if registered.
     */
    static std::optional<WCoefId> from_flha_key(int a,int b){ return Base::from_external({a,b}); }

    /**
     * @brief Returns the base FLHA pair (a,b) for a given WCoef.
     *
     * This uses the static mapping wcoef_flha_mapping().
     *
     * @param e WCoef enum value.
     * @return Pair of FLHA base indices.
     */
    static std::pair<int,int> flha_base(WCoef e){
        return wcoef_flha_mapping().at(e);
    }

    /**
     * @brief Returns the base FLHA pair (a,b) for a given WCoefId.
     *
     * This uses the registry's external_of() mapping. If no external key
     * is known for this id, an exception is thrown.
     *
     * @param id Wilson coefficient identifier.
     * @return Base FLHA indices (a,b).
     *
     * @throws std::runtime_error if @p id has no associated FLHA key.
     */
    static std::pair<int,int> flha_base(const WCoefId& id){
        auto k = external_of(id);
        if (!k) throw std::runtime_error("No FLHA key for " + id.str());
        return *k; // std::pair<int,int>
    }

    /**
     * @brief Builds a full FLHA LhaID for a given WCoef, order and contribution.
     *
     * The resulting LhaID has components:
     *   { base_a, base_b, qcd_order_index, contribution_type }
     * where:
     *   - (base_a, base_b) is from wcoef_flha_mapping(),
     *   - qcd_order_index = int(q) - 1,
     *   - contribution_type = int(c).
     *
     * @param e WCoef enum.
     * @param q QCD order.
     * @param c Contribution type (SM / BSM / TOTAL).
     * @return Corresponding LhaID.
     */
    static LhaID flha_full(WCoef e, QCDOrder q, ContributionType c){
        auto [x,y] = wcoef_flha_mapping().at(e);
        return LhaID{x,y, int(q)-1, int(c)};
    }

    /**
     * @brief Deserialize the id of a Wilson coefficient into (WCoefId, QCDOrder, ContributionType) using lha convention.
     *
     * The resulting object has components:
     *   std::pair<WCoefId, std::pair<QCDOrder, ContributionType>>
     * where:
     *   - WCoefId is the base id of the Wilson Coefficient (int,int)
     *   - QCDOrder the order of the coefficient (LO, NLO, NNLO)
     *   - ContributionType = SM, BSM or TOT (SM+BSM).
     *
     * @param id LhaID of the coefficient.
     */
    static std::pair<WCoefId, std::pair<QCDOrder, ContributionType>> lha_wilson_deserialize(LhaID id) {
        auto parts = id.get_parts();
        auto w_id = std::make_pair(parts[0], parts[1]);

        auto maybe = WCoefMapper::from_flha_key(w_id.first, w_id.second);
        if (!maybe) {
            LOG_ERROR("ValueError", "bad lha id for wilson conversion (unknown custom/base key)");
        }
        WCoefId coef = *maybe;
        QCDOrder order = parts[2] ? ((parts[2] -1) ? QCDOrder::NNLO : QCDOrder::NLO) : QCDOrder::LO;
        ContributionType part = parts[3] ? parts[3] -1 ? ContributionType::TOTAL : ContributionType::BSM : ContributionType::SM;

        std::pair<WCoefId, std::pair<QCDOrder, ContributionType>> ret;
        ret = {coef, {order, part}};

        return ret;
    }
    /**
     * @brief Builds a full FLHA LhaID for a given WCoefId, order and contribution.
     *
     * Uses the external base key attached to the id. If none is found,
     * an exception is thrown.
     *
     * @param id WCoef identifier.
     * @param q  QCD order.
     * @param c  Contribution type.
     * @return Corresponding LhaID.
     *
     * @throws std::runtime_error if @p id has no FLHA key.
     */
    static LhaID flha_full(const WCoefId& id, QCDOrder q, ContributionType c){
        auto k = external_of(id);
        if (!k) throw std::runtime_error("No FLHA key for "+id.str());
        return LhaID{k->first, k->second, int(q)-1, int(c)};
    }

    /**
     * @brief Returns the total number of Wilson coefficients in the mapping.
     */
    static size_t n_wilsons(){ return wcoef_mapping().size(); }

    /**
     * @brief Returns a static inverse mapping FLHA base pair -> WCoef.
     *
     * This is built once from wcoef_flha_mapping() and reused for all calls.
     *
     * @return Const reference to the inverse map.
     */
    static const std::map<std::pair<int,int>, WCoef>& inverse_flha_mapping() {
        static const std::map<std::pair<int,int>, WCoef> inv = []{
            std::map<std::pair<int,int>, WCoef> m;
            for (const auto& [e, ab] : wcoef_flha_mapping()) m.emplace(ab, e);
            return m;
        }();
        return inv;
    }

    /// @brief B-group (C1–C10) Wilson coefficients.
    static const std::vector<WCoef>& B_group(){
        static const std::vector<WCoef> g = {
            WCoef::C1,WCoef::C2,WCoef::C3,WCoef::C4,WCoef::C5,WCoef::C6,
            WCoef::C7,WCoef::C8,WCoef::C9,WCoef::C10
        }; return g;
    }

    /// @brief B'-group coefficients (primed + scalar primed).
    static const std::vector<WCoef>& B_prime_group(){
        static const std::vector<WCoef> g = {
            WCoef::CP1,WCoef::CP2,WCoef::CP3,WCoef::CP4,WCoef::CP5,WCoef::CP6,
            WCoef::CP7,WCoef::CP8,WCoef::CP9,WCoef::CP10,WCoef::CPQ1,WCoef::CPQ2
        }; return g;
    }

    /// @brief Scalar B-group coefficients.
    static const std::vector<WCoef>& B_scalar_group(){
        static const std::vector<WCoef> g = { WCoef::CQ1, WCoef::CQ2 }; return g;
    }

    /// @brief Charged-current Wilsons for b -> c l ν.
    static const std::vector<WCoef>& b_clnu_group(){
        static const std::vector<WCoef> g = { WCoef::C_V1_bc, WCoef::C_V2_bc, WCoef::C_S1_bc, WCoef::C_S2_bc, WCoef::C_T_bc }; return g;
    }

    /// @brief Charged-current Wilsons for b -> u l ν.
    static const std::vector<WCoef>& b_ulnu_group(){
        static const std::vector<WCoef> g = { WCoef::C_V1_bu, WCoef::C_V2_bu, WCoef::C_S1_bu, WCoef::C_S2_bu, WCoef::C_T_bu }; return g;
    }

    /// @brief Charged-current Wilsons for c -> s l ν.
    static const std::vector<WCoef>& c_slnu_group(){
        static const std::vector<WCoef> g = { WCoef::C_V1_cs, WCoef::C_V2_cs, WCoef::C_S1_cs, WCoef::C_S2_cs, WCoef::C_T_cs }; return g;
    }

    /// @brief Charged-current Wilsons for c -> d l ν.
    static const std::vector<WCoef>& c_dlnu_group(){
        static const std::vector<WCoef> g = { WCoef::C_V1_cd, WCoef::C_V2_cd, WCoef::C_S1_cd, WCoef::C_S2_cd, WCoef::C_T_cd }; return g;
    }

    /// @brief Charged-current Wilsons for s -> u l ν.
    static const std::vector<WCoef>& s_ulnu_group(){
        static const std::vector<WCoef> g = { WCoef::C_V1_su, WCoef::C_V2_su, WCoef::C_S1_su, WCoef::C_S2_su, WCoef::C_T_su }; return g;
    }

    /// @brief Charged-current Wilsons for d -> u l ν.
    static const std::vector<WCoef>& d_ulnu_group(){
        static const std::vector<WCoef> g = { WCoef::C_V1_du, WCoef::C_V2_du, WCoef::C_S1_du, WCoef::C_S2_du, WCoef::C_T_du}; return g;
    }

    /// @brief Kaon-sector Wilsons (K group).
    static const std::vector<WCoef>& k_group(){
        static const std::vector<WCoef> g = { WCoef::CK9, WCoef::CK10, WCoef::CKQ1, WCoef::CKQ2, WCoef::CK_L, WCoef::CPK9, WCoef::CPK10, WCoef::CPKQ1, WCoef::CPKQ2}; return g;
    }
    
    /// @brief Wilsons relevant for neutral meson mixing (B, K, D sectors).
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

    /**
     * @brief Returns the list of Wilson coefficients belonging to a WGroup.
     *
     * This is a convenience function on top of the various *_group() helpers.
     *
     * @param g Wilson group.
     * @return Vector of WCoef belonging to group @p g.
     *
     * @throws (via LOG_ERROR) if @p g is invalid.
     */
    static std::vector<WCoef> get_group(WGroup g){
        switch (g){
            case WGroup::B: return {B_group().begin(), B_group().end()};
            case WGroup::BPrime: return {B_prime_group().begin(), B_prime_group().end()};
            case WGroup::BScalar: return {B_scalar_group().begin(), B_scalar_group().end()};
            case WGroup::CC_bc: return {b_clnu_group().begin(), b_clnu_group().end()};
            case WGroup::CC_bu: return {b_ulnu_group().begin(), b_ulnu_group().end()};
            case WGroup::CC_cs: return {c_slnu_group().begin(), c_slnu_group().end()};
            case WGroup::CC_cd: return {c_dlnu_group().begin(), c_dlnu_group().end()};
            case WGroup::CC_su: return {s_ulnu_group().begin(), s_ulnu_group().end()};
            case WGroup::CC_du: return {d_ulnu_group().begin(), d_ulnu_group().end()};
            case WGroup::K:     return {k_group().begin(), k_group().end()};
            case WGroup::MESON_MIXING: return {meson_mixing_group().begin(), meson_mixing_group().end()};
            default: LOG_ERROR("Invalid WGroup","get_group couldn't find your group"); return {};
        }
    }

    /**
     * @brief Returns the WGroup to which a given Wilson coefficient belongs.
     *
     * Uses a static inverse map built once from the group definitions.
     *
     * @param c Wilson coefficient.
     * @return The corresponding WGroup.
     *
     * @throws std::out_of_range if the coefficient is not assigned to any group.
     */
    static WGroup group_of(WCoef c){
        static const std::map<WCoef, WGroup> inv = []{
            std::map<WCoef, WGroup> m;

            auto add = [&](WGroup g, const std::vector<WCoef>& v){
                for (auto e : v) {
                    // optional safety: detect duplicates across groups
                    if (m.find(e) != m.end())
                        throw std::runtime_error("WCoef appears in multiple WGroups");
                    m.emplace(e, g);
                }
            };

            add(WGroup::B,            B_group());
            add(WGroup::BPrime,       B_prime_group());
            add(WGroup::BScalar,      B_scalar_group());
            add(WGroup::CC_bc,        b_clnu_group());
            add(WGroup::CC_bu,        b_ulnu_group());
            add(WGroup::CC_cs,        c_slnu_group());
            add(WGroup::CC_cd,        c_dlnu_group());
            add(WGroup::CC_su,        s_ulnu_group());
            add(WGroup::CC_du,        d_ulnu_group());
            add(WGroup::K,            k_group());
            add(WGroup::MESON_MIXING, meson_mixing_group());

            return m;
        }();

        auto it = inv.find(c);
        if (it == inv.end()) throw std::out_of_range("WCoef not assigned to any WGroup");
        return it->second;
    }

    /**
     * @brief Optional version: returns std::nullopt if not in any group.
     */
    static std::optional<WGroup> group_of_opt(WCoef c){
        try { return group_of(c); }
        catch (...) { return std::nullopt; }
    }

    /**
     * @brief Convenience overload for identifiers.
     *
     * @throws std::runtime_error if id cannot be converted to an enum.
     */
    static WGroup group_of(const WCoefId& id){
        auto e = enum_of(id);
        if (!e) throw std::runtime_error("Unknown WCoefId: " + id.str());
        return group_of(*e);
    }
};

#endif