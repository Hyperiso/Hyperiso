// group_ids.hpp
#pragma once
#include "generic_mapper.hpp"
#include "Map.h"
#include "scaletype_ids.hpp"
#include "wilsonbasis_ids.hpp"

// Optional external-key map for groups (empty by default)
inline const std::map<WGroup, std::string>& group_external_mapping(){
    static const std::map<WGroup, std::string> empty{};
    return empty;
}

struct WGroupTag {};
using WGroupId = IdOf<WGroupTag>;

class GroupMapper
: public GenericMapperWithExt<
      WGroupTag,
      WGroup,
      std::string,
      std::hash<std::string>,
      group_mapping,
      group_external_mapping
  >
{
public:
    using Base = GenericMapperWithExt<
        WGroupTag, WGroup, std::string, std::hash<std::string>,
        group_mapping, group_external_mapping
    >;

    using Base::init_builtins;
    using Base::id_of;
    using Base::list_all;
    using Base::from_external;
    using Base::external_of;
    using Base::to_id;          // WGroup   -> WGroupId
    using Base::enum_of;        // WGroupId -> optional<WGroup>
    using Base::str;            // overloads: str(WGroupId) and str(WGroup)
    using Base::enum_elt;       // string -> WGroupId (keeps new runtime id API)
    using Base::enum_elt_legacy;// string -> WGroup  (legacy enum API)

    // Convenience: still allow registering with optional external string
    static bool register_custom(std::string canonical,
                                std::vector<std::string> aliases = {},
                                std::optional<std::string> ext = std::nullopt)
    {
        return Base::register_custom(std::move(canonical), std::move(aliases), std::move(ext));
    }

    // Domain-specific helpers (block name builders), preserved
    static std::string str(const WGroupId& gid, ScaleType s, WilsonBasis b = WilsonBasis::B_STANDARD){
        return str_impl(gid.str(), s, b);
    }
    static std::string str(WGroup g, ScaleType s, WilsonBasis b = WilsonBasis::B_STANDARD){
        return str_impl(Base::str(g), s, b);
    }

private:
    static std::string str_impl(const std::string& base, ScaleType s, WilsonBasis b){
        std::string out = base + "_" + ScaleTypeMapper::str(s);
        if (s == ScaleType::HADRONIC) out += "_" + WilsonBasisMapper::str(b);
        return out;
    }
};
