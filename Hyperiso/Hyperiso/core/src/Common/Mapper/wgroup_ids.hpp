#pragma once
#include "generic_mapper.hpp"
#include "Map.h"
#include "scaletype_ids.hpp"
#include "wilsonbasis_ids.hpp"

struct WGroupTag {};
using WGroupId = IdOf<WGroupTag>;

class GroupMapper {
    using Reg = DynamicRegistry<WGroupTag, std::string, std::hash<std::string>>;
    static Reg& reg(){ static Reg R; return R; }
    static bool& inited(){ static bool f=false; return f; }
    static void ensure_init(){
        if (inited()) return;
        for (auto& [e,name] : group_mapping())
            reg().register_id(WGroupId(name), {}, /*ext=*/std::nullopt, /*builtin=*/true);
        inited() = true;
    }
public:
    static void init_builtins(){ ensure_init(); }
    static WGroupId id_of(std::string_view s){ ensure_init(); auto r=reg().find(s); if(!r) throw std::runtime_error("unknown group"); return *r; }
    static WGroupId enum_elt(std::string_view s){ return id_of(s); }
    static bool register_custom(std::string canonical, std::vector<std::string> aliases={}, std::optional<std::string> ext=std::nullopt){
        ensure_init(); return reg().register_id(WGroupId(std::move(canonical)), std::move(aliases), std::move(ext), false);
    }
    static std::vector<WGroupId> list_all(){ ensure_init(); return reg().list_all(); }
    static std::string str(const WGroupId& id){ return id.str(); }
    static std::string str(WGroup e){ return group_mapping().at(e); }

    static std::string str(const WGroupId& gid, ScaleType s, WilsonBasis b = WilsonBasis::B_STANDARD){
        std::string out = gid.str() + "_" + ScaleTypeMapper::str(s);
        if (s == ScaleType::HADRONIC) out += "_" + WilsonBasisMapper::str(b);
        return out;
    }

    static WGroupId to_id(WGroup e){ return WGroupId(str(e)); }
    static std::optional<WGroup> enum_of(const WGroupId& id){
        auto key=nk(id.str()); for(auto&[e,name]:group_mapping()) if(nk(name)==key) return e; return std::nullopt;
    }
    static std::optional<WGroupId> from_external(std::string_view k){ ensure_init(); return reg().find_by_external(std::string(k)); }
    static std::optional<std::string> external_of(const WGroupId& id){ ensure_init(); return reg().external_of(id); }

    static std::string str(WGroup g, ScaleType s, WilsonBasis b=WilsonBasis::B_STANDARD){
        std::string out = str(g) + "_" + ScaleTypeMapper::str(s);
        if (s == ScaleType::HADRONIC) out += "_" + WilsonBasisMapper::str(b);
        return out;
    }
};
