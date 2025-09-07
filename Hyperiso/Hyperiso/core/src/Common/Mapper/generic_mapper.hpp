#pragma once
#include <map>
#include <optional>
#include <string_view>
#include <vector>
#include "dynamic_registry.hpp"
#include "Map.h"

// Utility
template<class Tag>
using IdOf = SymbolId<Tag>;

inline std::string nk(std::string_view s){ return normalize_key(s); }

// ------------------------------ Without External Keys ------------------------------
template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
class GenericMapperNoExt {
    using Reg = DynamicRegistry<Tag, void, void>;
    static Reg& reg(){ static Reg R; return R; }
    static bool& inited(){ static bool f=false; return f; }
    static void ensure_init(){
        if (inited()) return;
        for (auto& [e,name] : MapFn())
            reg().register_id(IdOf<Tag>(name), {}, /*builtin=*/true);
        inited() = true;
    }
public:
    // Explicit initialization (optional)
    static void init_builtins(){ ensure_init(); }

    // String to Id
    static IdOf<Tag> id_of(std::string_view s){ ensure_init(); auto r=reg().find(s); if(!r) throw std::runtime_error("unknown"); return *r; }
    static IdOf<Tag> enum_elt(std::string_view s){ return id_of(s); } // alias compat

    // For custom register
    static bool register_custom(std::string canonical, std::vector<std::string> aliases={}){
        ensure_init(); return reg().register_id(IdOf<Tag>(std::move(canonical)), std::move(aliases), /*builtin=*/false);
    }

    // List all member of enum
    static std::vector<IdOf<Tag>> list_all(){ ensure_init(); return reg().list_all(); }

    // Id/enum to str
    static std::string str(const IdOf<Tag>& id){ return id.str(); }
    static std::string str(EnumT e){ return MapFn().at(e); }

    // bridge enum <-> id
    static IdOf<Tag> to_id(EnumT e){ return IdOf<Tag>(str(e)); }
    static std::optional<EnumT> enum_of(const IdOf<Tag>& id){
        auto key = nk(id.str());
        for (auto& [e,name] : MapFn()) if (nk(name)==key) return e;
        return std::nullopt;
    }
};

// ------------------------------ With External Keys ------------------------------
template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
class GenericMapperWithExt {
    using Reg = DynamicRegistry<Tag, ExternalKey, Hash>;
    static Reg& reg(){ static Reg R; return R; }
    static bool& inited(){ static bool f=false; return f; }
    static void ensure_init(){
        if (inited()) return;
        for (auto& [e,name] : MapFn()){
            auto it = ExtFn().find(e);
            std::optional<ExternalKey> ext = (it==ExtFn().end() ? std::nullopt : std::optional<ExternalKey>(it->second));
            reg().register_id(IdOf<Tag>(name), {}, ext, /*builtin=*/true);
        }
        inited() = true;
    }
public:
    static void init_builtins(){ ensure_init(); }

    // --- string -> Id
    static IdOf<Tag> id_of(std::string_view s){ ensure_init(); auto r=reg().find(s); if(!r) throw std::runtime_error("unknown"); return *r; }
    static IdOf<Tag> enum_elt(std::string_view s){ return id_of(s); } // alias compat

    // --- ext <-> Id
    static std::optional<IdOf<Tag>> from_external(const ExternalKey& k){ ensure_init(); return reg().find_by_external(k); }
    static std::optional<ExternalKey> external_of(const IdOf<Tag>& id){ ensure_init(); return reg().external_of(id); }

    // External registering
    static bool register_custom(std::string canonical, std::vector<std::string> aliases={}, std::optional<ExternalKey> ext=std::nullopt){
        ensure_init(); return reg().register_id(IdOf<Tag>(std::move(canonical)), std::move(aliases), std::move(ext), /*builtin=*/false);
    }

    static bool set_external(const IdOf<Tag>& id, const ExternalKey& k){
        ensure_init();
        return reg().update_external(id, k);
    }

    // Listing
    static std::vector<IdOf<Tag>> list_all(){ ensure_init(); return reg().list_all(); }

    static std::string str(const IdOf<Tag>& id){ return id.str(); }
    static std::string str(EnumT e){ return MapFn().at(e); }

    static IdOf<Tag> to_id(EnumT e){ return IdOf<Tag>(str(e)); }
    static std::optional<EnumT> enum_of(const IdOf<Tag>& id){
        auto key = nk(id.str());
        for (auto& [e,name] : MapFn()) if (nk(name)==key) return e;
        return std::nullopt;
    }
};
