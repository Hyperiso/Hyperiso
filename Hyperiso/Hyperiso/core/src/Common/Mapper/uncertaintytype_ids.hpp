// #pragma once
// #include "dynamic_registry.hpp"
// #include <map>
// #include <optional>
// #include <string_view>

// ==== Tables compile-time (legacy, header-only) ====
inline const std::map<UncertaintyType, std::string>& uncertaintytype_mapping() {
    static const std::map<UncertaintyType, std::string> m = {
        {UncertaintyType::STAT,     "STATISTICAL"},
        {UncertaintyType::SYST,     "SYSTEMATICS"},
        {UncertaintyType::COMBINED, "COMBINED"},
    };
    return m;
}

// inline const std::map<std::string, UncertaintyType>& inverse_uncertaintytype_mapping() {
//     static const std::map<std::string, UncertaintyType> inv = []{
//         std::map<std::string, UncertaintyType> r;
//         for (auto& [k,v] : uncertaintytype_mapping()) r.emplace(v, k);
//         return r;
//     }();
//     return inv;
// }

// inline const std::map<UncertaintyType, DataType>& uncertainty_data_type_mapping() {
//     static const std::map<UncertaintyType, DataType> m = {
//         {UncertaintyType::STAT,     DataType::STD_STAT},
//         {UncertaintyType::SYST,     DataType::STD_SYST},
//         {UncertaintyType::COMBINED, DataType::STD_COMBINED},
//     };
//     return m;
// }

// // ==== Runtime IDs ====
// struct UncTag {};
// using UncertaintyTypeId = SymbolId<UncTag>;

// inline UncertaintyTypeId to_id(UncertaintyType e){
//     return UncertaintyTypeId(uncertaintytype_mapping().at(e));
// }

// class UncertaintyTypeMapper {
//     using Reg = DynamicRegistry<UncTag, void, void>;
//     static Reg& reg(){ static Reg R; return R; }

//     static void ensure_init() {
//         static bool ready = false;
//         if (!ready) { init_builtins(); ready = true; }
//     }

// public:
//     // Init (idempotent)
//     static void init_builtins(){
//         static bool done = false;
//         if (done) return; done = true;
//         for (auto& [e,name] : uncertaintytype_mapping())
//             reg().register_id(UncertaintyTypeId(name), {}, /*builtin=*/true);
//     }

//     // === API runtime ===
//     static UncertaintyTypeId enum_elt(std::string_view s){
//         ensure_init();
//         auto r = reg().find(s);
//         if (!r) throw std::runtime_error("Unknown UncertaintyType: " + std::string(s));
//         return *r;
//     }
//     static bool register_custom(const std::string& n, std::vector<std::string> a = {}) {
//         ensure_init();
//         return reg().register_id(UncertaintyTypeId(n), std::move(a), /*builtin=*/false);
//     }
//     static std::vector<UncertaintyTypeId> list_all(){ ensure_init(); return reg().list_all(); }

//     static std::string str(const UncertaintyTypeId& id){ return id.str(); }
//     static std::string str(UncertaintyType e){ return to_id(e).str(); }

//     // id runtime -> enum builtin (si possible)
//     static std::optional<UncertaintyType> enum_of(const UncertaintyTypeId& id){
//         auto it = inverse_uncertaintytype_mapping().find(id.str());
//         if (it != inverse_uncertaintytype_mapping().end()) return it->second;
//         return std::nullopt;
//     }

//     // ==== API legacy réintroduite ====
//     static const std::map<UncertaintyType, std::string>& mapping() {
//         return uncertaintytype_mapping();
//     }
//     static const std::map<std::string, UncertaintyType>& inverse_mapping() {
//         return inverse_uncertaintytype_mapping();
//     }
//     static const std::map<UncertaintyType, DataType>& data_type_mapping() {
//         return uncertainty_data_type_mapping();
//     }
//     static DataType d_type(UncertaintyType u_type) {
//         return data_type_mapping().at(u_type);
//     }
//     // Optionnel: version qui accepte un id runtime (fail si custom inconnu)
//     static DataType d_type(const UncertaintyTypeId& id) {
//         if (auto e = enum_of(id)) return d_type(*e);
//         throw std::runtime_error("No DataType mapping for custom uncertainty type: " + id.str());
//     }
// };

#pragma once
#include "generic_mapper.hpp"
#include "Map.h"
#include "GeneralEnum.h"

struct UncTag {};
using UncertaintyTypeId = IdOf<UncTag>;

class UncertaintyTypeMapper
: public GenericMapperNoExt<UncTag, UncertaintyType, uncertaintytype_mapping>
{
public:
    using Base = GenericMapperNoExt<UncTag, UncertaintyType, uncertaintytype_mapping>;
    using Base::str;

    static const std::map<UncertaintyType, DataType>& data_type_mapping(){
        static const std::map<UncertaintyType, DataType> m = {
            {UncertaintyType::STAT,     DataType::STD_STAT},
            {UncertaintyType::SYST,     DataType::STD_SYST},
            {UncertaintyType::COMBINED, DataType::STD_COMBINED},
        };
        return m;
    }
    static DataType d_type(UncertaintyType t){ return data_type_mapping().at(t); }
};
