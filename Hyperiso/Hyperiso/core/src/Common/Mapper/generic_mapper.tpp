

template<class Tag, class EnumT, const std::map<EnumT, std::string>& (*MapFn)()>
typename GenericMapperNoExt<Tag, EnumT, MapFn>::Reg&
GenericMapperNoExt<Tag, EnumT, MapFn>::reg(){ static Reg R; return R; }

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
bool& GenericMapperNoExt<Tag, EnumT,
          MapFn>::inited(){ static bool f=false; return f; }

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
void GenericMapperNoExt<Tag, EnumT,
          MapFn>::ensure_init(){
    if (inited()) return;
    for (auto& [e,name] : MapFn())
        reg().register_id(IdOf<Tag>(name), {}, /*builtin=*/true);
    inited() = true;
}

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
void GenericMapperNoExt<Tag, EnumT,
         MapFn>::init_builtins(){ ensure_init(); }

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
EnumT GenericMapperNoExt<Tag, EnumT,
         MapFn>::enum_elt_legacy(std::string_view s){
    auto key = nk(s);
    for (auto& [e,name] : MapFn()) if (nk(name)==key) return e;
    throw std::out_of_range("Unknown name: " + std::string(s));
}

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
std::vector<std::string> GenericMapperNoExt<Tag, EnumT,
         MapFn>::get_str(){
    std::vector<std::string> out;
    out.reserve(MapFn().size());
    for (const auto& [e,name] : MapFn()) out.push_back(name);
    return out;
}

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
std::vector<EnumT> GenericMapperNoExt<Tag, EnumT,
          MapFn>::get_enum(){
    std::vector<EnumT> out;
    out.reserve(MapFn().size());
    for (const auto& [e,name] : MapFn()) out.push_back(e);
    return out;
}

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
std::vector<std::string> GenericMapperNoExt<Tag, EnumT,
          MapFn>::get_str_all(){
    ensure_init();
    std::vector<std::string> out;
    for (auto& id : reg().list_all()) out.push_back(id.str());
    return out;
}

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
IdOf<Tag> GenericMapperNoExt<Tag, EnumT,
          MapFn>::id_of(std::string_view s){ ensure_init(); auto r=reg().find(s); if(!r) throw std::runtime_error("unknown"); return *r; }

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
IdOf<Tag> GenericMapperNoExt<Tag, EnumT,
          MapFn>::enum_elt(std::string_view s){ return id_of(s); } // alias compat

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
bool GenericMapperNoExt<Tag, EnumT,
         MapFn>::register_custom(std::string canonical, std::vector<std::string> aliases){
    ensure_init(); return reg().register_id(IdOf<Tag>(std::move(canonical)), std::move(aliases), /*builtin=*/false);
}

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
std::vector<IdOf<Tag>> GenericMapperNoExt<Tag, EnumT,
          MapFn>::list_all(){ ensure_init(); return reg().list_all(); }

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
std::string GenericMapperNoExt<Tag, EnumT,
          MapFn>::str(const IdOf<Tag>& id){ return id.str(); }

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
std::string GenericMapperNoExt<Tag, EnumT,
          MapFn>::str(EnumT e){ return MapFn().at(e); }

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
IdOf<Tag> GenericMapperNoExt<Tag, EnumT,
          MapFn>::to_id(EnumT e){ return IdOf<Tag>(str(e)); }

template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
std::optional<EnumT> enum_of(const IdOf<Tag>& id){
    auto key = nk(id.str());
    for (auto& [e,name] : MapFn()) if (nk(name)==key) return e;
    return std::nullopt;
}


template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
typename GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash, MapFn, ExtFn>::Reg&
GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash, MapFn, ExtFn>::reg(){ static Reg R; return R; }

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
bool& GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
         MapFn,
          ExtFn>::inited(){ static bool f=false; return f; }

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
void GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::ensure_init(){
    if (inited()) return;
    for (auto& [e,name] : MapFn()){
        auto it = ExtFn().find(e);
        std::optional<ExternalKey> ext = (it==ExtFn().end() ? std::nullopt : std::optional<ExternalKey>(it->second));
        reg().register_id(IdOf<Tag>(name), {}, ext, /*builtin=*/true);
    }
    inited() = true;
}

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
void GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::init_builtins(){ ensure_init(); }

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
std::vector<std::string> GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::get_str(){
    std::vector<std::string> out;
    out.reserve(MapFn().size());
    for (const auto& [e,name] : MapFn()) out.push_back(name);
    return out;
}

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
std::vector<EnumT> GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::get_enum(){
    std::vector<EnumT> out;
    out.reserve(MapFn().size());
    for (const auto& [e,name] : MapFn()) out.push_back(e);
    return out;
}

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
std::vector<std::string> GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
         ExtFn>::get_str_all(){
    ensure_init();
    std::vector<std::string> out;
    for (auto& id : reg().list_all()) out.push_back(id.str());
    return out;
}

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
EnumT GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
         MapFn,
          ExtFn>::enum_elt_legacy(std::string_view s){
    auto key = nk(s);
    for (auto& [e,name] : MapFn()) if (nk(name)==key) return e;
    throw std::out_of_range("Unknown name: " + std::string(s));
}


template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
IdOf<Tag> GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
         MapFn,
          ExtFn>::id_of(std::string_view s){ ensure_init(); auto r=reg().find(s); if(!r) throw std::runtime_error("unknown"); return *r; }

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
IdOf<Tag> GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::enum_elt(std::string_view s){ return id_of(s); } // alias compat

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
std::optional<IdOf<Tag>> GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::from_external(const ExternalKey& k){ ensure_init(); return reg().find_by_external(k); }

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
std::optional<ExternalKey> GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::external_of(const IdOf<Tag>& id){ ensure_init(); return reg().external_of(id); }

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
bool GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::register_custom(std::string canonical, std::vector<std::string> aliases, std::optional<ExternalKey> ext){
    ensure_init(); return reg().register_id(IdOf<Tag>(std::move(canonical)), std::move(aliases), std::move(ext), /*builtin=*/false);
}

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
bool GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::set_external(const IdOf<Tag>& id, const ExternalKey& k){
    ensure_init();
    return reg().update_external(id, k);
}

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
std::vector<IdOf<Tag>> GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::list_all(){ ensure_init(); return reg().list_all(); }

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
std::string GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::str(const IdOf<Tag>& id){ return id.str(); }

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
std::string GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::str(EnumT e){ return MapFn().at(e); }

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
IdOf<Tag> GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
         MapFn,
         ExtFn>::to_id(EnumT e){ return IdOf<Tag>(str(e)); }

template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
std::optional<EnumT> GenericMapperWithExt<Tag, EnumT, ExternalKey, Hash,
          MapFn,
          ExtFn>::enum_of(const IdOf<Tag>& id){
    auto key = nk(id.str());
    for (auto& [e,name] : MapFn()) if (nk(name)==key) return e;
    return std::nullopt;
}