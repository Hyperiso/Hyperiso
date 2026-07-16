#include "LhaParser.h"

namespace {

bool is_ignored_metadata_block(std::string_view name) {
    // Informational blocks do not contain numerical model parameters and may
    // have free-form string payloads.  Treating them as unsupported numerical
    // blocks only creates noisy warnings for otherwise valid spectrum files.
    return name == "SPINFO" || name == "spinfo" ||
           name == "DCINFO" || name == "dcinfo";
}

} // namespace

void LhaParser::addBlock(std::map<BlockName, std::shared_ptr<LhaBlock>>& blocks, const BlockName& id, const std::vector<std::vector<std::string>>& lines) const {
    auto block = std::make_shared<LhaBlock>(findPrototype(id));
    LOG_DEBUG(id);
    block->readData(lines);
    BlockName id_ci = block->getPrototype().blockName;
    id_ci.to_upper();
    blocks.insert(std::pair(id_ci, std::move(block)));
}

std::vector<Token> LhaParser::tokenize(const std::string &src) const {
    int cLine = 0;
    int cCol = 0;

    auto rit = std::sregex_iterator(src.begin(), src.end(), analyzer_rx);
    auto rend = std::sregex_iterator();

    std::vector<Token> tokens;
    while (rit != rend) {
        std::smatch m = *rit;
        size_t group_index = m.size();

        for (size_t idx = 1; idx < m.size(); ++idx) {
            if (m[idx].matched) {
                group_index = idx - 1;
                break;
            }
        }

        auto tokenType = static_cast<TokenType>(group_index);
        auto value = m[group_index + 1].str();

        if (tokenType == TokenType::NEWLINE) {
            ++cLine;
            cCol = 0;
        } else if (tokenType != TokenType::SKIP && value != "") {
            tokens.emplace_back(Token{tokenType, value, cLine, cCol});
            ++cCol;
        }

        ++rit;
    }

    return tokens;
}


std::map<BlockName, std::vector<std::vector<std::string>>>
LhaParser::parse_tokens(std::vector<Token> tokens, bool comments) const
{
    bool newBlock       = false;
    bool hasGlobalScale = false;
    bool isQ            = false;
    bool skipBlock      = false;
    bool decay          = false;

    std::string globalQ;
    BlockName   cBlock;
    int         cCol = INT_MAX;

    std::map<BlockName, std::vector<std::vector<std::string>>> rawBlocks;
    auto current_block_it = rawBlocks.end();

    for (const Token& t : tokens) {
        if (newBlock) {
            auto prototype = this->findPrototype(t.value);
            if (prototype.blockName != "") {
                LOG_DEBUG("LHA reader: Block ", prototype.blockName, " found.");

                cBlock        = prototype.blockName;
                hasGlobalScale= prototype.globalScale;
                skipBlock     = false;
                decay         = false;
                globalQ.clear();
                cCol = INT_MAX;

                auto [it, inserted] =
                    rawBlocks.emplace(cBlock, std::vector<std::vector<std::string>>{});
                current_block_it = it;
            }
            else if (decay) {
                LOG_DEBUG("LHA reader: Decay block found. Skipping.");
                skipBlock        = true;
                hasGlobalScale   = false;
                current_block_it = rawBlocks.end();
                decay            = false;
            }
            else {
                if (is_ignored_metadata_block(t.value)) {
                    LOG_DEBUG("LHA reader: Metadata block ", t.value, " skipped.");
                } else {
                    LOG_WARN("LHA reader: Unknown block " + t.value + " encountered. Skipping.");
                }
                skipBlock        = true;
                hasGlobalScale   = false;
                current_block_it = rawBlocks.end();
            }
            newBlock = false;
        }
        else if (t.type == TokenType::BLOCK) {
            newBlock = true;
        }
        else if (t.type == TokenType::DECAY) {
            decay    = true;
            newBlock = true;
        }
        else if (!comments && t.type == TokenType::COMMENT) {
            continue;
        }
        else if (hasGlobalScale && t.type == TokenType::WORD && !skipBlock) {
            isQ = (t.value == "Q=" || t.value == "q=");
        }
        else if (!skipBlock && current_block_it != rawBlocks.end()) {
            if (isQ) {
                globalQ = t.value;
                isQ     = false;
            }
            else {
                LOG_VERBOSE("Token : [", static_cast<int>(t.type), ", ", t.value, "]");
                auto& rows = current_block_it->second;

                if (t.col <= cCol) {
                    rows.emplace_back(std::vector<std::string>{});
                    if (hasGlobalScale)
                        rows.back().emplace_back(!globalQ.empty() ? globalQ : "-1");
                }
                rows.back().emplace_back(t.value);
                cCol = t.col;
            }
        }
    }

    return rawBlocks;
}

Prototype LhaParser::findPrototype(BlockName name) const {
    for (const auto& p : blockPrototypes) {
        name.to_upper();
        if (p.blockName == name)
            return p;
    }
    return Prototype{""};
}

std::shared_ptr<DBNode> LhaParser::readFromFile(const std::string& input_file) const {
    std::ifstream file(input_file.data());
    std::stringstream buffer;
    buffer << file.rdbuf();
    return this->parse(buffer.str());
}

std::shared_ptr<DBNode> LhaParser::parse(const std::string &src) const {
    auto rawBlocks = this->parse_tokens(this->tokenize(src));
    std::map<BlockName, std::shared_ptr<LhaBlock>> blocks;
    for (auto &[id, lines] : rawBlocks) {
        addBlock(blocks, id, lines);
    }     
    LOG_DEBUG("LHA file parsed.");

    return this->toDBNode(blocks);
}


void LhaParser::set_prototypes(const std::unordered_set<Prototype> &prototypes)
{
    this->blockPrototypes = prototypes;
}


std::shared_ptr<DBNode> LhaParser::toDBNode(std::map<BlockName, std::shared_ptr<LhaBlock>> blocks) const {
    DBNode root;
    
    for (const auto& [blockName, blockPtr] : blocks) {
        auto block_node = blockPtr->toDBNode();
        
        auto group = block_node->getGroup({blockName});
        root.setGroup({blockName}, group);
        
        if (block_node->contains("scale")) {
            auto sval = block_node->get("scale");
            if (std::holds_alternative<double>(sval)) {
                root.set(std::get<double>(sval), blockName, "scale");
            } else if (std::holds_alternative<int>(sval)) {
                root.set(static_cast<double>(std::get<int>(sval)), blockName, "scale");
            } else {
                LOG_WARN("Expected numeric 'scale' for block ", blockName, " but got variant index ", sval.index());
            }
        }
    }
    
    return std::make_shared<DBNode>(root);
}

static double value_to_double(const DBNode::Value& v, double fallback = 0.0)
{
    if (std::holds_alternative<double>(v)) return std::get<double>(v);
    if (std::holds_alternative<int>(v))    return static_cast<double>(std::get<int>(v));
    return fallback;
}

static bool value_is_node_ptr(const DBNode::Value& v)
{
    return std::holds_alternative<std::shared_ptr<DBNode>>(v);
}

static std::shared_ptr<DBNode> value_to_node_ptr(const DBNode::Value& v)
{
    if (!value_is_node_ptr(v)) return nullptr;
    return std::get<std::shared_ptr<DBNode>>(v);
}

static std::vector<std::string> split_tokens(const std::string& s)
{
    std::string tmp = s;
    std::replace(tmp.begin(), tmp.end(), '_', ' ');
    std::istringstream iss(tmp);

    std::vector<std::string> out;
    std::string tok;
    while (iss >> tok) out.push_back(tok);
    return out;
}

static std::vector<std::string> split_ws(const std::string& s)
{
    std::istringstream iss(s);
    std::vector<std::string> out;
    std::string tok;
    while (iss >> tok) out.push_back(tok);
    return out;
}

static bool read_ew_scale_Q(const std::shared_ptr<DBNode>& root, double& Q_out)
{
    if (!root || !root->contains(BlockName("EW_SCALE"))) return false;

    auto g = root->getGroup({BlockName("EW_SCALE")});

    auto it = g.find(BlockName("1"));
    if (it != g.end()) {
        auto node = value_to_node_ptr(it->second);
        if (node && node->contains("central_value")) {
            Q_out = value_to_double(node->get("central_value"), 0.0);
            return true;
        }
    }

    for (const auto& kv : g) {
        auto node = value_to_node_ptr(kv.second);
        if (node && node->contains("central_value")) {
            Q_out = value_to_double(node->get("central_value"), 0.0);
            return true;
        }
    }

    return false;
}

static bool starts_with(const std::string& s, const std::string& pre)
{
    return s.size() >= pre.size() && s.compare(0, pre.size(), pre) == 0;
}

static std::map<BlockName, DBNode::Value>
build_fwcoef_group_from_ewscale_blocks(const std::shared_ptr<DBNode>& root,
                                       const std::vector<std::string>& prefixes,
                                       bool& usedExternalBlocks)
{
    usedExternalBlocks = false;
    std::map<BlockName, DBNode::Value> out;

    if (!root) return out;

    for (const auto& k : root->get_keys()) {
        std::ostringstream oss; oss << k;
        const std::string ks = oss.str();

        if (!ends_with(ks, "_EW_SCALE")) continue;

        bool okPrefix = false;
        for (const auto& pre : prefixes) {
            if (starts_with(ks, pre)) { okPrefix = true; break; }
        }
        if (!okPrefix) continue;

        usedExternalBlocks = true;

        auto g = root->getGroup({k});
        for (const auto& kv : g) {
            if (kv.first == BlockName("scale")) continue;
            out[kv.first] = kv.second;
        }
    }

    return out;
}

static std::string blockname_primary(const BlockName& bn)
{
    std::ostringstream oss;
    oss << bn;
    return oss.str();
}

static const Prototype* find_proto(const std::unordered_set<Prototype>& protos, const BlockName& name)
{
    for (const auto& p : protos) {
        if (p.blockName == name) return &p;
    }
    return nullptr;
}

static void write_block_header(std::ostream& os, const Prototype& proto, double Q, bool hasQ)
{
    os << "BLOCK " << blockname_primary(proto.blockName);
    if (proto.globalScale && hasQ) {
        os << " Q= " << std::setprecision(8) << std::scientific << Q;
    }
    os << "\n";
}

static void write_value(std::ostream& os, double v)
{
    os << std::setprecision(8) << std::scientific << v;
}

static std::string pad_left_zeros(std::string s, size_t width)
{
    if (!s.empty() && std::all_of(s.begin(), s.end(), [](unsigned char c){ return std::isdigit(c); })) {
        if (s.size() < width) s = std::string(width - s.size(), '0') + s;
    }
    return s;
}

static void write_fwcoef_like(std::ostream& os,
                              const Prototype& proto,
                              const std::map<BlockName, DBNode::Value>& group,
                              double Q_from_ewscale,
                              bool hasQ_from_ewscale)
{
    write_block_header(os, proto, Q_from_ewscale, hasQ_from_ewscale);
    os << "#id                 order   M   value           comment\n";

    for (const auto& kv : group) {
        if (kv.first == BlockName("scale")) continue;

        auto node = value_to_node_ptr(kv.second);
        if (!node || !node->contains("central_value")) continue;

        const double val = value_to_double(node->get("central_value"), 0.0);

        const auto toks = split_tokens(blockname_primary(kv.first));
        if (toks.size() < 4) continue;

        std::string id1 = toks[0];
        std::string id2 = toks[1];
        std::string ord = toks[2];
        std::string cty = toks[3];

        if (id1.size() == 7) id1 = pad_left_zeros(id1, 8);

        int ord_i = 0;
        try { ord_i = std::stoi(ord); } catch (...) { ord_i = 0; }

        os << std::setw(8) << id1 << "  "
           << std::setw(6) << id2 << "  "
           << std::setfill('0') << std::setw(2) << ord_i << std::setfill(' ') << "  "
           << std::setw(2) << cty << "  ";
        write_value(os, val);
        os << "\n";
    }

    os << "\n";
}


static void write_fobs_like(std::ostream& os,
                            const Prototype& proto,
                            const std::map<BlockName, DBNode::Value>& group)
{
    write_block_header(os, proto, 0.0, false);
    os << "# ParentPDG type    value           q  bin_low  bin_high   NDA ID1 ID2 ID3 ... comment\n";

    for (const auto& kv : group) {
        if (kv.first == BlockName("scale")) continue;

        auto node = value_to_node_ptr(kv.second);
        if (!node || !node->contains("central_value")) continue;

        const double val = value_to_double(node->get("central_value"), 0.0);
        const double q   = node->contains("scale") ? value_to_double(node->get("scale"), 0.0) : 0.0;

        const double bin_low  = node->contains("bin_low")  ? value_to_double(node->get("bin_low"), 0.0)  : 0.0;
        const double bin_high = node->contains("bin_high") ? value_to_double(node->get("bin_high"), 0.0) : 0.0;

        const auto toks = split_tokens(blockname_primary(kv.first));
        if (toks.size() < 3) continue;

        os << std::setw(6) << toks[0] << "  "
           << std::setw(4) << toks[1] << "  ";
        write_value(os, val);

        os << "  " << std::setw(2) << static_cast<int>(q)
           << "  " << std::setw(7) << std::fixed << std::setprecision(0) << bin_low
           << "  " << std::setw(8) << std::fixed << std::setprecision(0) << bin_high
           << std::scientific;

        os << "  " << std::setw(3) << toks[2];
        for (size_t i = 3; i < toks.size(); ++i) {
            os << "  " << std::setw(4) << toks[i];
        }
        os << "\n";
    }

    os << "\n";
}

static void write_fmass_like(std::ostream& os,
                             const Prototype& proto,
                             const std::map<BlockName, DBNode::Value>& group)
{
    write_block_header(os, proto, 0.0, false);
    os << "# PDG_code  mass        scheme  Q   particle\n";
    for (const auto& kv : group) {
        const auto& key = kv.first;
        if (key == BlockName("scale")) continue;

        auto node = value_to_node_ptr(kv.second);
        if (!node || !node->contains("central_value")) continue;

        const double val = value_to_double(node->get("central_value"), 0.0);
        const double q   = node->contains("scale") ? value_to_double(node->get("scale"), 0.0) : 0.0;

        auto toks = split_tokens(blockname_primary(key));
        std::string pdg = toks.size() > 0 ? toks[0] : "0";
        std::string scheme = toks.size() > 1 ? toks[1] : "0";

        os << std::setw(6) << pdg << "  ";
        write_value(os, val);
        os << "  " << std::setw(2) << scheme
           << "  " << std::setw(2) << static_cast<int>(q)
           << "\n";
    }

    os << "\n";
}

static void write_fconst_like(std::ostream& os,
                              const Prototype& proto,
                              const std::map<BlockName, DBNode::Value>& group)
{
    write_block_header(os, proto, 0.0, false);
    os << "# PDG_code  number  decay_constant  scheme  scale   particle\n";

    for (const auto& kv : group) {
        const auto& key = kv.first;
        if (key == BlockName("scale")) continue;

        auto node = value_to_node_ptr(kv.second);
        if (!node || !node->contains("central_value")) continue;

        const double val   = value_to_double(node->get("central_value"), 0.0);
        const double scale = node->contains("scale") ? value_to_double(node->get("scale"), 1.0) : 1.0;

        auto toks = split_tokens(blockname_primary(key));
        std::string pdg    = toks.size() > 0 ? toks[0] : "0";
        std::string number = toks.size() > 1 ? toks[1] : "1";
        std::string scheme = toks.size() > 2 ? toks[2] : "0";

        os << std::setw(6) << pdg << "  "
           << std::setw(3) << number << "  ";
        write_value(os, val);
        os << "  " << std::setw(2) << scheme << "  "
           << std::fixed << std::setprecision(1) << scale
           << std::scientific << "\n";
    }

    os << "\n";
}

static void write_fconstratio_like(std::ostream& os,
                                   const Prototype& proto,
                                   const std::map<BlockName, DBNode::Value>& group)
{
    write_block_header(os, proto, 0.0, false);
    os << "# PDG_code1 code2   nb1 nb2 ratio           scheme  scale   comment\n";

    for (const auto& kv : group) {
        const auto& key = kv.first;
        if (key == BlockName("scale")) continue;

        auto node = value_to_node_ptr(kv.second);
        if (!node || !node->contains("central_value")) continue;

        const double val   = value_to_double(node->get("central_value"), 0.0);
        const double scale = node->contains("scale") ? value_to_double(node->get("scale"), 1.0) : 1.0;

        auto toks = split_tokens(blockname_primary(key));
        if (toks.size() < 4) continue;

        std::string scheme = toks.size() > 4 ? toks[4] : "0";

        os << std::setw(6) << toks[0] << "  "
           << std::setw(6) << toks[1] << "  "
           << std::setw(2) << toks[2] << "  "
           << std::setw(2) << toks[3] << "  ";
        write_value(os, val);
        os << "  " << std::setw(2) << scheme << "  "
           << std::fixed << std::setprecision(1) << scale
           << std::scientific << "\n";
    }

    os << "\n";
}

static void write_fbag_like(std::ostream& os,
                            const Prototype& proto,
                            const std::map<BlockName, DBNode::Value>& group)
{
    write_fconst_like(os, proto, group);
}

static void write_generic_block(std::ostream& os,
                                const Prototype& proto,
                                const std::map<BlockName, DBNode::Value>& group)
{
    bool hasQ = false;
    double Q = 0.0;

    auto itScale = group.find(BlockName("scale"));
    if (itScale != group.end()) {
        Q = value_to_double(itScale->second, 0.0);
        hasQ = true;
    }

    write_block_header(os, proto, Q, hasQ);

    for (const auto& kv : group) {
        if (kv.first == BlockName("scale")) continue;

        auto node = value_to_node_ptr(kv.second);
        if (!node || !node->contains("central_value")) continue;

        const double val = value_to_double(node->get("central_value"), 0.0);

        auto toks = split_tokens(blockname_primary(kv.first));

        for (const auto& t : toks) {
            os << std::setw(6) << t << "  ";
        }
        write_value(os, val);
        os << "\n";
    }

    os << "\n";
}

static std::map<BlockName, DBNode::Value> build_imaginary_group(
    const std::map<BlockName, DBNode::Value>& group
)
{
    std::map<BlockName, DBNode::Value> imaginary_group;
    bool has_imaginary_entries = false;

    for (const auto& [key, value] : group) {
        if (key == BlockName("scale")) {
            continue;
        }

        auto node = value_to_node_ptr(value);
        if (!node || !node->contains("imaginary_value")) {
            continue;
        }

        const double imaginary = value_to_double(node->get("imaginary_value"), 0.0);
        if (imaginary == 0.0) {
            continue;
        }

        auto imaginary_node = std::make_shared<DBNode>();
        imaginary_node->set(imaginary, "central_value");
        for (const auto& metadata_key : {
                 BlockName("scale"),
                 BlockName("bin_low"),
                 BlockName("bin_high")
             }) {
            if (node->contains(metadata_key)) {
                imaginary_node->set(node->get(metadata_key), metadata_key);
            }
        }

        imaginary_group[key] = imaginary_node;
        has_imaginary_entries = true;
    }

    if (has_imaginary_entries) {
        const auto scale = group.find(BlockName("scale"));
        if (scale != group.end()) {
            imaginary_group[BlockName("scale")] = scale->second;
        }
    }

    return imaginary_group;
}

static void write_flife_like(std::ostream& os,
                             const Prototype& proto,
                             const std::map<BlockName, DBNode::Value>& group)
{
    write_block_header(os, proto, 0.0, false);
    os << "# PDG_code  lifetime        particle\n";

    for (const auto& kv : group) {
        if (kv.first == BlockName("scale")) continue;

        auto node = value_to_node_ptr(kv.second);
        if (!node || !node->contains("central_value")) continue;

        const double val = value_to_double(node->get("central_value"), 0.0);

        auto toks = split_tokens(blockname_primary(kv.first));
        if (toks.empty()) continue;

        os << std::setw(6) << toks[0] << "  ";
        write_value(os, val);
        os << "\n";
    }
    os << "\n";
}

void LhaParser::writeToFile(const std::string &filename,
                            const std::shared_ptr<DBNode> &root) const
{
    if (!root) { LOG_ERROR("LhaParser", "writeToFile: root is null"); return; }

    std::unordered_set<Prototype> allowed = LHA_BLOCKS;
    if (filename.size() >= 5 && filename.substr(filename.size() - 5) == ".flha") {
        allowed.merge(std::unordered_set<Prototype>(FLHA_BLOCKS));
    } else if (filename.size() >= 5 && filename.substr(filename.size() - 5) == ".slha") {
        allowed.merge(std::unordered_set<Prototype>(SLHA_BLOCKS));
    }

    std::ofstream out(filename);
    if (!out) { LOG_ERROR("LhaParser", "Cannot open output file:", filename); return; }

    double Qew = 0.0;
    bool hasQew = read_ew_scale_Q(root, Qew);

    bool usedExternal = false;
    std::map<BlockName, DBNode::Value> fwcoef_from_external;
    {
        const auto prefixes = GroupMapper().get_str();
        fwcoef_from_external = build_fwcoef_group_from_ewscale_blocks(root, prefixes, usedExternal);
    }

    std::vector<BlockName> written_blocks;

    for (const auto& proto : allowed) {
        const auto bn = proto.blockName;

        bool present = false;
        for (const auto& k : root->get_keys()) {
            if (BlockName(k) == bn) { present = true; break; }
        }
        if (!present) continue;

        written_blocks.push_back(bn);

        auto group = root->getGroup({bn});

        const std::string name = blockname_primary(bn);

        if (name == "FWCOEF" || name == "IMFWCOEF") {
            if (usedExternal) write_fwcoef_like(out, proto, fwcoef_from_external, Qew, hasQew);
            else              write_fwcoef_like(out, proto, group,               Qew, hasQew);
        } else if (name == "FOBS" || name == "FOBSSM" || name == "FOBSERR") {
            write_fobs_like(out, proto, group);
        } else if (name == "FMASS") {
            write_fmass_like(out, proto, group);
        } else if (name == "FCONST") {
            write_fconst_like(out, proto, group);
        } else if (name == "FCONSTRATIO") {
            write_fconstratio_like(out, proto, group);
        } else if (name == "FBAG") {
            write_fbag_like(out, proto, group);
        } else if (name == "FLIFE") {
            write_flife_like(out, proto, group);
        } else {
            write_generic_block(out, proto, group);
        }
    }

    // Preserve custom/runtime blocks as generic LHA blocks instead of silently
    // dropping them because no built-in parsing prototype exists. Their entry
    // ids are already serialized as underscore-separated DBNode keys by
    // ParamBlockWriter, which write_generic_block expands into LHA columns.
    auto root_keys = root->get_keys();
    std::sort(root_keys.begin(), root_keys.end(), [](const BlockName& lhs, const BlockName& rhs) {
        return lhs.canonical() < rhs.canonical();
    });

    for (const auto& block_name : root_keys) {
        const bool already_written = std::any_of(
            written_blocks.begin(),
            written_blocks.end(),
            [&](const BlockName& written) { return written == block_name; }
        );
        if (already_written) {
            continue;
        }

        auto group = root->getGroup({block_name});
        Prototype generic{block_name};
        generic.globalScale = group.contains(BlockName("scale"));
        write_generic_block(out, generic, group);
    }

    // LHA-family formats represent complex entries through companion IM...
    // blocks. Generate those blocks when a parameter is stored directly as a
    // complex scalar and no explicit companion block is already present.
    for (const auto& block_name : root_keys) {
        const std::string name = blockname_primary(block_name);
        if (starts_with(name, "IM")) {
            continue;
        }

        const BlockName imaginary_name("IM" + name);
        const bool explicit_companion = std::any_of(
            root_keys.begin(),
            root_keys.end(),
            [&](const BlockName& candidate) { return candidate == imaginary_name; }
        );
        if (explicit_companion) {
            continue;
        }

        const auto imaginary_group = build_imaginary_group(root->getGroup({block_name}));
        if (imaginary_group.empty()) {
            continue;
        }

        Prototype imaginary_proto{imaginary_name};
        if (const Prototype* known = find_proto(allowed, imaginary_name)) {
            imaginary_proto = *known;
        } else {
            imaginary_proto.globalScale =
                imaginary_group.contains(BlockName("scale"));
        }

        if (blockname_primary(imaginary_name) == "IMFWCOEF") {
            write_fwcoef_like(out, imaginary_proto, imaginary_group, Qew, hasQew);
        } else {
            write_generic_block(out, imaginary_proto, imaginary_group);
        }
    }
}
