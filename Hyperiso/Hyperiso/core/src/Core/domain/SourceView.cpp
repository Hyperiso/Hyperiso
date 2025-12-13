#include "SourcesView.h"

std::string join(const std::vector<std::string>& v) {
    std::ostringstream oss;
    for (size_t i = 0; i < v.size(); ++i) { if (i) oss << ", "; oss << v[i]; }
    return oss.str();
}

std::string to_string(const LhaID& id) {
    return id.to_string();
}

std::string to_string(const ParamId& id) {
    std::ostringstream oss;
    oss << ParameterTypeMapper::str(id.type.value()) << "::" << id.block << "::" << id.code.to_string();
    return oss.str();
}

BlockSrc::BlockSrc(const std::unordered_map<std::string, std::shared_ptr<Block>>& m,
                    std::string ctx)
    : m_(&m), ctx_(std::move(ctx)) {}

bool BlockSrc::has_block(std::string_view name) const {
    return m_->find(std::string(name)) != m_->end();
}

const std::shared_ptr<Block>& BlockSrc::block(std::string_view name) const {
    auto it = m_->find(std::string(name));
    if (it == m_->end()) {
        std::vector<std::string> avail;
        avail.reserve(m_->size());
        for (auto& [k,_] : *m_) avail.push_back(k);
        std::ostringstream oss;
        oss << "Missing source Block '" << name << "'";
        if (!ctx_.empty()) oss << " in " << ctx_;
        oss << ". Available blocks: [" << join(avail) << "].";
        throw std::invalid_argument(oss.str());
    }
    return it->second;
}

std::shared_ptr<Parameter> BlockSrc::get_param(std::string_view blk, const LhaID& code) const {
    // auto& b = block(blk);
    // if (!b->contains(code)) {
    //     std::cout << "here ?" << std::endl;
    //     std::vector<std::string> ids;
    //     for (auto& id : b->getAllIDs()) ids.push_back(to_string(id));
    //     std::ostringstream oss;
    //     oss << "Missing parameter " << to_string(code)
    //         << " in Block '" << blk << "'";
    //     if (!ctx_.empty()) oss << " (context: " << ctx_ << ")";
    //     oss << ". Available IDs: [" << join(ids) << "].";
    //     throw std::invalid_argument(oss.str());
    // }
    // return b->retrieve(code);
    auto& b = block(blk);
    return b->retrieve(code);
}

std::shared_ptr<Parameter> BlockSrc::get_param(std::string_view blk, int code) const {
    return get_param(blk, LhaID(code));
}
std::shared_ptr<Parameter> BlockSrc::get_param(std::string_view blk, std::pair<int,int> code) const {
    return get_param(blk, LhaID(code.first, code.second));
}

scalar_t BlockSrc::get_val(std::string_view blk, std::initializer_list<long> code) const {
    return get_param(blk, code)->get_val();
}

double BlockSrc::get_val(std::string_view blk, const LhaID& code) const {
    return get_param(blk, code)->get_val();
}

double BlockSrc::get_val(std::string_view blk, int code) const {
    return get_param(blk, LhaID(code))->get_val();
}

double BlockSrc::get_val(std::string_view blk, std::pair<int,int> code) const {
    return get_param(blk, LhaID(code.first, code.second))->get_val();
}

const std::unordered_map<std::string, std::shared_ptr<Block>>& BlockSrc::raw() const { return *m_; }



ParamSrc::ParamSrc(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& m,
                    std::string ctx)
    : m_(&m), ctx_(std::move(ctx)) {}

bool ParamSrc::has(const ParamId& id) const { return m_->find(id) != m_->end(); }

const std::shared_ptr<Parameter>& ParamSrc::require(const ParamId& id) const {
    auto it = m_->find(id);
    if (it == m_->end()) {
        std::vector<std::string> avail;
        avail.reserve(m_->size());
        for (auto& [k,_] : *m_) avail.push_back(to_string(k));
        std::ostringstream oss;
        oss << "Missing source Parameter '" << to_string(id) << "'";
        if (!ctx_.empty()) oss << " in " << ctx_;
        oss << ". Available: [" << join(avail) << "].";
        throw std::invalid_argument(oss.str());
    }
    return it->second;
}

std::shared_ptr<Parameter> ParamSrc::get_param(const ParamId& id) const { return require(id); }
scalar_t ParamSrc::get_val(const ParamId& id) const { return require(id)->get_val(); }

scalar_t ParamSrc::get_val(ParameterType t, std::string_view block, std::initializer_list<long> code) const {
    return get_val(ParamId{t, std::string(block), code});
}
scalar_t ParamSrc::get_val(ParameterType t, std::string_view block, const LhaID& code) const {
    return get_val(ParamId{t, std::string(block), code});
}
scalar_t ParamSrc::get_val(ParameterType t, std::string_view block, int code) const {
    return get_val(ParamId{t, std::string(block), LhaID(code)});
}
scalar_t ParamSrc::get_val(ParameterType t, std::string_view block, std::pair<int,int> code) const {
    return get_val(ParamId{t, std::string(block), LhaID(code.first, code.second)});
}

const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& ParamSrc::raw() const { return *m_; }

BlockSrc as_block_src(const std::unordered_map<std::string, std::shared_ptr<Block>>& m, std::string ctx) {
    return BlockSrc(m, std::move(ctx));
}
ParamSrc as_param_src(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& m, std::string ctx) {
    return ParamSrc(m, std::move(ctx));
}
