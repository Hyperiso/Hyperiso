#ifndef DEFAULTMAPPINGDATABASEADAPTER_H
#define DEFAULTMAPPINGDATABASEADAPTER_H

#include "IMappingDatabaseAdapter.h"
#include "MappingDatabase.h"
#include <utility>

class DefaultMappingDatabaseAdapter : public IMappingDatabaseAdapter {
public:
    explicit DefaultMappingDatabaseAdapter(std::string instanceName,
                                           std::shared_ptr<MappingDatabase> db)
        : name_(std::move(instanceName)), db_(std::move(db)) {}

    std::unordered_map<std::string, InterpretedParam>
    getParams() const override {
        return db_->getParams();
    }

    std::optional<InterpretedParam>
    getParam(const std::string& name) const override {
        auto m = getParams();
        if (auto it = m.find(name); it != m.end()) return it->second;
        return std::nullopt;
    }

    std::string instanceName() const override { return name_; }

private:
    std::string name_;
    std::shared_ptr<MappingDatabase> db_;
};

#endif
