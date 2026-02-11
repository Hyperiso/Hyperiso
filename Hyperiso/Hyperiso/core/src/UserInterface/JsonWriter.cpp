#include "JsonWriter.h"

DBNode::Value to_dbnode_value(const Value& v) {
    if (std::holds_alternative<std::monostate>(v)) {
        return std::shared_ptr<DBNode>(nullptr);
    }
    if (auto p = std::get_if<double>(&v)) return *p;
    if (auto p = std::get_if<int64_t>(&v)) return static_cast<int>(*p);
    if (auto p = std::get_if<bool>(&v)) return *p;
    if (auto p = std::get_if<std::string>(&v)) return BlockName(*p);

    return std::shared_ptr<DBNode>(nullptr);
}

std::shared_ptr<DBNode> make_object_node() {
    return std::make_shared<DBNode>();
}

std::shared_ptr<DBNode> make_scalar_node(const Value& v) {
    auto n = std::make_shared<DBNode>();
    n->set(to_dbnode_value(v), "");
    return n;
}

DBNode::Value make_array_value(const std::vector<std::shared_ptr<DBNode>>& elems) {
    return elems; 
}

std::shared_ptr<DBNode> dataset_to_dbnode(const DataSet& ds, const OutputSpec& spec) {
    auto root = make_object_node();

    {
        auto meta = make_object_node();
        for (const auto& [k,v] : ds.meta) {
            meta->set(to_dbnode_value(v), k);
        }
        root->set(meta, "meta");
    }

    {
        std::vector<std::shared_ptr<DBNode>> varsArr;
        varsArr.reserve(ds.vars.size());

        for (const auto& v : ds.vars) {
            auto o = make_object_node();
            o->set(BlockName(v.name), "name");
            o->set(BlockName(v.block), "block");
            o->set(static_cast<int>(v.pdg), "pdg");
            o->set(v.min, "min");
            o->set(v.max, "max");
            o->set(v.step, "step");
            varsArr.push_back(o);
        }

        auto scan = make_object_node();
        scan->set(make_array_value(varsArr), "vars");
        root->set(scan, "scan");
    }

    if (!ds.outputs_schema.empty()) {
        std::vector<std::shared_ptr<DBNode>> schemaArr;
        schemaArr.reserve(ds.outputs_schema.size());
        for (const auto& k : ds.outputs_schema) schemaArr.push_back(make_scalar_node(std::string(k)));
        root->set(make_array_value(schemaArr), "outputs_schema");
    }

    {
        std::unordered_set<std::string> wanted;
        const bool filter = !spec.keys.empty();
        if (filter) wanted = {spec.keys.begin(), spec.keys.end()};

        std::vector<std::shared_ptr<DBNode>> pointsArr;
        pointsArr.reserve(ds.points.size());

        for (const auto& p : ds.points) {
            auto jp = make_object_node();

            auto xobj = make_object_node();
            for (size_t i = 0; i < ds.vars.size(); ++i) {
                const auto& var = ds.vars[i];
                const double val = (i < p.x.size()) ? p.x[i] : 0.0;
                xobj->set(val, var.name);
            }
            jp->set(xobj, "x");

            auto yobj = make_object_node();
            for (const auto& [k,v] : p.y) {
                if (filter && !wanted.count(k)) continue;
                yobj->set(to_dbnode_value(v), k);
            }
            jp->set(yobj, "y");

            pointsArr.push_back(jp);
        }

        root->set(make_array_value(pointsArr), "points");
    }

    return root;
}


void JsonWriter::write(const DataSet& ds, const OutputSpec& spec) {
    auto root = dataset_to_dbnode(ds, spec);

    std::ofstream out(path_);
    if (!out) throw std::runtime_error("Cannot open JSON file: " + path_);

    root->printJSONToStream(out, 0);
    out << "\n";
}


