#ifndef JSON_WRITER_H
#define JSON_WRITER_H

#include "DBNode.h"
#include "IDataSetWriter.h"

DBNode::Value to_dbnode_value(const Value& v);

std::shared_ptr<DBNode> make_object_node();

std::shared_ptr<DBNode> make_scalar_node(const Value& v);

DBNode::Value make_array_value(const std::vector<std::shared_ptr<DBNode>>& elems);

std::shared_ptr<DBNode> dataset_to_dbnode(const DataSet& ds, const OutputSpec& spec);

class JsonWriter : public IDataSetWriter {
public:
    explicit JsonWriter(std::string path) : path_(std::move(path)) {}

    void write(const DataSet& ds, const OutputSpec& spec) override;

private:
    std::string path_;
};


#endif