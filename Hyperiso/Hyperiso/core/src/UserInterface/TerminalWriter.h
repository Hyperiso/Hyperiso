#ifndef TERMINAL_WRITER_H
#define TERMINAL_WRITER_H

#include "JsonWriter.h"


class TerminalWriter : public IDataSetWriter {
public:
    void write(const DataSet& ds, const OutputSpec& spec) override {
        auto root = dataset_to_dbnode(ds, spec);
        root->printJSONToStream(std::cout, 0);
        std::cout << "\n";
    }
};

#endif