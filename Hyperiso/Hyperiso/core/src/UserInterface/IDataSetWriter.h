#ifndef IDATASET_WRITER_H
#define IDATASET_WRITER_H

#include "ObjectsOutputs.h"

class IDataSetWriter {
public:
    virtual ~IDataSetWriter() = default;
    virtual void write(const DataSet& ds, const OutputSpec& spec) = 0;
};

#endif