
#ifndef CSV_WRITER
#define CSV_WRITER

#include <unordered_set>
#include <fstream>
#include <stdexcept>
#include <algorithm>

#include "ObjectsOutputs.h"
#include "IDataSetWriter.h"
#include "UtilOutput.h"

class CsvWriter : public IDataSetWriter {
public:
    explicit CsvWriter(std::string path) : path_(std::move(path)) {}

    void write(const DataSet& ds, const OutputSpec& spec) override;

private:
    std::string path_;

    static std::string double_to_string(double x);

    static void write_row(std::ostream& os,
                          const std::vector<std::string>& cells,
                          char sep);

    static std::string escape_csv(std::string s, char sep);
};

#endif