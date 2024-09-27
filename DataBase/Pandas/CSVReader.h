#ifndef CSVREADER_H
#define CSVREADER_H

#include "DataFrame.h"
#include <string>
#include <vector>
#include <typeindex>
#include <unordered_map>
#include "CSVOptions.h"


class CSVReader {
public:
    DataFrame read_csv(const std::string& filename, const CSVOptions& options = CSVOptions());
};

#endif // CSVREADER_H
