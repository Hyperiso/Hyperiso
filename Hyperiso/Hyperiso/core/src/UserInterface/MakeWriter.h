#include <memory>
#include "CsvWriter.h"
#include "TerminalWriter.h"
#include "JsonWriter.h"

inline std::unique_ptr<IDataSetWriter> make_writer(OutputFormat fmt,
                                                   const std::string& path_if_file)
{
    switch (fmt) {
        case OutputFormat::CSV:  return std::make_unique<CsvWriter>(path_if_file);
        case OutputFormat::JSON: return std::make_unique<JsonWriter>(path_if_file);
        case OutputFormat::TERMINAL: return std::make_unique<TerminalWriter>();
    }
    throw std::runtime_error("Unknown OutputFormat");
}