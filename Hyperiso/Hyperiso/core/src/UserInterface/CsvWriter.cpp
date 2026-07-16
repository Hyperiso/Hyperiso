#include "CsvWriter.h"

void CsvWriter::write(const DataSet& ds, const OutputSpec& spec) {
    std::ofstream out(path_);
    if (!out) throw std::runtime_error("Cannot open CSV file: " + path_);

    const char sep = spec.csv_sep;

    std::vector<std::string> header;
    header.reserve(ds.vars.size() + ds.outputs_schema.size());
    for (const auto& v : ds.vars) header.push_back(v.name);

    std::vector<std::string> outCols = ds.outputs_schema;
    if (outCols.empty()) {
        std::unordered_set<std::string> uniq;
        for (const auto& p : ds.points)
            for (const auto& kv : p.y) uniq.insert(kv.first);

        outCols.assign(uniq.begin(), uniq.end());
        std::sort(outCols.begin(), outCols.end());
    }

    if (!spec.keys.empty()) {
        std::unordered_set<std::string> wanted(spec.keys.begin(), spec.keys.end());
        std::vector<std::string> filtered;
        filtered.reserve(outCols.size());
        for (auto& k : outCols) if (wanted.count(k)) filtered.push_back(k);
        outCols = std::move(filtered);
    }

    header.insert(header.end(), outCols.begin(), outCols.end());

    if (spec.csv_write_header) {
        write_row(out, header, sep);
    }

    for (const auto& p : ds.points) {
        std::vector<std::string> row;
        row.reserve(header.size());

        for (size_t i = 0; i < ds.vars.size(); ++i) {
            double val = (i < p.x.size()) ? p.x[i] : 0.0;
            row.push_back(double_to_string(val));
        }

        for (const auto& col : outCols) {
            auto it = p.y.find(col);
            if (it == p.y.end()) row.push_back("");
            else row.push_back(escape_csv(to_string_value(it->second), sep));
        }

        write_row(out, row, sep);
    }
}

std::string CsvWriter::double_to_string(double x) {
    std::ostringstream oss;
    oss << std::setprecision(17) << x;
    return oss.str();
}

void CsvWriter::write_row(std::ostream& os,
                        const std::vector<std::string>& cells,
                        char sep)
{
    for (size_t i = 0; i < cells.size(); ++i) {
        if (i) os << sep;
        os << cells[i];
    }
    os << "\n";
}

std::string CsvWriter::escape_csv(std::string s, char sep) {
    // RFC4180 (simplified)
    bool need = false;
    for (char c : s) {
        if (c == sep || c == '"' || c == '\n' || c == '\r') { need = true; break; }
    }
    if (!need) return s;

    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('"');
    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out.push_back(c);
    }
    out.push_back('"');
    return out;
}