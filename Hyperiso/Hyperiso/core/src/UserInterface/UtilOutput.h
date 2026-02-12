#include <string>
#include <variant>
#include <iomanip>

#include "ObjectsOutputs.h"

inline std::string to_string_value(const Value& v) {
    struct {
        std::string operator()(std::monostate) const { return ""; }
        std::string operator()(double x) const {
            std::ostringstream oss;
            oss << std::setprecision(17) << x;
            return oss.str();
        }
        std::string operator()(int64_t x) const { return std::to_string(x); }
        std::string operator()(bool b) const { return b ? "true" : "false"; }
        std::string operator()(const std::string& s) const { return s; }
    } visitor;
    return std::visit(visitor, v);
}