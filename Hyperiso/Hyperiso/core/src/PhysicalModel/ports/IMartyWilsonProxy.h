#ifndef IMARTY_WILSON_PROXY_H
#define IMARTY_WILSON_PROXY_H

#include <set>
#include <unordered_set>
#include <string>

template <typename T>
class IMartyWilsonProxy {
public:

    virtual void calculate(std::string wilson, std::string model, double Q_match, std::string model_path, bool new_params = false) = 0;
    virtual std::set<std::string>  get_special_blocks() = 0;
    virtual std::unordered_set<T> get_dependencies(std::string wilson) = 0;

};

#endif