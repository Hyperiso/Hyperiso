#ifndef IMARTY_WILSON_ADAPTER_H
#define IMARTY_WILSON_ADAPTER_H

template<typename T>
class IMartyWilsonAdapter {
public:
    virtual void calculate(std::string wilson, std::string model, double Q_match, std::string model_path, bool new_params = false) = 0;

    virtual std::set<std::string>  get_special_blocks() = 0;
    virtual std::unordered_set<T> get_dependencies(std::string wilson) = 0;


};

#endif