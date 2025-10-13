#ifndef MARTY_WILSON_PROXY_H
#define MARTY_WILSON_PROXY_H

#include "IMartyWilsonProxy.h"
#include "MartyWilsonAdapter.h"

class MartyWilsonProxy : public IMartyWilsonProxy<Interpreter::InterpretedParam> {
public:
    MartyWilsonProxy() {martyAdapter = MartyWilsonAdapter();}

    void calculate(std::string wilson, std::string model, double Q_match, std::string model_path, bool new_params = false) override {martyAdapter.calculate(wilson, model, Q_match, model_path, new_params);}

    std::set<std::string>  get_special_blocks() override {return martyAdapter.get_special_blocks();}
    std::unordered_set<Interpreter::InterpretedParam> get_dependencies(std::string wilson) override {return martyAdapter.get_dependencies(wilson);}

private:
    MartyWilsonAdapter martyAdapter;
};

#endif