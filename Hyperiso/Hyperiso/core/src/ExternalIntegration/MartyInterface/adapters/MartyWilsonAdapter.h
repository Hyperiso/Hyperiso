#ifndef MARTY_WILSON_ADAPTER_H
#define MARTY_WILSON_ADAPTER_H

#include "MartyInterface.h"
#include "IMartyWilsonAdapter.h"

class MartyWilsonAdapter : public IMartyWilsonAdapter<Interpreter::InterpretedParam> {
public:
    MartyWilsonAdapter() {martyInterface = MartyInterface();}

    void calculate(std::string wilson, std::string model, double Q_match, std::string model_path, bool new_params = false) override {martyInterface.calculate(wilson, model, Q_match, model_path, new_params);}

    std::set<std::string>  get_special_blocks() override {return martyInterface.get_special_blocks();}
    std::unordered_set<Interpreter::InterpretedParam> get_dependencies(std::string wilson) override {return martyInterface.get_dependencies(wilson);}

private:
    MartyInterface martyInterface;
};

#endif