#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <string>
#include <iostream>
#include "Logger.h"

typedef std::pair<std::string, int> ParamId;

enum class ParameterMode {
    FIXED,
    SHIFTABLE
};

class Parameter {

private:

    ParamId id;
    double expected;
    double deviation;
    double value;
    ParameterMode mode; 

public:
    
    inline Parameter() : id({"NullBlock", 0}), expected(0), deviation(0), mode(ParameterMode::FIXED) {}
    Parameter(std::string block, int pdgCode, double mean, double std);

    void set_mode(ParameterMode mode);
    double get_val() const;
    double get_std() const;
    ParamId get_id() const;
    void shift(double shift);

    Parameter& operator=(const Parameter& other) {
        this->id = other.id;
        this->expected = other.expected;
        this->deviation = other.deviation;
        this->mode = other.mode;
        this->value = other.value;
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const Parameter& p) {
        os << "Parameter " << p.id.first << "," << p.id.second << "=" << p.expected << "+-" << p.deviation << std::endl;
        return os;
    }

};


#endif // __PARAMETER_H__

