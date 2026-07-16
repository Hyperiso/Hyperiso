#ifndef IFFCALCULATOR_H
#define IFFCALCULATOR_H

template<typename FFType>
class IFFCalculator {
public:
    virtual ~IFFCalculator() = default;	
    virtual double get(FFType, double) = 0;
};

#endif // __IFFCALCULATOR_H__
