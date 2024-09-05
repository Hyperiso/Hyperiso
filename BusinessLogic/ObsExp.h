#ifndef OBSEXP_H
#define OBSEXP_H

class ObsExp {

    double mean;
    double stat_error;
    double sys_error;

public:
    ObsExp(double center, double stat_error, double sys_error) : mean(center), stat_error(stat_error), sys_error(sys_error) {}

};

#endif // OBSEXP_H