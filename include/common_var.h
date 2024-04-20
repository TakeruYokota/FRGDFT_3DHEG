//
// Created by Takeru Yokota on 2019/12/25.
//

#ifndef COULOMB3D_2_COMMON_VAR_H
#define COULOMB3D_2_COMMON_VAR_H

#include <cmath>

const int Ns{2};

struct param
{
    double u;
    double p;
    double pz1;
};

struct param_wo_u
{
    double p;
    double pz1;
};

struct param_pf
{
    double u;
    double p;
    double pf;
};

inline double xlogy(double x, double y)
{
    return y == 0. ? 0. : x * log(y);
}

inline long double xlogyl(long double x, long double y)
{
    return y == 0. ? 0. : x * logl(y);
}

#endif // COULOMB3D_2_COMMON_VAR_H
