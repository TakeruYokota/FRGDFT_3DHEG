#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <time.h>
#include <cmath>
#include "integ_prep.h"
#include "common_var.h"
#include "calc_c.h"

double beta{1.};
double lam{1.};

double log1x_x(double x)
{
    double thr{1e-5};
    if (fabs(x) < thr)
    {
        return x * x * (-0.5 + x * (1. / 3. + x * (-0.25 + x * (0.2 - x / 6.))));
    }
    else
    {
        return log(1. + x) - x;
    }
}

double logfunc1(double x)
{
    double x1{x + 1.};
    double thr{1e-5};
    if (fabs(x) < thr)
    {
        return gsl_pow_3(x) * (-1. / 3. + x * (1. / 2. + x * (-3. / 5 + 2. * x / 3.)));
    }
    else
    {
        return 2. * log(x1) + 1. / x1 - x1;
    }
}

double calc_G20(double u, double p)
{
    double thr{1e-7};
    double factor{(double)Ns / (4. * M_PI * M_PI)};
    double u2{u * u};
    double u2_1{1. + u2};
    double p2{p * p};
    double p4{p2 * p2};

    if (p < thr)
    {
        double cor_p2{p2 / (12. * gsl_pow_2(u2_1))};
        double cor_p4{(-1. + 5. * u2) * p4 / (240. * gsl_pow_4(u2_1))};
        return 2. * factor * (1. - u * atan2(1., u) + cor_p2 + cor_p4);
    }
    else
    {
        double tp{1. + 0.5 * p};
        double tm{1. - 0.5 * p};
        if (tp / u < 1e-3)
        {
            double v{1. / u};
            double v2{v * v};
            double v4{v2 * v2};
            double v6{v4 * v2};

            return factor * (2. * v2 / 3. + (-12. - 5. * p2) * v4 / 30. + (48. + 56. * p2 + 7. * p4) * v6 / 168.);
        }
        else
        {
            return factor * (1. - u * atan2(tp, u) - u * atan2(tm, u) - (u2 + tp * tm) * log((u2 + tm * tm) / (u2 + tp * tp)) / (2. * p));
        }
    }
}

double tanh_x(double x)
{
    double thr{1e-5};
    if (x > thr)
    {
        return tanh(x) / x;
    }
    else
    {
        return 1. + x * x * (-1. + 0.4 * x * x) / 3.;
    }
}

double tanh_x_1(double x)
{
    double thr{1e-5};

    if (x > thr)
    {
        return tanh(x) / x - 1.;
    }
    else
    {
        return x * x * (-1 / 3. + 2. * x * x / 15.);
    }
}

// double calc_Erpa(double rs)
// {
//     double pf{pow(2.25 * M_PI, 1. / 3.) / rs};

//     double result{lam * lam * (1. - log(2.)) / (M_PI * M_PI) * log(4. * lam / (beta * beta * M_PI * pf))};
//     result += -0.03266559211908057 * lam * lam;

//     double epsabs = 1e-8;
//     double epsrel = 1e-8;
//     size_t limit = 100;
//     double result1, result2;
//     double abserr1, abserr2;

//     IntegrationWorkspace wsp1(limit);
//     IntegrationWorkspace wsp2(limit);
//     IntegrationWorkspace wsp3(limit);
//     IntegrationWorkspace wsp4(limit);
//     IntegrationWorkspace wsp5(limit);

//     auto integ1 = make_gsl_function([&](double u) {
//         double f;
//         if (fabs(u) > 1e4)
//         {
//             double ui2{1. / (u * u)};
//             f = ui2 * (1. / 3. - ui2 * (0.2 - ui2 / 7.));
//         }
//         else if (fabs(u) < 1e-4)
//         {
//             double u2{u * u};
//             f = 1. - 0.5 * M_PI * u + u2 * (1. - u2 * (1. / 3. - u2 / 5.));
//         }
//         else
//         {
//             f = 1. - u * atan2(1., u);
//         }

//         double x{4. * lam * f / (M_PI * beta * beta * pf)};

//         if (x < 1e-4)
//         {
//             return x * x * x * (2. - x * (0.75 - x * (0.4 - x * (0.25 - x * (6. / 35. - 0.125 * x))))) / 3.;
//         }
//         else
//         {
//             return x - 0.5 * x * x - (1. - x * x) * log(1. + x);
//         }
//     });

//     gsl_integration_qagiu(integ1, 0., epsabs, epsrel, limit, wsp1, &result1, &abserr1);
//     result += -(3. * gsl_pow_4(beta) * pf * pf) / (16. * M_PI) * result1;

//     if (pf > 1.)
//     {
//         auto integ2 = make_gsl_function([&](double u) {
//             auto integ2_1 = make_gsl_function([&](double p) {
//                 if (p == 0.)
//                 {
//                     return 0.;
//                 }
//                 double p2{p * p};
//                 double g0{calc_G20(u, p)};
//                 double g0_p0{calc_G20(u, 0)};
//                 double fg0{4. * M_PI * lam * g0 / (pf * p2)};
//                 double fg0_p0{4. * M_PI * lam * g0_p0 / (pf * p2)};
//                 return p2 * p * (log1x_x(fg0) - log1x_x(fg0_p0));
//             });
//             gsl_integration_qag(integ2_1, 0., beta, epsabs, epsrel, limit, GSL_INTEG_GAUSS31, wsp2, &result2, &abserr2);
//             return result2;
//         });

//         gsl_integration_qagiu(integ2, 0., epsabs, epsrel, limit, wsp3, &result1, &abserr1);
//         result += (3. * pf * pf) / (4. * M_PI) * result1;

//         auto integ3 = make_gsl_function([&](double u) {
//             auto integ3_1 = make_gsl_function([&](double p) {
//                 if (p == 0.)
//                 {
//                     return 0.;
//                 }
//                 double p2{p * p};
//                 double g0{calc_G20(u, p)};
//                 double fg0{4. * M_PI * lam * g0 / (pf * p2)};
//                 return p2 * p * log1x_x(fg0);
//             });
//             gsl_integration_qagiu(integ3_1, beta, epsabs, epsrel, limit, wsp4, &result2, &abserr2);
//             return result2;
//         });

//         gsl_integration_qagiu(integ3, 0., epsabs, epsrel, limit, wsp5, &result1, &abserr1);
//         result += (3. * pf * pf) / (4. * M_PI) * result1;
//     }
//     else
//     {
//         auto integ2 = make_gsl_function([&](double u) {
//             auto integ2_1 = make_gsl_function([&](double p) {
//                 if (p == 0.)
//                 {
//                     return 0.;
//                 }
//                 double p2{p * p};
//                 double g0{pf * calc_G20(u / pf, p / pf)};
//                 double g0_p0{pf * calc_G20(u / pf, 0)};
//                 double fg0{4. * M_PI * lam * g0 / p2};
//                 double fg0_p0{4. * M_PI * lam * g0_p0 / p2};
//                 return p2 * p * (log1x_x(fg0) - log1x_x(fg0_p0));
//             });
//             gsl_integration_qag(integ2_1, 0., pf * beta, epsabs, epsrel, limit, GSL_INTEG_GAUSS31, wsp2, &result2, &abserr2);
//             return result2;
//         });

//         gsl_integration_qagiu(integ2, 0., epsabs, epsrel, limit, wsp3, &result1, &abserr1);
//         result += (3. * pf * pf) / (4. * M_PI) * result1;

//         auto integ3 = make_gsl_function([&](double u) {
//             auto integ3_1 = make_gsl_function([&](double p) {
//                 if (p == 0.)
//                 {
//                     return 0.;
//                 }
//                 double p2{p * p};
//                 double g0{pf * calc_G20(u / pf, p / pf)};
//                 double fg0{4. * M_PI * lam * g0 / p2};
//                 return p2 * p * log1x_x(fg0);
//             });
//             gsl_integration_qagiu(integ3_1, pf * beta, epsabs, epsrel, limit, wsp4, &result2, &abserr2);
//             return result2;
//         });

//         gsl_integration_qagiu(integ3, 0., epsabs, epsrel, limit, wsp5, &result1, &abserr1);
//         result += 3. / (4. * M_PI * gsl_pow_3(pf)) * result1;
//     }

//     return result;
// }

double calc_R(double u, double fu, double fd)
{
    double thr_u{1e-3};
    double tu{0.};
    double td{0.};
    if (fabs(u) < thr_u * fd)
    {
        double ufd{u / fd};
        td = fd * (1. - ufd * (0.5 * M_PI - ufd * (1. - ufd * ufd / 3.)));
    }
    else
    {
        td = fd - u * atan2(fd, u);
    }

    if (fabs(u) < thr_u * fu)
    {
        double ufu{u / fu};
        tu = fu * (1. - ufu * (0.5 * M_PI - ufu * (1. - ufu * ufu / 3.)));
    }
    else
    {
        tu = fu - u * atan2(fu, u);
    }

    return 0.5 * (tu + td);
}

double calc_Erpa(double rs, double zeta)
{
    double pf{pow(2.25 * M_PI, 1. / 3.) / rs};
    double fu{pow(1. + zeta, 1. / 3.)};
    double fd{pow(1. - zeta, 1. / 3.)};
    double thr_fd_A{1.e-3};

    double A{0.};
    if (fd > thr_fd_A)
    {
        A = 0.5 / (M_PI * M_PI) * (2. * (1. - log(2.)) + fu * fd * (fu + fd) - fu * fu * fu * log((fu + fd) / fu) - fd * fd * fd * log((fu + fd) / fd));
    }
    else
    {
        A = 0.5 / (M_PI * M_PI) * (2. * (1. - log(2.)) + fu * fd * (fu + fd) - fu * fu * fu * log((fu + fd) / fu));
    }

    double result{0.5 * lam * lam * A * (log(4. * lam / (beta * beta * M_PI * pf)) - 0.5)};

    double epsabs = 1e-8;
    double epsrel = 1e-8;
    size_t limit = 100;
    double result1, result2;
    double abserr1, abserr2;

    IntegrationWorkspace wsp1(limit);
    IntegrationWorkspace wsp2(limit);
    IntegrationWorkspace wsp3(limit);
    IntegrationWorkspace wsp4(limit);
    IntegrationWorkspace wsp5(limit);

    auto integ0 = make_gsl_function([&](double u) {
        double R{calc_R(u, fu, fd)};
        return R * R * log(R);
    });

    gsl_integration_qagiu(integ0, 0., epsabs, epsrel, limit, wsp1, &result1, &abserr1);
    result += lam * lam * 3. * result1 / (M_PI * M_PI * M_PI);

    auto integ1 = make_gsl_function([&](double u) {
        double x{4. * lam * calc_R(u, fu, fd) / (M_PI * beta * beta * pf)};

        if (x < 1e-4)
        {
            return x * x * x * (2. - x * (0.75 - x * (0.4 - x * (0.25 - x * (6. / 35. - 0.125 * x))))) / 3.;
        }
        else
        {
            return x - 0.5 * x * x - (1. - x * x) * log(1. + x);
        }
    });

    gsl_integration_qagiu(integ1, 0., epsabs, epsrel, limit, wsp1, &result1, &abserr1);
    result += -(3. * gsl_pow_4(beta) * pf * pf) / (16. * M_PI) * result1;

    if (pf > 1.)
    {
        auto integ2 = make_gsl_function([&](double u) {
            auto integ2_1 = make_gsl_function([&](double p) {
                if (p == 0.)
                {
                    return 0.;
                }
                double p2{p * p};
                // double g0{calc_G20(u, p)};
                // double g0_p0{calc_G20(u, 0)};

                double g0{0.};
                double g0_p0{0.};

                if (zeta == 0.)
                {
                    g0 = calc_G20(u, p);
                    g0_p0 = calc_G20(u, 0);
                }
                else
                {
                    double thr_fd{0.01};
                    double g0u{calc_G20(u / fu, p / fu)};
                    double g0d{fd > thr_fd ? calc_G20(u / fd, p / fd) : 0.};
                    double g0u_p0{calc_G20(u / fu, 0.)};
                    double g0d_p0{fd > thr_fd ? calc_G20(u / fd, 0.) : 0.};

                    g0 = 0.5 * (fu * g0u + fd * g0d);
                    g0_p0 = 0.5 * (fu * g0u_p0 + fd * g0d_p0);
                }

                double fg0{4. * M_PI * lam * g0 / (pf * p2)};
                double fg0_p0{4. * M_PI * lam * g0_p0 / (pf * p2)};
                return p2 * p * (log1x_x(fg0) - log1x_x(fg0_p0));
            });
            gsl_integration_qag(integ2_1, 0., beta, epsabs, epsrel, limit, GSL_INTEG_GAUSS31, wsp2, &result2, &abserr2);
            return result2;
        });

        gsl_integration_qagiu(integ2, 0., epsabs, epsrel, limit, wsp3, &result1, &abserr1);
        result += (3. * pf * pf) / (4. * M_PI) * result1;

        auto integ3 = make_gsl_function([&](double u) {
            auto integ3_1 = make_gsl_function([&](double p) {
                if (p == 0.)
                {
                    return 0.;
                }
                double p2{p * p};

                double g0{0.};

                if (zeta == 0.)
                {
                    g0 = calc_G20(u, p);
                }
                else
                {
                    double thr_fd{0.01};
                    double g0u{calc_G20(u / fu, p / fu)};
                    double g0d{fd > thr_fd ? calc_G20(u / fd, p / fd) : 0.};
                    g0 = 0.5 * (fu * g0u + fd * g0d);
                }

                double fg0{4. * M_PI * lam * g0 / (pf * p2)};
                return p2 * p * log1x_x(fg0);
            });
            gsl_integration_qagiu(integ3_1, beta, epsabs, epsrel, limit, wsp4, &result2, &abserr2);
            return result2;
        });

        gsl_integration_qagiu(integ3, 0., epsabs, epsrel, limit, wsp5, &result1, &abserr1);
        result += (3. * pf * pf) / (4. * M_PI) * result1;
    }
    else
    {
        auto integ2 = make_gsl_function([&](double u) {
            auto integ2_1 = make_gsl_function([&](double p) {
                if (p == 0.)
                {
                    return 0.;
                }
                double p2{p * p};

                // double g0{pf * calc_G20(u / pf, p / pf)};
                // double g0_p0{pf * calc_G20(u / pf, 0)};

                double g0{0.};
                double g0_p0{0.};

                if (zeta == 0.)
                {
                    g0 = calc_G20(u / pf, p / pf);
                    g0_p0 = calc_G20(u / pf, 0);
                }
                else
                {
                    double thr_fd{0.01};
                    double g0u{calc_G20(u / (fu * pf), p / (fu * pf))};
                    double g0d{fd > thr_fd ? calc_G20(u / (fd * pf), p / (fd * pf)) : 0.};
                    double g0u_p0{calc_G20(u / (fu * pf), 0.)};
                    double g0d_p0{fd > thr_fd ? calc_G20(u / (fd * pf), 0.) : 0.};

                    g0 = 0.5 * (fu * g0u + fd * g0d);
                    g0_p0 = 0.5 * (fu * g0u_p0 + fd * g0d_p0);
                }

                g0 *= pf;
                g0_p0 *= pf;

                double fg0{4. * M_PI * lam * g0 / p2};
                double fg0_p0{4. * M_PI * lam * g0_p0 / p2};
                return p2 * p * (log1x_x(fg0) - log1x_x(fg0_p0));
            });
            gsl_integration_qag(integ2_1, 0., pf * beta, epsabs, epsrel, limit, GSL_INTEG_GAUSS31, wsp2, &result2, &abserr2);
            return result2;
        });

        //TODO:check effect of coefficient
        gsl_integration_qagiu(integ2, 0., epsabs, epsrel, limit, wsp3, &result1, &abserr1);
        result += 0.75 / (M_PI * gsl_pow_3(pf)) * result1;
        // result += (3. * pf * pf) / (4. * M_PI) * result1;

        auto integ3 = make_gsl_function([&](double u) {
            auto integ3_1 = make_gsl_function([&](double p) {
                if (p == 0.)
                {
                    return 0.;
                }
                double p2{p * p};
                // double g0{pf * calc_G20(u / pf, p / pf)};

                double g0{0.};
                if (zeta == 0.)
                {
                    g0 = calc_G20(u / pf, p / pf);
                }
                else
                {
                    double thr_fd{0.01};
                    double g0u{calc_G20(u / (fu * pf), p / (fu * pf))};
                    double g0d{fd > thr_fd ? calc_G20(u / (fd * pf), p / (fd * pf)) : 0.};
                    g0 = 0.5 * (fu * g0u + fd * g0d);
                }

                g0 *= pf;

                double fg0{4. * M_PI * lam * g0 / p2};
                return p2 * p * log1x_x(fg0);
            });
            gsl_integration_qagiu(integ3_1, pf * beta, epsabs, epsrel, limit, wsp4, &result2, &abserr2);
            return result2;
        });

        gsl_integration_qagiu(integ3, 0., epsabs, epsrel, limit, wsp5, &result1, &abserr1);
        result += 0.75 / (M_PI * gsl_pow_3(pf)) * result1;
    }

    return result;
}

// double energy_int(double u, double p, double thr_u, double thr_p, double pf)
// {
//     if (p == 0.)
//     {
//         return 0.;
//     }

//     double Ci{calc_c(u, p, thr_u, thr_p)};
//     double G2{calc_G20(u, p)};

//     double p2{p * p};

//     double p_f_c{lam * sqrt(4. * M_PI * Ci) / pf};
//     double f_c{p_f_c / p};
//     double p2_f_g{lam * 4. * M_PI * G2 / pf};
//     double f_g{p2_f_g / p2};

//     double p2_log_1;
//     double p2_log_2;

//     if (fabs(f_c) < 1e-4)
//     {
//         double c2{f_c * f_c};
//         p2_log_1 = p2_f_g * c2 * (-1. / 3. + c2 * ((12. + 7. * f_g) / (90. * (1. + f_g)) + c2 * (-153. - 180. * f_g - 62. * f_g * f_g) / (2835. * gsl_pow_2(1. + f_g)))) / (1. + f_g);
//     }
//     else
//     {
//         p2_log_1 = p2 * log((p2 + p2_f_g * tanh_x(f_c)) / (p2 + p2_f_g));
//     }

//     if (fabs(f_c) < 1e-4)
//     {
//         double c2{f_c * f_c};
//         p2_log_2 = p2 * c2 * (1. / 2. - c2 * (1. / 12. - c2 / 45.));
//     }

//     else if (fabs(f_c) < 10)
//     {
//         p2_log_2 = p2 * log(cosh(f_c));
//     }
//     else
//     {
//         p2_log_2 = p2 * log(0.5 * (1. + exp(-2. * f_c))) + p * p_f_c;
//     }

//     return p2_log_1 + p2_log_2;
// }

double energy_int(double u, double p, double thr_u, double thr_p, double pf, double zeta)
{
    if (p == 0.)
    {
        return 0.;
    }

    double Ci{0.};
    double G2{0.};

    if (zeta == 0.)
    {
        Ci = calc_c(u, p, thr_u, thr_p);
        G2 = calc_G20(u, p);
    }
    else
    {
        double thr_fd{0.01};
        double fu{pow(1. + zeta, 1. / 3.)};
        double fd{pow(1. - zeta, 1. / 3.)};
        double Ciu{calc_c(u / fu, p / fu, thr_u, thr_p)};
        double Cid{fd > thr_fd ? calc_c(u / fd, p / fd, thr_u, thr_p) : 0.};
        double G2u{calc_G20(u / fu, p / fu)};
        double G2d{fd > thr_fd ? calc_G20(u / fd, p / fd) : 0.};
        Ci = 0.5 * (Ciu + Cid);
        G2 = 0.5 * (fu * G2u + fd * G2d);
    }

    double p2{p * p};

    double p_f_c{lam * sqrt(4. * M_PI * Ci) / pf};
    double f_c{p_f_c / p};
    double p2_f_g{lam * 4. * M_PI * G2 / pf};
    double f_g{p2_f_g / p2};

    double p2_log_1;
    double p2_log_2;

    if (fabs(f_c) < 1e-4)
    {
        double c2{f_c * f_c};
        p2_log_1 = p2_f_g * c2 * (-1. / 3. + c2 * ((12. + 7. * f_g) / (90. * (1. + f_g)) + c2 * (-153. - 180. * f_g - 62. * f_g * f_g) / (2835. * gsl_pow_2(1. + f_g)))) / (1. + f_g);
    }
    else
    {
        p2_log_1 = p2 * log((p2 + p2_f_g * tanh_x(f_c)) / (p2 + p2_f_g));
    }

    if (fabs(f_c) < 1e-4)
    {
        double c2{f_c * f_c};
        p2_log_2 = p2 * c2 * (1. / 2. - c2 * (1. / 12. - c2 / 45.));
    }

    else if (fabs(f_c) < 10)
    {
        p2_log_2 = p2 * log(cosh(f_c));
    }
    else
    {
        p2_log_2 = p2 * log(0.5 * (1. + exp(-2. * f_c))) + p * p_f_c;
    }

    return p2_log_1 + p2_log_2;
}

double calc_energy_gb(double rs, double zeta)
{
    double Erpa{calc_Erpa(rs, zeta)};
    return Erpa + lam * lam * 0.0241791589181444;
}

double calc_energy(double rs, double zeta)
{
    double pf{pow(2.25 * M_PI, 1. / 3.) / rs};

    double Erpa{calc_Erpa(rs, zeta)};

    double thr_p{1e-4};
    double thr_u{1e-4};
    double epsabs = 1e-7;
    double epsrel = 1e-6;
    size_t limit = 100;
    double result1, result2;
    double abserr1, abserr2;

    IntegrationWorkspace wsp1(limit);
    IntegrationWorkspace wsp2(limit);

    auto integ = make_gsl_function([&](double p) {
        auto integ_inner = make_gsl_function([&](double u) {
            return energy_int(u, p, thr_u, thr_p, pf, zeta);
        });
        gsl_integration_qagiu(integ_inner, 0, epsabs, epsrel, limit, wsp2, &result2, &abserr2);
        return p * result2;
    });

    gsl_integration_qagiu(integ, 0, epsabs, epsrel, limit, wsp1, &result1, &abserr1);

    double Ecorr{Erpa + result1 * 3. * pf * pf / (4. * M_PI)};
    return Ecorr;
}

///energy diff////////////////////////////////////////////////////////////////////////
double calc_energy_gb_diff(double rs)
{
    double pf{pow(2.25 * M_PI, 1. / 3.) / rs};

    double dedpf{-lam * lam * (1. - log(2.)) / (M_PI * M_PI * pf)};

    double epsabs = 1e-8;
    double epsrel = 1e-8;
    size_t limit = 100;
    double result1, result2;
    double abserr1, abserr2;

    IntegrationWorkspace wsp1(limit);
    IntegrationWorkspace wsp2(limit);
    IntegrationWorkspace wsp3(limit);
    IntegrationWorkspace wsp4(limit);
    IntegrationWorkspace wsp5(limit);

    auto integ1 = make_gsl_function([&](double u) {
        double Ru;
        if (fabs(u) > 1e4)
        {
            double ui2{1. / (u * u)};
            Ru = ui2 * (1. / 3. - ui2 * (0.2 - ui2 / 7.));
        }
        else if (fabs(u) < 1e-4)
        {
            double u2{u * u};
            Ru = 1. - 0.5 * M_PI * u + u2 * (1. - u2 * (1. / 3. - u2 / 5.));
        }
        else
        {
            Ru = 1. - u * atan2(1., u);
        }

        double x{4. * lam * Ru / (M_PI * beta * beta * pf)};

        if (x < 1e-4)
        {
            return gsl_pow_3(x) * (-1. + x * (0.75 + x * (-0.6 + 0.5 * x))) / 3.;
        }
        else
        {
            return x * (1. - 0.5 * x) - log(1. + x);
        }
    });

    gsl_integration_qagiu(integ1, 0., epsabs, epsrel, limit, wsp1, &result1, &abserr1);
    dedpf += -(3. * gsl_pow_4(beta) * pf) / (8. * M_PI) * result1;

    if (pf > 1.)
    {
        auto integ2 = make_gsl_function([&](double u) {
            auto integ2_1 = make_gsl_function([&](double p) {
                if (p == 0.)
                {
                    return 0.;
                }
                double p2{p * p};
                double g0{calc_G20(u, p)};
                double g0_p0{calc_G20(u, 0)};
                double fg0{4. * M_PI * lam * g0 / (pf * p2)};
                double fg0_p0{4. * M_PI * lam * g0_p0 / (pf * p2)};
                return p2 * p * (logfunc1(fg0) - logfunc1(fg0_p0));
            });
            gsl_integration_qag(integ2_1, 0., beta, epsabs, epsrel, limit, GSL_INTEG_GAUSS31, wsp2, &result2, &abserr2);
            return result2;
        });

        gsl_integration_qagiu(integ2, 0., epsabs, epsrel, limit, wsp3, &result1, &abserr1);
        dedpf += 0.75 * pf * result1 / M_PI;

        auto integ3 = make_gsl_function([&](double u) {
            auto integ3_1 = make_gsl_function([&](double p) {
                if (p == 0.)
                {
                    return 0.;
                }
                double p2{p * p};
                double g0{calc_G20(u, p)};
                double fg0{4. * M_PI * lam * g0 / (pf * p2)};
                return p2 * p * logfunc1(fg0);
            });
            gsl_integration_qagiu(integ3_1, beta, epsabs, epsrel, limit, wsp4, &result2, &abserr2);
            return result2;
        });

        gsl_integration_qagiu(integ3, 0., epsabs, epsrel, limit, wsp5, &result1, &abserr1);
        dedpf += 0.75 * pf * result1 / M_PI;
    }
    else
    {
        auto integ2 = make_gsl_function([&](double u) {
            auto integ2_1 = make_gsl_function([&](double p) {
                if (p == 0.)
                {
                    return 0.;
                }
                double p2{p * p};
                double g0{pf * calc_G20(u / pf, p / pf)};
                double g0_p0{pf * calc_G20(u / pf, 0)};
                double fg0{4. * M_PI * lam * g0 / p2};
                double fg0_p0{4. * M_PI * lam * g0_p0 / p2};
                return p2 * p * (logfunc1(fg0) - logfunc1(fg0_p0));
            });
            gsl_integration_qag(integ2_1, 0., pf * beta, epsabs, epsrel, limit, GSL_INTEG_GAUSS31, wsp2, &result2, &abserr2);
            return result2;
        });

        gsl_integration_qagiu(integ2, 0., epsabs, epsrel, limit, wsp3, &result1, &abserr1);
        dedpf += 0.75 * result1 / (M_PI * gsl_pow_4(pf));

        auto integ3 = make_gsl_function([&](double u) {
            auto integ3_1 = make_gsl_function([&](double p) {
                if (p == 0.)
                {
                    return 0.;
                }
                double p2{p * p};
                double g0{pf * calc_G20(u / pf, p / pf)};
                double fg0{4. * M_PI * lam * g0 / p2};
                return p2 * p * logfunc1(fg0);
            });
            gsl_integration_qagiu(integ3_1, pf * beta, epsabs, epsrel, limit, wsp4, &result2, &abserr2);
            return result2;
        });

        gsl_integration_qagiu(integ3, 0., epsabs, epsrel, limit, wsp5, &result1, &abserr1);
        dedpf += 0.75 * result1 / (M_PI * gsl_pow_4(pf));
    }

    double dedr{-pf * dedpf / rs};

    return dedr;
}

double energy_int_diff(double u, double p, double thr_u, double thr_p, double pf)
{
    if (p == 0.)
    {
        return 0.;
    }

    double Ci{calc_c(u, p, thr_u, thr_p)};
    double G2{calc_G20(u, p)};

    double p2{p * p};

    double p_y{lam * sqrt(4. * M_PI * Ci) / pf};
    double y{p_y / p};
    double p2_x{lam * 4. * M_PI * G2 / pf};
    double x{p2_x / p2};

    double f;

    if (fabs(x) < 1e-4 && fabs(y) < 1e-4)
    {
        double y2{y * y};
        double x_1 = x + 1.;
        double y2_x1 = y2 / x_1;
        f = (y2_x1 / (3. * x_1)) * (x + y2_x1 * ((15. + x * (9. + x * (5. + x))) / 30. - y2_x1 * (252. + x * (243. + x * (126. + x * (34. + 4. * x)))) / 945.));
    }
    else
    {
        double p2_log_1;
        double p2_log_2;

        if (fabs(y) < 1e-4)
        {
            double y2{y * y};
            p2_log_1 = p2_x * y2 * (-1. / 3. + y2 * ((12. + 7. * x) / (90. * (1. + x)) + y2 * (-153. - 180. * x - 62. * x * x) / (2835. * gsl_pow_2(1. + x)))) / (1. + x);
        }
        else
        {
            double p2{p * p};
            p2_log_1 = p2 * log((p2 + p2_x * tanh_x(y)) / (p2 + p2_x));
        }

        if (fabs(y) < 1e-4)
        {
            double y2{y * y};
            p2_log_2 = p * p * y2 * (1. / 2. - y2 * (1. / 12. - y2 / 45.));
        }
        else if (fabs(y) < 10)
        {
            p2_log_2 = p * p * log(cosh(y));
        }
        else
        {
            p2_log_2 = p * p * log(0.5 * (1. + exp(-2. * y))) + p * p_y;
        }

        double p2_d{p * p * (x * x * tanh_x_1(y) - (1. + x) * y * y * tanh_x(y)) / ((1. + x) * (1. + x * tanh_x(y)))};

        f = 2. * p2_log_1 + 2. * p2_log_2 + p2_d;
    }

    return f;
}

double calc_energy_diff(double rs)
{
    double pf{pow(2.25 * M_PI, 1. / 3.) / rs};

    double dEgbdr{calc_energy_gb_diff(rs)};

    double thr_p{1e-4};
    double thr_u{1e-4};
    double epsabs = 1e-7;
    double epsrel = 1e-6;
    size_t limit = 100;
    double result1, result2;
    double abserr1, abserr2;

    IntegrationWorkspace wsp1(limit);
    IntegrationWorkspace wsp2(limit);

    auto integ = make_gsl_function([&](double p) {
        auto integ_inner = make_gsl_function([&](double u) {
            return energy_int_diff(u, p, thr_u, thr_p, pf);
        });
        gsl_integration_qagiu(integ_inner, 0, epsabs, epsrel, limit, wsp2, &result2, &abserr2);
        return p * result2;
    });

    gsl_integration_qagiu(integ, 0, epsabs, epsrel, limit, wsp1, &result1, &abserr1);

    double dEdpf{result1 * 3. * pf / (4. * M_PI)};
    double dEdr{dEdpf * (-pf) / rs};
    return dEgbdr + dEdr;
}