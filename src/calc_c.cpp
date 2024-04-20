#include "common_var.h"
#include "integ_prep.h"
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <fstream>

std::vector<double> set_npt(double a, double b, std::vector<double> npt_cand)
{
    std::vector<double> v{a, b};

    for (auto pt : npt_cand)
    {
        if (a < pt && pt < b)
        {
            v.push_back(pt);
        }
    }

    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());

    return v;
}

double integrand1(double pz1, double pz2, double p, double u)
{
    double eps{pz2 - pz1};
    double aeps{fabs(eps)};

    double l12;
    double log_1;
    double e2G;

    double b2{fabs(1 - (pz1 - 0.5 * p) * (pz1 - 0.5 * p))};

    if (eps >= 0)
    {
        l12 = (1. - pz2 + 0.5 * p) * (1. + pz1 - 0.5 * p);
        log_1 = xlogy(eps * eps * b2, aeps / (1. - pz1 + 0.5 * p));
        e2G = eps * eps * l12 + 2. * log_1;
    }
    else
    {
        l12 = (1. - pz1 + 0.5 * p) * (1. + pz2 - 0.5 * p);
        log_1 = xlogy(eps * eps * b2, aeps / (1. + pz1 - 0.5 * p));
        e2G = eps * eps * l12 + 2. * log_1;
    }

    double u2{u * u};
    double ker{(gsl_pow_2(u2 - pz1 * pz2) - u2 * gsl_pow_2(pz1 + pz2)) / (gsl_pow_2((u2 + pz1 * pz1) * (u2 + pz2 * pz2)))};

    return ker * e2G;
}

double integrand2(double pz1, double pz2, double p, double u)
{
    double eps{pz2 + pz1};
    double eps2{eps * eps};

    double t11{(eps - p) * (2. * pz1 - eps)};
    double l12{0.5 * (2. + p * eps + 2. * pz1 * pz2 - p * p * 0.5 - sqrt(t11 * t11 + 2. * eps2 * (2. - gsl_pow_2(pz1 - 0.5 * p) - gsl_pow_2(pz2 - 0.5 * p) + 0.5 * eps2)))};
    double log_1;
    double e2G;

    double b2{fabs(1. - gsl_pow_2(pz1 - 0.5 * p))};

    double t1{t11 + eps2};
    if (t1 >= 0)
    {
        log_1 = xlogy(eps2 * b2, 2. * eps2 / (sqrt(t11 * t11 + 2. * eps2 * (2. - gsl_pow_2(pz1 - 0.5 * p) - gsl_pow_2(pz1 + 0.5 * p) + 2. * (pz1 + 0.5 * p) * eps - eps2 * 0.5)) + t1));
    }
    else
    {
        log_1 = xlogy(eps2 * b2, (sqrt(t11 * t11 + 2. * eps2 * (2. - gsl_pow_2(pz1 - 0.5 * p) - gsl_pow_2(pz1 + 0.5 * p) + 2. * (pz1 + 0.5 * p) * eps - eps2 * 0.5)) - t1) / (2. * (1. - gsl_pow_2(pz1 - 0.5 * p))));
    }

    e2G = eps2 * l12 + 2. * log_1;

    double u2{u * u};
    double ker{(gsl_pow_2(u2 + pz1 * pz2) - u2 * gsl_pow_2(pz1 - pz2)) / (gsl_pow_2((u2 + pz1 * pz1) * (u2 + pz2 * pz2)))};
    return ker * e2G;
}

double calc_c(double u, double p, double thr_u, double thr_p)
{
    if (p < 1e-4)
    {
        if (u == 0)
        {
            return (Ns / 2.) * gsl_pow_3(M_1_PI);
        }
        if (u > 1.e2)
        {
            double u2{u * u};
            double uinv2{1. / u2};
            double xatanx{1. + uinv2 * (-1. / 3. + uinv2 / 5.)};
            double xatanx_m1{uinv2 * (-1. / 3. + uinv2 / 5.)};
            return (Ns / 2.) * M_1_PI * M_1_PI * M_1_PI * (1. - 2. * xatanx - 3. * xatanx_m1 * u2) / (1. + u2);
        }
        return (Ns / 2.) * M_1_PI * M_1_PI * M_1_PI * (1. + 3. * u * u - (2. + 3. * u * u) * u * atan(1. / u)) / (1. + u * u);
    }

    u = (u > thr_u) ? u : thr_u;

    double epsabs{1e-9};
    double epsrel{1e-8};
    size_t limit{10000};

    double result{0.};
    double result1, result2, result3, result4;
    double abserr1, abserr2, abserr3, abserr4;

    IntegrationWorkspace wsp1(limit);
    IntegrationWorkspace wsp2(limit);
    IntegrationWorkspace wsp3(limit);
    IntegrationWorkspace wsp4(limit);

    {
        auto integ1 = make_gsl_function([&](double pz1) {
            auto integ1_1 = make_gsl_function([&](double pz2) {
                return integrand1(pz1, pz2, p, u);
            });

            std::vector<double> npt_cand{-pz1, 0, pz1};
            std::vector<double> pts_vec{set_npt(-1. + 0.5 * p, 1. + 0.5 * p, npt_cand)};
            size_t npts{pts_vec.size()};
            double *pts{&pts_vec[0]};

            try
            {
                int err{gsl_integration_qagp(&integ1_1, pts, npts, epsabs, epsrel, limit, wsp1, &result1, &abserr1)};
                if (err)
                    throw err;
            }
            catch (int err)
            {
                std::cout << "integ1_1: " << err << std::endl;
                std::cout << "u=" << u << ", p=" << p << std::endl;
                std::cout << "result=" << result1 << std::endl;
                std::cout << "abserr=" << abserr1 << std::endl;
            }
            return result1;
        });

        std::vector<double> npt_cand{0, 1. - 0.5 * p};
        std::vector<double> pts_vec{set_npt(-1. + 0.5 * p, 1. + 0.5 * p, npt_cand)};
        size_t npts{pts_vec.size()};
        double *pts{&pts_vec[0]};

        try
        {
            int err{gsl_integration_qagp(&integ1, pts, npts, epsabs, epsrel, limit, wsp2, &result2, &abserr2)};
            if (err)
                throw err;
        }
        catch (int err)
        {
            std::cout << "integ1: " << err << std::endl;
            std::cout << "u=" << u << ", p=" << p << std::endl;
            std::cout << "result=" << result2 << std::endl;
            std::cout << "abserr=" << abserr2 << std::endl;
        }
        result += result2;
    }

    {
        auto integ2 = make_gsl_function([&](double pz1) {
            auto integ2_1 = make_gsl_function([&](double pz2) {
                return integrand2(pz1, pz2, p, u);
            });

            std::vector<double> npt_cand{-pz1, 0, pz1};
            std::vector<double> pts_vec{set_npt(-1. + 0.5 * p, 1. + 0.5 * p, npt_cand)};
            size_t npts{pts_vec.size()};
            double *pts{&pts_vec[0]};
            try
            {
                int err{gsl_integration_qagp(&integ2_1, pts, npts, epsabs, epsrel, limit, wsp3, &result3, &abserr3)};
                if (err)
                    throw err;
            }
            catch (int err)
            {
                std::cout << "integ2_1: " << err << std::endl;
                std::cout << "u=" << u << ", p=" << p << std::endl;
                std::cout << "result=" << result3 << std::endl;
                std::cout << "abserr=" << abserr3 << std::endl;
            }
            return result3;
        });

        std::vector<double> npt_cand{0, 1. - 0.5 * p};
        std::vector<double> pts_vec{set_npt(-1. + 0.5 * p, 1. + 0.5 * p, npt_cand)};
        size_t npts{pts_vec.size()};
        double *pts{&pts_vec[0]};

        try
        {
            int err{gsl_integration_qagp(&integ2, pts, npts, epsabs, epsrel, limit, wsp4, &result4, &abserr4)};
            if (err)
                throw err;
        }
        catch (int err)
        {
            std::cout << "integ2: " << err << std::endl;
            std::cout << "u=" << u << ", p=" << p << std::endl;
            std::cout << "result=" << result4 << std::endl;
            std::cout << "abserr=" << abserr4 << std::endl;
        }

        result -= result4;
    }

    result *= (double)Ns / (16. * gsl_pow_3(M_PI) * gsl_pow_2(p));

    return (result > 0) ? result : 0.;
}
