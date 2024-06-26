

#ifndef COULOMB3D_2_INTEG_PREP
#define COULOMB3D_2_INTEG_PREP

#include <iostream>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <vector>
#include <algorithm>

// Simple RAII wrapper
class IntegrationWorkspace
{
    gsl_integration_workspace *wsp;

public:
    IntegrationWorkspace(const size_t n = 1000) : wsp(gsl_integration_workspace_alloc(n)) {}
    ~IntegrationWorkspace() { gsl_integration_workspace_free(wsp); }

    operator gsl_integration_workspace *() { return wsp; }
};

// Build gsl_function from lambda
template <typename F>
class gsl_function_pp : public gsl_function
{
    const F func;
    static double invoke(double x, void *params)
    {
        return static_cast<gsl_function_pp *>(params)->func(x);
    }

public:
    gsl_function_pp(const F &f) : func(f)
    {
        function = &gsl_function_pp::invoke; // inherited from gsl_function
        params = this;                       // inherited from gsl_function
    }
    operator gsl_function *() { return this; }
};

// Helper function for template construction
template <typename F>
gsl_function_pp<F> make_gsl_function(const F &func)
{
    return gsl_function_pp<F>(func);
}

std::vector<double> set_npt(double a, double b, std::vector<double> npt_cand);

#endif // COULOMB3D_2_INTEG_PREP
