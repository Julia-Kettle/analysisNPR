#ifndef FITTER_H
#define FITTER_H


#include "AnalysisNPR.h"
#include <iostream>
#include <vector>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit_nlin.h>

////////////////////////////////////////////////////////////////////////
// Kind of hacky wrapper class for gsl minimiser. 
// have to set fitter.fdf.f and fitter.fdf.df in the main
// data attached in construction. Can then fit multiple times with same
// fitter as long as f and df set correctly. 
// Plan to rewrite to make less hack but will do for now

struct DataSet
{
    std::vector<double> y;
    std::vector<double> x;
    std::vector<double> sigma;
    size_t              n;
};

////////////////////////////////

////////////////////////////////////////////////////
// wrapper class for GSL minimiser for least sq
// requires user 
class Fitter
{
    private:
        DataSet data;
    public:
        Fitter(DataSet d);
        int minimise();
        void printResults();
        void free();

        // variables
        size_t  n, np;
        double chi, chi0;
        int status, info;

        std::vector<double> p_init, weights;

        const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
        gsl_multifit_fdfsolver *s;
        gsl_multifit_function_fdf fdf;
        gsl_vector *res_f;
        gsl_matrix *J, *covar;
        
        gsl_vector_view sgm, wght,pars;

        // set tolerances 
        const double xtol = 1e-8;
        const double gtol = 1e-8;
        const double ftol = 0.0;

};

Fitter::Fitter(DataSet d)
{
    this->data      = d;
    this->n         = d.n;
    this->weights   = std::vector<double>(n,1.0);
    
    //set up weights and errors
    sgm = gsl_vector_view_array(&(data.sigma[0]), n);
    wght = gsl_vector_view_array(&weights[0], n);
    
}

void Fitter::free()
{
    // free up the solvers and covariance and Jacobian matrices 
    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    gsl_matrix_free (J);

}

int Fitter::minimise()
{
    J = gsl_matrix_alloc(n, np);
    covar = gsl_matrix_alloc (np, np);
    pars = gsl_vector_view_array (&p_init[0], np);

    fdf.n = n;
    fdf.p = np;
    fdf.params = &data;
    s = gsl_multifit_fdfsolver_alloc (T, n, np);


    /* initialize solver with starting point and weights */
    gsl_multifit_fdfsolver_wset (s, &fdf, &pars.vector, &wght.vector);

    /* compute initial residual norm */
    res_f = gsl_multifit_fdfsolver_residual(s);

    chi0=0;
    for (size_t i=0;i<res_f->size;i++)
    {
        chi0 += gsl_vector_get(res_f,i)*gsl_vector_get(res_f,i);
    }
    chi0 = Grid::sqrt(chi0);

    /* solve the system with a maximum of 20 iterations */
    status = gsl_multifit_fdfsolver_driver(s, 20, xtol, gtol, ftol, &info);

    gsl_multifit_fdfsolver_jac(s, J);
    gsl_multifit_covar (J, 0.0, covar);

    /* compute final residual norm */
    chi=0;
    for (size_t i=0;i<res_f->size;i++)
    {
        chi += gsl_vector_get(res_f,i)*gsl_vector_get(res_f,i);
    }
    chi = Grid::sqrt(chi);

}

void Fitter::printResults()
{

    #define FIT(i) gsl_vector_get(s->x, i)
    #define ERR(i) Grid::sqrt(gsl_matrix_get(covar,i,i))

    fprintf(stderr, "summary from method '%s'\n",
                    gsl_multifit_fdfsolver_name(s));
    fprintf(stderr, "number of iterations: %zu\n",
                      gsl_multifit_fdfsolver_niter(s));
    fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n",
                    (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %g\n", chi0);
    fprintf(stderr, "final   |f(x)| = %g\n", chi);

    { 
        double dof = n - np;
        double c = GSL_MAX_DBL(1, chi / Grid::sqrt(dof)); 

        fprintf(stderr, "chisq/dof = %g\n",  pow(chi, 2.0) / dof);
        
        for( size_t i=0; i<np; i++)
        {
            fprintf (stderr, "p(%i)      = %.5f +/- %.5f\n", i, FIT(i), c*ERR(i));
        }
    }

fprintf (stderr, "status = %s\n", gsl_strerror (status));
}



// hacky jackknife fitter surely a more efficient way
class DistributionFitter
{
    private:
        std::vector<Fitter>                 fitters;
        int                                 nSamples;
        int                                 nParams;
        Distribution<std::vector<double>>   y;
        std::vector<double>                 x;
        std::vector<double>                 sigma;
        std::vector<DataSet>                data;
        std::string                         resamplingType;

        std::vector<std::vector<double>>    params;
        std::vector<double>                 chi;

        double          (*function)(const gsl_vector *, double);
    public:
        DistributionFitter(Distribution<std::vector<double>> y, std::vector<double> x);
        void assignFitFunction( int (*f)(const gsl_vector *, void *, gsl_vector *),  int(*df)(const gsl_vector *, void *, gsl_matrix *), double(*func)(const gsl_vector *, double), std::vector<double> p_init );
        void fitAll();
        void fit(int i);

        // get functions
        std::vector<std::vector<double>>    get_params(){ return params; }
        std::vector<double>                 get_chi(){ return chi; }
        
        //extrapolate function
        std::vector<double>                 extrapolate(double x0); 

};


DistributionFitter::DistributionFitter(Distribution<std::vector<double>> y, std::vector<double> x)
{
    // Create vector of structs and fitters
    nSamples = y.get_values().size();
    data.resize(nSamples);

    for(int i=0; i<nSamples; i++)
    {
        data[i].y       = y.get_value(i);
        data[i].x       = x;
        data[i].sigma   = y.get_std();
        data[i].n       = data[i].y.size();
           
        //fitter
        fitters.push_back(Fitter(data[i]));
    }

    resamplingType = y.get_resamplingType();

    chi.resize(nSamples);
    params.resize(nSamples);
}

void DistributionFitter::assignFitFunction( int (*f)(const gsl_vector *, void *, gsl_vector *), int(*df)(const gsl_vector *, void *, gsl_matrix *), double(*func)(const gsl_vector *, double), std::vector<double> p_init )
{
    nParams = p_init.size();
    for (int i=0;i<nSamples;i++)
    {
        fitters[i].fdf.f  = f;
        fitters[i].fdf.df = df;
        fitters[i].np = nParams;
        fitters[i].p_init = p_init;
    }
    function = func;
}

void DistributionFitter::fit(int i)
{
    fitters[i].minimise();
}

void DistributionFitter::fitAll()
{
    for (int i=0; i<nSamples; i++)
    {
        fitters[i].minimise();    
        chi[i] = (fitters[i].chi);
        chi[i] = (fitters[i].chi);
        for (int j=0; j<nParams; j++)
        {
            params[i].push_back(gsl_vector_get(fitters[i].s->x, j));
        }
    }
}

std::vector<double>  DistributionFitter::extrapolate(double x0)
{
    std::vector<double> extrap_vector;
    for (int i=0; i<nSamples; i++)
    {
        extrap_vector.push_back(function(fitters[i].s->x,x0));
    }
    return extrap_vector;
}

#endif
