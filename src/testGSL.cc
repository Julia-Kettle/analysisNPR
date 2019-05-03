#include <iostream>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <functional>

#define N   10 /* number of data-points to fit*/
#define TMAX    (40.0)


// x is x data
// data is the struct

struct data{
    size_t n;
    double *x;
    double *y;
};

//
// Here we define the residual
int psq_inv_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    //unpack data from struct
    size_t n    =   ((struct data*)data)->n;
    double *x   =   ((struct data*)data)->x;
    double *y   =   ((struct data*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    //calculate model value and put residual into gsl vector
    for ( size_t i=0; i < n; i++)
    {
        double Yi = p0 + p1/(x[i]*x[i]);
        gsl_vector_set (f, i, Yi - y[i]);
    }
    return GSL_SUCCESS;
}

//derivative
int psq_inv_df(const gsl_vector *p, void *data, gsl_matrix *J)
{
    //unpack data from struct
    size_t n    =   ((struct data*)data)->n;
    double *x   =   ((struct data*)data)->x;
    double *y   =   ((struct data*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    //Jacobian matrix J(i,j) = dfi / dpj
    // fi = (Yi - yi)/sigma[i]
    for ( size_t i=0; i<n; i++ )
    {
        gsl_matrix_set(J,i,0,1.0);
        gsl_matrix_set(J,i,1,1.0/(x[i]*x[i]));
    }
    return GSL_SUCCESS;
}


// print function to update user on fit
void callback(const size_t iter, void *params,  const gsl_multifit_nlinear_workspace *w)
{
        gsl_vector *f = gsl_multifit_nlinear_residual(w);
        gsl_vector *p = gsl_multifit_nlinear_position(w);
        double rcond;

        /* compute reciprocal condition number of J(x) */
        gsl_multifit_nlinear_rcond(&rcond, w);

        fprintf(stderr, "iter %2zu: p0 = %.4f, p1 = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
                   iter, gsl_vector_get(p, 0), gsl_vector_get(p, 1), 1.0 / rcond, gsl_blas_dnrm2(f));
}



int main()
{
    // how do we write a functions to do this?
    // we need to do a jackknife fit.
    // file with the functions and dfs
    // we have a distribution of datav( resampled )
    // bootstrap/jackknife fit - ideally pass a distribution
    // distribution is vector of y values.
    // we will have Npoints distributions. So ideally we want to pass Npoints distribution
    // but first lets just pass the actual data.  
    //
    //


}

