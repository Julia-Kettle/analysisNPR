#ifndef FIT_FUNCTIONS_H
#define FIT_FUNCTIONS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "AnalysisNPR.h"
#include <iostream>
#include <vector>


/////////////////////////// A + B/p^2 /////////////////////////////
int psq_inv_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

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
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

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
//plain function
double psq_inv(const gsl_vector *p, double x0)
{
    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    return p0 + p1/pow(x0,2);
}

/////////////////////////// A + B/p^6 /////////////////////////////
int psix_inv_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    //calculate model value and put residual into gsl vector
    for ( size_t i=0; i < n; i++)
    {
        double Yi = p0 + p1/(pow(x[i],6));
        gsl_vector_set (f, i, Yi - y[i]);
    }
    return GSL_SUCCESS;
}

//derivative
int psix_inv_df(const gsl_vector *p, void *data, gsl_matrix *J)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    //Jacobian matrix J(i,j) = dfi / dpj
    // fi = (Yi - yi)/sigma[i]
    for ( size_t i=0; i<n; i++ )
    {
        gsl_matrix_set(J,i,0,1.0);
        gsl_matrix_set(J,i,1,1.0/pow(x[i],6));
    }
    return GSL_SUCCESS;
}
//plain function
double psix_inv(const gsl_vector *p, double x0)
{
    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    return p0 + p1/pow(x0,6);
}

/////////////////////////// A + +B/p^2 + C/p^6 /////////////////////////////
int psq_psix_inv_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);
    double p2 = gsl_vector_get(p,2);

    //calculate model value and put residual into gsl vector
    for ( size_t i=0; i < n; i++)
    {
        double Yi = p0 + p1/pow(x[i],2) + p2/pow(x[i],6);
        gsl_vector_set (f, i, Yi - y[i]);
    }
    return GSL_SUCCESS;
}

//derivative
int psq_psix_inv_df(const gsl_vector *p, void *data, gsl_matrix *J)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);
    double p2 = gsl_vector_get(p,2);

    //Jacobian matrix J(i,j) = dfi / dpj
    // fi = (Yi - yi)/sigma[i]
    for ( size_t i=0; i<n; i++ )
    {
        gsl_matrix_set(J,i,0,1.0);
        gsl_matrix_set(J,i,1,1.0/pow(x[i],2));
        gsl_matrix_set(J,i,2,1.0/pow(x[i],6));
    }
    return GSL_SUCCESS;
}

//plain function
double psq_psix_inv(const gsl_vector *p, double x0)
{
    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);
    double p2 = gsl_vector_get(p,2);

    return p0 + p1/pow(x0,2) + p2/pow(x0,6);
}


/*************************** A + B*p^2 *****************************/
int psq_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    //calculate model value and put residual into gsl vector
    for ( size_t i=0; i < n; i++)
    {
        double Yi = p0 + p1*(x[i]*x[i]);
        gsl_vector_set (f, i, Yi - y[i]);
    }
    return GSL_SUCCESS;
}

//derivative
int psq_df(const gsl_vector *p, void *data, gsl_matrix *J)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    //Jacobian matrix J(i,j) = dfi * dpj
    // fi = (Yi - yi)*sigma[i]
    for ( size_t i=0; i<n; i++ )
    {
        gsl_matrix_set(J,i,0,1.0);
        gsl_matrix_set(J,i,1,1.0*(x[i]*x[i]));
    }
    return GSL_SUCCESS;
}
//plain function
double psq(const gsl_vector *p, double x0)
{
    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    return p0 + p1*pow(x0,2);
}

/*************************** A + B*p^6 *****************************/
int psix_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    //calculate model value and put residual into gsl vector
    for ( size_t i=0; i < n; i++)
    {
        double Yi = p0 + p1*(pow(x[i],6));
        gsl_vector_set (f, i, Yi - y[i]);
    }
    return GSL_SUCCESS;
}

//derivative
int psix_df(const gsl_vector *p, void *data, gsl_matrix *J)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    //Jacobian matrix J(i,j) = dfi * dpj
    // fi = (Yi - yi)*sigma[i]
    for ( size_t i=0; i<n; i++ )
    {
        gsl_matrix_set(J,i,0,1.0);
        gsl_matrix_set(J,i,1,1.0*pow(x[i],6));
    }
    return GSL_SUCCESS;
}
//plain function
double psix(const gsl_vector *p, double x0)
{
    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    return p0 + p1*pow(x0,6);
}

/*************************** A + +B*p^2 + C*p^6 *****************************/
int psq_psix_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);
    double p2 = gsl_vector_get(p,2);

    //calculate model value and put residual into gsl vector
    for ( size_t i=0; i < n; i++)
    {
        double Yi = p0 + p1*pow(x[i],2) + p2*pow(x[i],6);
        gsl_vector_set (f, i, Yi - y[i]);
    }
    return GSL_SUCCESS;
}

//derivative
int psq_psix_df(const gsl_vector *p, void *data, gsl_matrix *J)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);

    //Jacobian matrix J(i,j) = dfi * dpj
    // fi = (Yi - yi)*sigma[i]
    for ( size_t i=0; i<n; i++ )
    {
        gsl_matrix_set(J,i,0,1.0);
        gsl_matrix_set(J,i,1,1.0*pow(x[i],2));
        gsl_matrix_set(J,i,2,1.0*pow(x[i],6));
    }
    return GSL_SUCCESS;
}

//plain function
double psq_psix(const gsl_vector *p, double x0)
{
    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);
    double p2 = gsl_vector_get(p,2);

    return p0 + p1*pow(x0,2) + p2*pow(x0,6);
}

//Should really change these to take psq as input

/*************************** A + B*p^2 + C*p^4 + D*p^6 *****************************/
int polysq_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);
    double p2 = gsl_vector_get(p,2);
    double p3 = gsl_vector_get(p,3);

    //calculate model value and put residual into gsl vector
    for ( size_t i=0; i < n; i++)
    {
        double Yi = p0 + p1*pow(x[i],2) + p2*pow(x[i],4) + p3*pow(x[i],6);
        gsl_vector_set (f, i, Yi - y[i]);
    }
    return GSL_SUCCESS;
}

//derivative
int polysq_df(const gsl_vector *p, void *data, gsl_matrix *J)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;


    //Jacobian matrix J(i,j) = dfi * dpj
    // fi = (Yi - yi)*sigma[i]
    for ( size_t i=0; i<n; i++ )
    {
        gsl_matrix_set(J,i,0,1.0);
        gsl_matrix_set(J,i,1,1.0*pow(x[i],2));
        gsl_matrix_set(J,i,2,1.0*pow(x[i],4));
        gsl_matrix_set(J,i,3,1.0*pow(x[i],6));
    }
    return GSL_SUCCESS;
}

//plain function
double polysq(const gsl_vector *p, double x0)
{
    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);
    double p2 = gsl_vector_get(p,2);
    double p3 = gsl_vector_get(p,3);

    return p0 + p1*pow(x0,2) + p2*pow(x0,4) + p3*pow(x0,6);
}


//Should really change these to take psq as input

/*************************** A + B/p^2 + C/p^4 + D/p^6 *****************************/
int inv_polysq_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);
    double p2 = gsl_vector_get(p,2);
    double p3 = gsl_vector_get(p,3);

    //calculate model value and put residual into gsl vector
    for ( size_t i=0; i < n; i++)
    {
        double Yi = p0 + p1/pow(x[i],2) + p2/pow(x[i],4) + p3/pow(x[i],6);
        gsl_vector_set (f, i, Yi - y[i]);
    }
    return GSL_SUCCESS;
}

//derivative
int inv_polysq_df(const gsl_vector *p, void *data, gsl_matrix *J)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;


    //Jacobian matrix J(i,j) = dfi * dpj
    // fi = (Yi - yi)*sigma[i]
    for ( size_t i=0; i<n; i++ )
    {
        gsl_matrix_set(J,i,0,1.0);
        gsl_matrix_set(J,i,1,1.0/pow(x[i],2));
        gsl_matrix_set(J,i,2,1.0/pow(x[i],4));
        gsl_matrix_set(J,i,3,1.0/pow(x[i],6));
    }
    return GSL_SUCCESS;
}

//plain function
double inv_polysq(const gsl_vector *p, double x0)
{
    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);
    double p2 = gsl_vector_get(p,2);
    double p3 = gsl_vector_get(p,3);

    return p0 + p1/pow(x0,2) + p2/pow(x0,4) + p3/pow(x0,6);
}

/*************************** A + B*p^2 + C/p^2  *****************************/
int psq_psq_inv_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;

    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);
    double p2 = gsl_vector_get(p,2);

    //calculate model value and put residual into gsl vector
    for ( size_t i=0; i < n; i++)
    {
        double Yi = p0 + p1*pow(x[i],2) + p2/pow(x[i],2);
        gsl_vector_set (f, i, Yi - y[i]);
    }
    return GSL_SUCCESS;
}

//derivative
int psq_psq_inv_df(const gsl_vector *p, void *data, gsl_matrix *J)
{
    //unpack data from struct
    size_t n    =   ((struct DataSet*)data)->n;
    std::vector<double> x   =   ((struct DataSet*)data)->x;
    std::vector<double> y   =   ((struct DataSet*)data)->y;


    //Jacobian matrix J(i,j) = dfi * dpj
    // fi = (Yi - yi)*sigma[i]
    for ( size_t i=0; i<n; i++ )
    {
        gsl_matrix_set(J,i,0,1.0);
        gsl_matrix_set(J,i,1,1.0*pow(x[i],2));
        gsl_matrix_set(J,i,2,1.0/pow(x[i],2));
    }
    return GSL_SUCCESS;
}

//plain function
double psq_psq_inv(const gsl_vector *p, double x0)
{
    // get parameters
    double p0 = gsl_vector_get(p,0);
    double p1 = gsl_vector_get(p,1);
    double p2 = gsl_vector_get(p,2);

    return p0 + p1*pow(x0,2) + p2/pow(x0,2);
}


//Should really change these to take psq as input


#endif
