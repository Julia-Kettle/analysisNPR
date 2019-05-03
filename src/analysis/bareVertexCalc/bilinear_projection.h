#ifndef BILINEAR_PROJECTION_H
#define BILINEAR_PROJECTION_H

#include "Grid/Grid.h"
#include <string>
#include <iostream>
#include <Grid/Eigen/Core>
#include <Grid/Eigen/Dense>
#include <Grid/Eigen/SVD>
#include <complex>
#include "distribution/distribution.h"
#include "maths/maths.h"

//////////////////////////////////////////////////////////////////////////
// Project amputated vertex functions in the gamma scheme
// S  - 1/12 Tr ( LambdaS * I ) 
// P  - 1/12 Tr ( LambdaP * g5 )
// V  - 1/48 Tr ( LambdaV^mu * gmu )
// A  - 1/48 Tr ( LambdaA^mu * gmu * g5 )
/////////////////////////////////////////////////////////////////////////
template <typename T>
Distribution<Grid::Real> project_gamma(std::vector<Distribution<T>> amp_vertex, std::vector<Grid::QCD::Gamma::Algebra> gamma_indices)
{
    
    Distribution<Grid::QCD::ComplexD> tr;
    for (int mu=0;mu<gamma_indices.size();mu++)
    {
        auto gi = gamma_indices[mu];
        if( mu == 0)
        {
            tr = trace( amp_vertex[gi] * Grid::QCD::Gamma(gi) ); 
        }
        else
        {
            tr = tr +  trace( amp_vertex[gi] * Grid::QCD::Gamma(gi) );
        }
    }
    return real(tr*(1/(12.0*gamma_indices.size())));
}

template <typename T>
Grid::RealD project_gamma(std::vector<T> amp_vertex, std::vector<Grid::QCD::Gamma::Algebra> gamma_indices)
{
    
    Grid::QCD::ComplexD tr;
    for (int mu=0;mu<gamma_indices.size();mu++)
    {
        auto gi = gamma_indices[mu];
        if( mu == 0)
        {
            tr = trace( amp_vertex[gi] * Grid::QCD::Gamma(gi) ); 
        }
        else
        {
            tr = tr +  trace( amp_vertex[gi] * Grid::QCD::Gamma(gi) );
        }
    }
    return real(tr*(1/(12.0*gamma_indices.size())));
}


//////////////////////////////////////////////////////////////////////////
// Project amputated vertex functions in the qslash scheme
// S  - 1/12 Tr ( Lambda S * I ) 
// P  - 1/12 Tr ( Lambda P * g5 )
// V  - 1/12q^2 Tr ( qmu * LambdaV^mu  * qslash ) = 1/12q^2 Tr ( q_mu * LambdaV_mu * gamma_nu * q_nu )
// A  - 1/12q^2 Tr ( qmu * LambdaA^mu * g5 * qslash ) = 1/12q^2 Tr ( q_mu * LambdaV_mu * gamma5 & gamma_nu * q_nu )
/////////////////////////////////////////////////////////////////////////
template <typename T>
Distribution<Grid::Real> project_qslash(std::vector<Distribution<T>> amp_vertex, std::vector<double> q, std::vector<Grid::QCD::Gamma::Algebra> gamma_indices)
{
   
    //// qsq ///////////////////////////////////
    double qsq = 0;
    for(int mu=0;mu<gamma_indices.size();mu++)
    {
        qsq += q[mu]*q[mu];
    }


    //// Trace ///////////////////////////////////
    // needed only for V and A
    // q[mu] and qsq s cancel for P and S, but can just use the gamma scheme
    Distribution<Grid::QCD::ComplexD> tr;
    Distribution<Grid::QCD::ComplexD> tr_mu;

    for(int mu=0;mu<gamma_indices.size();mu++)
    {
        //convention for gamma indexing and dimension indexing differ
        // Gamma: T X Y Z   vs   X Y Z T
        // 0 1 2 3    ->    3 0 1 2  ( txyz )
        int mu_dim = mu; //(mu+3)%4; 
        auto gmu = gamma_indices[mu];
            
        for (int nu=0;nu<gamma_indices.size();nu++)
        {

            //indexing for qslash 
            auto gnu = gamma_indices[nu];
            int  nu_dim = nu;//(nu+3)%4; // 0 -> 3, 1->0, etc
            

            // Final result is :  trace( q[mu]*PI[mu] * sum (Gamma[nu]*q[nu]) ) = trace( q[mu]*PI[mu]*qslash )
            // here we 
            auto proj_nu = q[mu_dim] * amp_vertex[gmu] * Grid::QCD::Gamma(gnu) * q[nu_dim]; 

            // take trace over nu for mu component 
            (nu == 0) ? tr_mu = trace( proj_nu ) : tr_mu = tr_mu + trace(proj_nu);
        }

        // trace over mu 
        (mu == 0 ) ? tr = tr_mu : tr = tr +tr_mu;
    }
    return real(tr*(1/(12.0*qsq)));
}

template <typename T>
Grid::Real project_qslash(std::vector<T> amp_vertex, std::vector<double> q, std::vector<Grid::QCD::Gamma::Algebra> gamma_indices)
{
   
    //// qsq ///////////////////////////////////
    double qsq = 0;
    for(int mu=0;mu<gamma_indices.size();mu++)
    {
        qsq += q[mu]*q[mu];
    }

    //// Trace ///////////////////////////////////
    // needed only for V and A
    // q[mu] and qsq s cancel for P and S, but can just use the gamma scheme
    Grid::QCD::ComplexD tr;
    for(int mu=0;mu<gamma_indices.size();mu++)
    {
        auto gi = gamma_indices[mu];
        if(mu==0)
        {
            tr = trace( q[mu] * amp_vertex[gi] * Grid::QCD::Gamma(gi) * q[mu] );
        } 
        else 
        {
            tr = tr +  trace( q[mu] * amp_vertex[gi] * Grid::QCD::Gamma(gi) * q[mu] );
        }
    }
    return real(tr*(1/(12.0*qsq)));
}

#endif
