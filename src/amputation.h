#ifndef AMPUTATION_H
#define AMPUTATION_H

#include "Grid/Grid.h"
#include <string>
#include <iostream>
#include <Grid/Eigen/Core>
#include <Grid/Eigen/Dense>
#include <Grid/Eigen/SVD>
#include <complex>
#include "distribution.h"

using namespace Grid;
using namespace QCD;


////////////////////////////////////////////////////
// Invert the volume averaged propogators
// Rewrite as 12x12 matrix for inverting 
// then restructure as spin-colour matrix after
////////////////////////////////////////////////////
SpinColourMatrix invert(SpinColourMatrix sc_matrix)
{
    SpinColourMatrix sc_matrixInv;
    Eigen::MatrixXcd matrix = Eigen::MatrixXcd::Zero(Ns*Nc,Ns*Nc);

    // convet spincolourmatrix into eigen matrix     
    for(int c1=0;c1<Nc;c1++)
    for(int c2=0;c2<Nc;c2++)
    for(int s1=0;s1<Ns;s1++)
    for(int s2=0;s2<Ns;s2++)
    {
        int i1=s1*Nc+c1;
        int i2=s2*Nc+c2;
        matrix(i1,i2) = sc_matrix()(s1,s2)(c1,c2);
    }
    
    // SVD
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(matrix,Eigen::ComputeThinU|Eigen::ComputeThinV);

    // Inversion
    Eigen::MatrixXcd matrixInv = matrix.inverse();

    // convert eigen matrix into spincolourmatrix
    for(int c1=0;c1<Nc;c1++)
    for(int c2=0;c2<Nc;c2++)
    for(int s1=0;s1<Ns;s1++)
    for(int s2=0;s2<Ns;s2++)
    {
        int i1=s1*Nc+c1;
        int i2=s2*Nc+c2;
        sc_matrixInv()(s1,s2)(c1,c2) = matrixInv(i1,i2);
    }

    return sc_matrixInv;
}

///////// Inversion for distributions //////////////
///////// wrapper for invert above    /////////////
Distribution<SpinColourMatrix> invert(Distribution<SpinColourMatrix> dist_scmat)
{
    std::vector<SpinColourMatrix> vec_Inv;
    vec_Inv.reserve(dist_scmat.get_values().size());
    for ( auto scmat : dist_scmat.get_values() )
    {
        vec_Inv.push_back(invert(scmat));
    }
    return Distribution<SpinColourMatrix>(vec_Inv);
}

///////////////////////////////////////////////
//Amputation code ( no projection in this function )
// S-1 V S
//////////////////////////////////////////////
std::vector<Distribution<SpinColourMatrix>> amputate(Distribution<SpinColourMatrix> prop1, Distribution<SpinColourMatrix> prop2, std::vector<Distribution<SpinColourMatrix>> vertex)
{
    //invert propagators
    auto propInv1 = invert(prop1);
    auto propInv2 = invert(prop2);
    // set up vector to hold amputated result
    std::vector<Distribution<SpinColourMatrix>> amputated;

    // loop through the gammas in the vertex
    for ( auto gamma : vertex )
    {
        amputated.push_back( propInv1 * gamma * propInv2 );
    }
    return amputated;
}

//////////////////////////////////////////////////////////////////////////
// Project amputated vertex functions in the gamma scheme
// S  - 1/12 Tr ( LambdaS * I ) 
// P  - 1/12 Tr ( LambdaP * g5 )
// V  - 1/48 Tr ( LambdaV^mu * gmu )
// A  - 1/48 Tr ( LambdaA^mu * gmu * g5 )
/////////////////////////////////////////////////////////////////////////
template <typename T>
auto project_gamma(std::vector<T> amp_vertex, std::vector<Gamma::Algebra> gamma_indices)
{
    
    //std::vector<Gamma::Algebra> gmu = {Gamma::Algebra::GammaT, Gamma::Algebra::GammaX, Gamma::Algebra::GammaY, Gamma::Algebra::GammaZ};
    //std::vector<Gamma::Algebra> gmug5 = {Gamma::Algebra::GammaTGamma5, Gamma::Algebra::GammaXGamma5, Gamma::Algebra::GammaYGamma5, Gamma::Algebra::GammaZGamma5}


    Distribution<ComplexD> tr;
    for (int mu=0;mu<gamma_indices.size();mu++)
    {
        auto gi = gamma_indices[mu];
        if( mu == 0)
        {
            tr = trace( amp_vertex[gi] * Gamma(gi) ); 
        }
        else
        {
            tr = tr +  trace( amp_vertex[gi] * Gamma(gi) );
        }
    }
    return tr*(1/(12.0*gamma_indices.size()));
}


//////////////////////////////////////////////////////////////////////////
// Project amputated vertex functions in the qslash scheme
// S  - 1/12 Tr ( Lambda S * I ) 
// P  - 1/12 Tr ( Lambda P * g5 )
// V  - 1/12q^2 Tr ( qmu * LambdaV^mu  * qslash )
// A  - 1/12q^2 Tr ( qmu * LambdaA^mu * g5 * qslash )
/////////////////////////////////////////////////////////////////////////
template <typename T>
auto project_qslash(std::vector<T> amp_vertex, std::vector<double> q, std::vector<Gamma::Algebra> gamma_indices)
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
    Distribution<ComplexD> tr;
    for(int mu=0;mu<gamma_indices.size();mu++)
    {
        auto gi = gamma_indices[mu];
        if(mu==0)
        {
            tr = trace( q[mu] * amp_vertex[gi] * Gamma(gi) * q[mu] );
        } 
        else 
        {
            tr = tr +  trace( q[mu] * amp_vertex[gi] * Gamma(gi) * q[mu] );
        }
    }
    return tr*(1/(12.0*qsq));
}

#endif
