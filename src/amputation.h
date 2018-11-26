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

////////////////////////////////////////////////////
// Invert the volume averaged propogators
// Rewrite as 12x12 matrix for inverting 
// then restructure as spin-colour matrix after
////////////////////////////////////////////////////
Grid::QCD::SpinColourMatrix invert(Grid::QCD::SpinColourMatrix sc_matrix)
{
    int Ns(Grid::QCD::Ns);
    int Nc(Grid::QCD::Nc);
    
    Grid::QCD::SpinColourMatrix sc_matrixInv;
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
Distribution<Grid::QCD::SpinColourMatrix> invert(Distribution<Grid::QCD::SpinColourMatrix> dist_scmat)
{
    std::vector<Grid::QCD::SpinColourMatrix> vec_Inv;
    vec_Inv.reserve(dist_scmat.get_values().size());
    for ( auto scmat : dist_scmat.get_values() )
    {
        vec_Inv.push_back(invert(scmat));
    }
    return Distribution<Grid::QCD::SpinColourMatrix>(vec_Inv);
}

///////////////////////////////////////////////
//Amputation code ( no projection in this function )
// S-1 V S
//////////////////////////////////////////////
std::vector<Distribution<Grid::QCD::SpinColourMatrix>> amputate(Distribution<Grid::QCD::SpinColourMatrix> prop1, Distribution<Grid::QCD::SpinColourMatrix> prop2, std::vector<Distribution<Grid::QCD::SpinColourMatrix>> vertex)
{
    //invert propagators
    auto propInv1 = invert(prop1);
    auto propInv2 = invert(prop2);
    // set up vector to hold amputated result
    std::vector<Distribution<Grid::QCD::SpinColourMatrix>> amputated;

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
// V  - 1/12q^2 Tr ( qmu * LambdaV^mu  * qslash )
// A  - 1/12q^2 Tr ( qmu * LambdaA^mu * g5 * qslash )
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
    std::cout << "tr - " << tr.get_values() << std::endl;
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
