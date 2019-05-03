#ifndef AMPUTATION_H
#define AMPUTATION_H

#include <Grid/Grid.h>
#include <string>
#include <iostream>
#include <Grid/Eigen/Core>
#include <Grid/Eigen/Dense>
#include <Grid/Eigen/SVD>
#include <complex>
#include "distribution/distribution.h"


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


#endif
