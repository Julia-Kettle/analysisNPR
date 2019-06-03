#ifndef DISTRIBUTION_UTILS_H
#define DISTRINUTION_UTILS_H

#include <Grid/Grid.h>
#include <string>
#include <iostream>
#include <Grid/Eigen/Core>
#include <Grid/Eigen/Dense>
#include <Grid/Eigen/SVD>
#include <complex>
#include "distribution.h"
#include "maths/maths.h"

////////////////////// trace for distributions /////////////////////////
auto trace(Distribution<Grid::QCD::SpinColourMatrix> dist)
{
    std::vector<Grid::QCD::ComplexD> vec_tr;
    for (auto value : dist.get_values() )
    {
        vec_tr.push_back(Grid::trace(value));
    }
    Distribution<Grid::ComplexD> tr(vec_tr);
    return tr;
}


////////////////////////////////////// Define conversion to Real for distribution ///////////////////////
template<typename T>
Distribution<Grid::Real> real(Distribution<T> dist){ return Distribution<Grid::Real>(real(dist.get_values())); } 



/////////////////////////////// Zeros for different data types ///////////////////////////////////
Grid::QCD::ComplexD zero(Grid::QCD::ComplexD value){ return Grid::QCD::ComplexD(0.0,0.0); }

Grid::Real zero(Grid::Real value){ return 0; }

template <typename T>
T zero(T value){ return Grid::QCD::zero; }

Eigen::MatrixXd zero(Eigen::MatrixXd value)
{
    // hacky my hackface. Surely a better way?
    return value-value;
}

template <typename T>
std::vector<T> zero(std::vector<T> value){ return std::vector<T>(value.size(),Grid::QCD::zero); }

std::vector<double> zero(std::vector<double> value){ return std::vector<double>(value.size(),0.0); }



/////////////////////////////Vector of Distributions///////////////////////////////////
// method to create vector of distribution from vector (configs) of vector (gammas)
// the vectors are switched so that we get vector (gammas) of a distribution (configs)
// typically used invertex functions where there has been looping over gamma functions
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<Distribution<T>> get_vector_distributions(std::vector<std::vector<T>> data)
{
    std::vector<std::vector<T>>  transposed = transpose(data);
    std::vector<Distribution<T>> distributions;
    distributions.reserve(transposed.size());

    for ( auto subvec : transposed )
    {
        Distribution<T> dist(subvec);
        distributions.push_back(dist);
    }
    return distributions;
}

/////////////////////////////////////////////////////////////////////////////////
// get a "matrix" of distributions. mainly only used for outputting to file
std::vector<std::vector<Distribution<double>>> get_matrix_distributions(Distribution<Eigen::MatrixXd> dist)
{
    int nMeas = dist.get_values().size();
    int nRows = dist.get_value(0).rows();
    int nCols = dist.get_value(0).cols();

    std::vector<std::vector<std::vector<double>>> data(nRows,std::vector<std::vector<double>>(nCols,std::vector<double>(nMeas)));
    std::vector<std::vector<Distribution<double>>> result(nRows,std::vector<Distribution<double>>(nCols));

    for( int i=0;i<nRows;i++)
    for( int j=0;j<nCols;j++)
    {
        for( int n=0;n<dist.get_values().size();n++ )
        {
            data[i][j][n]=dist.get_value(n)(i,j);
        }
        result[i][j] = Distribution<double>(data[i][j],dist.get_resamplingType());
    }
    return result;
}

///////////////////jackknifing for vector of distributions///////////////////////
template <typename T>
std::vector<Distribution<T>> get_vector_resample(std::vector<Distribution<T>> vec_dist, std::string resampling, int nboot=0)
{
    int count = 0;
    std::vector<Distribution<T>> resample;
    for ( auto elem : vec_dist )
    {
        if(resampling == "jackknife")
        {
            resample.push_back(elem.jackknife());
        }
        else if(resampling == "bootstrap")
        {
            resample.push_back(elem.bootstrap(nboot));
        }
            count++;
    }
    return resample;
}

//////////////////////////////wrapper for eigen inversion//////////////////////////////
auto invert(Distribution<Eigen::MatrixXd> dist)
{
    std::vector<Eigen::MatrixXd> result;
    for( auto elem : dist.get_values() )
    {
        result.push_back(elem.inverse());
    }
    return Distribution<Eigen::MatrixXd>(result);
}

#endif
