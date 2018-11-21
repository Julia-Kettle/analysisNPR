#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "Grid/Grid.h"
#include "arithmetic.h"
#include "utils.h"

/* Distribution.h
 * Julia Kettle Nov 2018
 * Defines a distribution class for holding statistical distribution of various data types
 * Properties: 
 *  - values : vector<T> 
 *  - mean   : T
 *  - Nmeas  : T
 *  These are protected and cannot be set other than in constructor
 *
 *  Methods:
 *  - jackknife : returns the means of jackknifed resamples
 *
 *  Operator overloading: +,-,,* for Distribution op {Distribution or T or some type V}
 *   requires T*T and T*V to be defined 
*/

/////////////////// Distribution Class ////////////////////
template<class T>
class Distribution
{
    private:
        // Properties
        std::vector<T>  values;
        size_t          Nmeas;
        T               mean;

    public:
        // Constructors
        Distribution(){}; 
        Distribution(std::vector<T> values);
        
        // return functions
        std::vector<T>  get_values(){ return this->values; }
        T               get_value(int index){ return this->values[index]; }
        T               get_mean(){ return this->mean; }
        size_t          get_Nmeas(){ return this->Nmeas; }

        // resampling
        Distribution<T> jackknife();

        //operator overloading 
        //Dist<T> op Dist<T>
        Distribution <T> operator + (Distribution<T> obj2){ return Distribution(this->values+obj2.values); }
        Distribution <T> operator - (Distribution<T> obj2){ return Distribution(this->values-obj2.values); }
        Distribution <T> operator * (Distribution<T> obj2){ return Distribution(this->values*obj2.values); }
        //Dist<T> op V 
        template <typename V>
        auto operator + (V obj2){ return Distribution(this->values + obj2); }
        template <typename V>
        auto operator - (V obj2){ return Distribution(this->values - obj2); }
        template <typename V>
        auto operator * (V obj2){ return Distribution(this->values * obj2); }
        
};

//////////////////////Constructor/////////////////////
template<class T>
Distribution<T>::Distribution(std::vector<T> dist)
{
    values  = dist;
    Nmeas   = values.size();

    mean = zero(values);  // Need to make this zero
    for(auto value : values){ mean=mean+value*(1.0/values.size()); }
}


////////////////jackknife resamples/////////////////////
// Sample = value[i] ; i=0->N
// Resample[j] = { values[i] } i!=j
// -> N resamples of distributions size N-1
// Save only the means to the distribution
///////////////////////////////////////////////////////
template<class T>
Distribution<T> Distribution<T>::jackknife()
{
    std::vector<T> resampled_means;
    resampled_means.resize(Nmeas+1);

    for (int i=0;i<Nmeas;i++)
    {
        resampled_means[i] = mean - values[i];
        resampled_means[i] = resampled_means[i]*(double(Nmeas)/double(Nmeas-1));
    }
    resampled_means.back() = mean;
    return Distribution<T>(resampled_means);
    //return jk;
}

/* Non class functions 
 * some operator overloading
 * some methods to create vectors of Distributions
*/ 


/////////////////////////////// Zeros for different data types ///////////////////////////////////
Grid::QCD::ComplexD zero(std::vector<Grid::QCD::ComplexD> values){ return Grid::QCD::ComplexD(0.0,0.0); }

template <typename T>
T zero(std::vector<T> values){ return Grid::QCD::zero; }

template <typename T>
std::vector<T> zero(std::vector<std::vector<T>> values){ return std::vector<T>(values[0].size(),Grid::QCD::zero); }


/////////////////////////////////// Operator Overloading ( left-wise ) /////////////////////////////////////////

// V op Distribution<T>
template<typename T, typename V>
Distribution<T> operator * (V lhs, Distribution<T> rhs){ return Distribution<T>(lhs*rhs.get_values()); }

template<typename T, typename V>
Distribution<T> operator + (V lhs, Distribution<T> rhs){ return Distribution<T>(lhs+rhs.get_values()); }

template<typename T, typename V>
Distribution<T> operator - (V lhs, Distribution<T> rhs){ return Distribution<T>(lhs+rhs.get_values()); }

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

template <typename T>
std::vector<Distribution<T>> get_vector_jackknife(std::vector<Distribution<T>> vec_dist)
{
    int count = 0;
    std::vector<Distribution<T>> jackknife;
    for ( auto elem : vec_dist )
    {
        jackknife.push_back(elem.jackknife());
        count++;
    }
    return jackknife;
}

//overload functions with distribution as argument

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


#endif
