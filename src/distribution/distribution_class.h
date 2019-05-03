#ifndef DISTRIBUTION_CLASS_H
#define DISTRIBUTION_CLASS_H

#include <algorithm>
#include "Grid/Grid.h"
#include <Grid/Eigen/Core>
#include "maths/maths.h"
#include <cmath>
#include <random>
#include <ctime>

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
 *  - bootstrap : returns randomly resample distribution
 *
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
        std::string     resampling; // "jackknife", "bootstrap" or "none"

    public:
        // Constructors
        Distribution(){}; 
        Distribution(std::vector<T> values);
        Distribution(std::vector<T> values, std::string resampling);
        
        // return functions
        std::vector<T>  get_values(){ return this->values; }
        T               get_value(int index){ return this->values[index]; }
        T               get_mean(){ return this->mean; }
        size_t          get_Nmeas(){ return this->Nmeas; }
        std::string     get_resamplingType(){ return this->resampling; }
        T               get_central();
        // These may not work for certain types
        T               get_std();

        // resampling
        Distribution<T> jackknife();
        Distribution<T> bootstrap(int Nboot);

        //operator overloading 
        //Dist<T> op Dist<T>
        Distribution <T> operator + (Distribution<T> obj2){ return Distribution(this->values+obj2.values,this->resampling); }
        Distribution <T> operator - (Distribution<T> obj2){ return Distribution(this->values-obj2.values,this->resampling); }
        Distribution <T> operator * (Distribution<T> obj2){ return Distribution(this->values*obj2.values,this->resampling); }
        //Dist<T> op V 
        template <typename V>
        auto operator + (V obj2){ return Distribution(this->values + obj2,this->resampling); }
        template <typename V>
        auto operator - (V obj2){ return Distribution(this->values - obj2,this->resampling); }
        template <typename V>
        auto operator * (V obj2){ return Distribution(this->values * obj2,this->resampling); }
        
};

/////////////////////////////////////////Constructors///////////////////////////////////////////
template<class T>
Distribution<T>::Distribution(std::vector<T> dist)
{
    values  = dist;
    Nmeas   = values.size();
    mean = zero(values[0]);  // Need to make this zero
    for(auto value : values){ mean=mean+value*(1.0/values.size()); }
    resampling = "none";
}

template<class T>
Distribution<T>::Distribution(std::vector<T> dist, std::string sampleType)
{
    values  = dist;
    resampling = sampleType;
    // if resampled then the last value of distribution is the central value - not really part of the dist. 
    (resampling == "none") ? Nmeas   = values.size() : Nmeas = values.size()-1; 
    mean = zero(values[0]);  // Need to make this zero
    for(int i=0; i<Nmeas; i++){ mean=mean+values[i]*(1.0/Nmeas); }
}


/////////////////////////////////////////Resampling//////////////////////////////////////////////

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
        resampled_means[i] = mean*(double(Nmeas)) - values[i];
        resampled_means[i] = resampled_means[i]*(1/double(Nmeas-1));
    }
    resampled_means.back() = mean;
    return Distribution<T>(resampled_means,"jackknife");
    //return jk;
}

////////////////////bootstrap resamples////////////////////
// Sample = values[i] ; i=0->N
// Resample[j] = { values[k] } ; values[k] = values[rand(i=0->N)] ; k=0->Nboot
// Nboot not neccesarily == N
// save the means of each distribution to form new distribution
///////////////////////////////////////////////////////
template <class T>
Distribution<T> Distribution<T>::bootstrap(int Nboot)
{
    std::default_random_engine generator;
    generator.seed(time(NULL));
    std::uniform_int_distribution<int> uniform(0,this->Nmeas-1);
    std::vector<T> bootstrap_means;
    for ( int iboot=0; iboot<Nboot; iboot++ )
    {
        T mean_b = zero(this->values[0]);
        for (int i=0;i<this->Nmeas;i++)
        {
            int index = uniform(generator);
            T sample = values[index];
            mean_b = mean_b + sample;
        }
        bootstrap_means.push_back(mean_b*(1.0/this->Nmeas));
    }
    bootstrap_means.push_back(this->mean);

    return Distribution<T>(bootstrap_means,"bootstrap");
}

//////////////////////////////////////////////////stats/////////////////////////////////////////////////
template <class T>
T Distribution<T>::get_central()
{
    // DEfines a central value - last value if resampled where we store central
    // otherwise take the mean
    T central;
    (resampling == "none") ? central = this->mean : central = this->get_values().back();
    return central;
}

template <class T>
T Distribution<T>::get_std()
{
    int factor;
    (resampling == "jackknife") ? factor = this->Nmeas - 1 : factor = 1;
    T std;
    std=zero(values[0]);
    for ( auto elem : this->values )
    {
        // technically should skip last element for jk and bs 
        // but should = 0 with this anyway
        std = std + (elem - get_central())*(elem - get_central());
    }
    std  = std*( double(factor) / this->Nmeas);
    std  = sqrt(std);
    return std;
}




#endif

