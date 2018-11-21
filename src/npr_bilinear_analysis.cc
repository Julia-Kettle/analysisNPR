#include <iostream>
#include <string>
#include <math.h>
#include <sstream>
#include <tuple>
#include "Grid/Grid.h"
#include "IO.h"
#include "distribution.h"
#include "amputation.h"
#include "arithmetic.h"

/////////////////////////////////////////
//  Jackknifing function
//
//     nconf * vertex functions
//      nconf - 1  resamples
//
//  Bootstrapping function
//  So can get averages + std dev + boots. 
//  Need to consider the extrapolation too
//  Write test functions here then organise code mode sensibly
//  File with resampling classes
//  Extrapolation class
//
//  Use Grid as library for the data structures. 
//
///////////////////////////////////////

//Let's write a code for reading in a vector of props
//




// A Distribuion Class?
// average, high, low
//
// Plus a resample class which bootstraps or jackknifes
// central distr, resampled distributions
//

using namespace Grid;
using namespace QCD;

int main()
{
    //////////////////////// Read parameter info from xml //////////////////////////////////
    std::string parameterFileName="../test.xml";
    Grid::XmlReader reader(parameterFileName);

    // read all the inputs 
    // later - check if way to loop through. But then we'd have to use only strings
    // But otherwise might get super lengthy code if many input params
    std::vector<int> latt_size = parseParam<std::vector<int>>(reader,"latt_size");
    std::vector<int> momentum  = parseParam<std::vector<int>>(reader,"momentum");
    int conf_start             = parseParam<int>(reader,"conf_start");
    int conf_inc               = parseParam<int>(reader,"conf_inc");
    int conf_end               = parseParam<int>(reader,"conf_end");
    std::string prop1_file     = parseParam<std::string>(reader,"prop1_file");
    std::string prop2_file     = parseParam<std::string>(reader,"prop2_file");
    std::string vertex_file    = parseParam<std::string>(reader,"vertex_file");
    
    // write to template in case failure
    // Need to put in some error handling here
    std::vector<std::string> par_list = {"latt_size","momentum","conf_start", "conf_inc","conf_end","prop1_file","prop2_file","vertex_file"};
    Grid::XmlWriter writer("template.xml");
    
    for ( auto par_name : par_list )
    {
        write(writer,par_name,"");
    }
    ///////////////////////////////////////////////////////////////////////////////////////////

    // data types for vertex and props //
    std::vector<SpinColourMatrix>               propin,propout;
    std::vector<std::vector<SpinColourMatrix>>  bilin;

    // get configs from start, stop, inc //
    std::vector<int>    configs;
    for(int ic=conf_start; ic < conf_end; ic += conf_inc)
    {
        configs.push_back(ic);
    }
   
    readConfigs(prop1_file, "SinAve", configs, propin);
    readConfigs(prop2_file, "SoutAve", configs, propout);
    readConfigs(vertex_file, "bilinear", configs, bilin);

    // form distribution
    Distribution<SpinColourMatrix>               Sin(propin);
    Distribution<SpinColourMatrix>               Sout(propout);
    std::vector<Distribution<SpinColourMatrix>>  vertex_funcs;
    vertex_funcs = get_vector_distributions(bilin);

    //std::cout << Sin.get_values() << std::endl; 
    //std::cout << Sin.get_mean() << std::endl; 

    //get jackknifes    
    Distribution<SpinColourMatrix>              Sin_jk    =   Sin.jackknife();
    Distribution<SpinColourMatrix>              Sout_jk   =   Sout.jackknife();
    std::vector<Distribution<SpinColourMatrix>> vf_jk     =   get_vector_jackknife(vertex_funcs);


    //amputate the vertices
    auto amp  = amputate(Sout_jk,Sin_jk,vf_jk);
    
    // set gamma indices for projection - S,P,V,A
    std::vector<Gamma::Algebra> I     = {Gamma::Algebra::Identity};
    std::vector<Gamma::Algebra> g5    = {Gamma::Algebra::Gamma5};
    std::vector<Gamma::Algebra> gmu   = {Gamma::Algebra::GammaT,Gamma::Algebra::GammaX,Gamma::Algebra::GammaY,Gamma::Algebra::GammaZ};
    std::vector<Gamma::Algebra> gmug5 = {Gamma::Algebra::GammaTGamma5,Gamma::Algebra::GammaXGamma5,Gamma::Algebra::GammaYGamma5,Gamma::Algebra::GammaZGamma5};

    Distribution<Real> LambdaS = project_gamma(amp,I);
    Distribution<Real> LambdaP = project_gamma(amp,g5);
    Distribution<Real> LambdaV = project_gamma(amp,gmu);
    Distribution<Real> LambdaA = project_gamma(amp,gmug5);
    
    // qslash scheme + projection
    std::vector<double> q(4);
    for (int mu=0;mu<q.size();mu++){ q[mu] = 2*M_PI*momentum[mu]/latt_size[mu]; }

    Distribution<Real> LambdaVq = project_qslash(amp,q,gmu);
    Distribution<Real> LambdaAq = project_qslash(amp,q,gmug5);
    
    std::cout << "gg S" << LambdaS.get_values() << std::endl;
    std::cout << "gg P" << LambdaP.get_values() << std::endl;
    std::cout << "gg V" << LambdaV.get_values() << std::endl;
    std::cout << "gg A" << LambdaA.get_values() << std::endl;
    
    std::cout << "qq V" << LambdaVq.get_values() << std::endl;
    std::cout << "qq A" << LambdaAq.get_values() << std::endl;

    return 0;
}

