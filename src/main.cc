#include <iostream>
#include <string>
#include <math.h>
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

struct parameters {
    std::string         S1_filename;
    std::string         S2_filename;
    std::string         vertex_filename;
    int                 conf_start;
    int                 conf_end;
    int                 conf_inc;
    std::vector<int>    momentum;
    std::vector<int>    latt_dim;
};

void parseParams()
{
    ////////////////////////////////////////////////
    // let's think about getting params
    // we need: prop1_filename, prop2_filename, vertex_functions ( strings )
    //          config start, inc & end (ints)
    //          mom latt_dim ( vector<ints> )
    //          outputdir (string)
    //
    //

    std::string parameterFileName="../test.xml";
    Grid::XmlReader reader(parameterFileName);
    std::vector<int> output(4);

    
    read(reader,"test", output);

    std::cout << output << std::endl;

    Grid::XmlWriter writer("template.xml");
    write(writer,"test","output");
}


int main()
{
    //LatticeSpinColourMatrix Sin;
    std::vector<SpinColourMatrix>               propin,propout;
    std::vector<std::vector<SpinColourMatrix>>  bilin;
    std::vector<int>                            configs = {4690,4730,4770};

    std::string dirname  = "/tessfs1/home/dp008/dp008/dc-kett1/NPR/GridTest/24cubedTests/vertex_functions/";
    std::string filestem = "../prop_m0.005_p0330";

    //readConfigs(filestem, "SinAve", configs, obj);
    readConfigs(dirname+"props/prop_m0.005_p0330", "SinAve", configs, propin);
    readConfigs(dirname+"props/prop_m0.005_p3300", "SoutAve", configs, propout);
    readConfigs(dirname+"bilinears/Lambda_m0.005_m0.005_p0330_p3300", "bilinear", configs, bilin);

    // form distribution
    Distribution<SpinColourMatrix>               Sin(propin);
    Distribution<SpinColourMatrix>               Sout(propout);
    std::vector<Distribution<SpinColourMatrix>>  vertex_funcs;
    vertex_funcs = get_vector_distributions(bilin); 


    //get jackknifes    
    Distribution<SpinColourMatrix>              Sin_jk    =   Sin.jackknife();
    Distribution<SpinColourMatrix>              Sout_jk   =   Sout.jackknife();
    std::vector<Distribution<SpinColourMatrix>> vf_jk     =   get_vector_jackknife(vertex_funcs);

    //amputate the vertices
    auto amp  = amputate(Sout_jk,Sin_jk,vf_jk);
    
    // get gamma indices for projection
    std::vector<Gamma::Algebra> I     = {Gamma::Algebra::Identity};
    std::vector<Gamma::Algebra> g5    = {Gamma::Algebra::Gamma5};
    std::vector<Gamma::Algebra> gmu   = {Gamma::Algebra::GammaT,Gamma::Algebra::GammaX,Gamma::Algebra::GammaY,Gamma::Algebra::GammaZ};
    std::vector<Gamma::Algebra> gmug5 = {Gamma::Algebra::GammaTGamma5,Gamma::Algebra::GammaXGamma5,Gamma::Algebra::GammaYGamma5,Gamma::Algebra::GammaZGamma5};

    Distribution<ComplexD> LambdaS = project_gamma(amp,I);
    Distribution<ComplexD> LambdaP = project_gamma(amp,g5);
    Distribution<ComplexD> LambdaV = project_gamma(amp,gmu);
    Distribution<ComplexD> LambdaA = project_gamma(amp,gmug5);
    
    // qslash

    std::vector<double> mom = {3.0, 3.0, 0.0, 0.0};
    std::vector<double> q(4);
    std::vector<int>   latt_dim = {24, 24, 24, 64};
    for (int mu=0;mu<q.size();mu++)
    {
        q[mu] = 2*M_PI*mom[mu]/latt_dim[mu];
    }

    Distribution<ComplexD> LambdaVq = project_qslash(amp,q,gmu);
    Distribution<ComplexD> LambdaAq = project_qslash(amp,q,gmug5);
    
    std::cout << LambdaVq.get_mean() << std::endl;
    
    parseParams();

    return 0;
}

