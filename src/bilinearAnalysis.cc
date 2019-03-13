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

using namespace Grid;
using namespace QCD;

int main(int argc, char *argv[])
{
    //////////////////////// Read parameter info from xml //////////////////////////////////
    std::cout << "Reading parameters from xml" << std::endl;
    std::string parameterFileName;

    if (argc <= 1)
    {
        std::cout << "Usage - " << argv[0] << " <input xml filename>" << std::endl;
        // write to template in case failure
        std::vector<std::string> par_list = {"latt_size","momentum","conf_start", "conf_inc","conf_end","prop1_file","prop2_file","vertex_file","output_dir"};
        Grid::XmlWriter writer("template.xml");
        for ( auto par_name : par_list ){ write(writer,par_name,""); }
        
        // exit with err
        return -1;
    }
    else
    {
        parameterFileName=argv[1];
    }

    // set up xml reader
    Grid::XmlReader reader(parameterFileName);

    // read all the inputs 
    std::vector<int> latt_size = parseParam<std::vector<int>>(reader,"latt_size");
    std::vector<int> momentum  = parseParam<std::vector<int>>(reader,"momentum");
    std::vector<int> twist     = parseParam<std::vector<int>>(reader,"twist");
    int conf_start             = parseParam<int>(reader,"conf_start");
    int conf_inc               = parseParam<int>(reader,"conf_inc");
    int conf_end               = parseParam<int>(reader,"conf_end");
    std::string prop1_file     = parseParam<std::string>(reader,"prop1_file");
    std::string prop2_file     = parseParam<std::string>(reader,"prop2_file");
    std::string vertex_file    = parseParam<std::string>(reader,"vertex_file");
    std::string output_dir     = parseParam<std::string>(reader,"output_dir");
    ///////////////////////////////////////////////////////////////////////////////////////////

    // data types for vertex and props //
    std::vector<SpinColourMatrix>               propin,propout;
    std::vector<std::vector<SpinColourMatrix>>  bilin;

    // get vector configs from start, stop, inc //
    std::vector<int>    configs;
    for(int ic=conf_start; ic < conf_end; ic += conf_inc)
    {
        configs.push_back(ic);
    }
   
    // read in the data for all specified configurations
    readDataByConfig(prop1_file, "SinAve", configs, propin);
    readDataByConfig(prop2_file, "SoutAve", configs, propout);
    readDataByConfig(vertex_file, "bilinear", configs, bilin);

    // form distribution
    Distribution<SpinColourMatrix>               Sin(propin);
    Distribution<SpinColourMatrix>               Sout(propout);
    std::vector<Distribution<SpinColourMatrix>>  vertex_funcs;
    // get vector of distributions
    vertex_funcs = get_vector_distributions(bilin);

    //get jackknifes    
    Distribution<SpinColourMatrix>              Sin_jk    =   Sin.jackknife();
    Distribution<SpinColourMatrix>              Sout_jk   =   Sout.jackknife();
    std::vector<Distribution<SpinColourMatrix>> vf_jk     =   get_vector_jackknife(vertex_funcs);

    //amputate the vertices
    auto amp  = amputate(Sout,Sin,vertex_funcs);
    
    // set gamma indices for projection - S,P,V,A
    std::vector<Gamma::Algebra> I     = {Gamma::Algebra::Identity};
    std::vector<Gamma::Algebra> g5    = {Gamma::Algebra::Gamma5};
    std::vector<Gamma::Algebra> gmu   = {Gamma::Algebra::GammaT,Gamma::Algebra::GammaX,Gamma::Algebra::GammaY,Gamma::Algebra::GammaZ};
    std::vector<Gamma::Algebra> gmug5 = {Gamma::Algebra::GammaTGamma5,Gamma::Algebra::GammaXGamma5,Gamma::Algebra::GammaYGamma5,Gamma::Algebra::GammaZGamma5};

    // Gamma scheme projections
    std::cout << "Projecting" << std::endl;
    Distribution<Real> LambdaS = project_gamma(amp,I);
    Distribution<Real> LambdaP = project_gamma(amp,g5);
    Distribution<Real> LambdaV = project_gamma(amp,gmu);
    Distribution<Real> LambdaA = project_gamma(amp,gmug5);
    
    // qslash scheme + projection
    std::vector<double> q(4);
    for (int mu=0;mu<q.size();mu++){ q[mu] = 2*M_PI*(momentum[mu])/latt_size[mu]; }

    Distribution<Real> LambdaVq = project_qslash(amp,q,gmu);
    Distribution<Real> LambdaAq = project_qslash(amp,q,gmug5);

    // get jackknives
    LambdaS  = LambdaS.jackknife();
    LambdaP  = LambdaP.jackknife();
    LambdaV  = LambdaV.jackknife();
    LambdaA  = LambdaA.jackknife();
    LambdaVq = LambdaVq.jackknife();
    LambdaAq = LambdaAq.jackknife();
    
   

    // Print out values + save to file    
    std::cout << "gg S" << LambdaS.get_values() << std::endl;
    std::cout << LambdaS.get_mean() << " " << LambdaS.get_std() << std::endl;
    save_result<std::vector<double>>(output_dir+"/Sg.h5","LambdaSg",LambdaS.get_values());
    
    std::cout << "gg P" << LambdaP.get_values() << std::endl;
    std::cout << LambdaP.get_mean() << " " << LambdaP.get_std() << std::endl;
    save_result<std::vector<double>>(output_dir+"/Pg.h5","LambdaPg",LambdaP.get_values());
    
    std::cout << "gg V" << LambdaV.get_values() << std::endl;
    std::cout << LambdaV.get_mean() << " " << LambdaV.get_std() << std::endl;
    save_result<std::vector<double>>(output_dir+"/Vg.h5","LambdaVg",LambdaV.get_values());
    
    std::cout << "gg A" << LambdaA.get_values() << std::endl;
    std::cout << LambdaA.get_mean() << " " << LambdaA.get_std() << std::endl;
    save_result<std::vector<double>>(output_dir+"/Ag.h5","LambdaAg",LambdaA.get_values());
    
    
    std::cout << "qq V" << LambdaVq.get_values() << std::endl;
    std::cout << LambdaVq.get_mean() << " " << LambdaVq.get_std() << std::endl;
    save_result<std::vector<double>>(output_dir+"/Vq.h5","LambdaVq",LambdaVq.get_values());
    
    std::cout << "qq A" << LambdaAq.get_values() << std::endl;
    std::cout << LambdaAq.get_mean() << " " << LambdaAq.get_std() << std::endl;
    save_result<std::vector<double>>(output_dir+"/Aq.h5","LambdaAq",LambdaAq.get_values());

    return 0;
}

