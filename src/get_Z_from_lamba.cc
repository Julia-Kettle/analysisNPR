#include <iostream>
#include <string>
#include <math.h>
#include <sstream>
#include <tuple>
#include "Grid/Grid.h"
#include "AnalysisNPR.h"
#include <random>

int main(int argc, char *argv[])
{
    
    //////////////////////// Read parameter info from xml //////////////////////////////////
    std::string parameterFileName;

    std::cout << "starting xml read" << std::endl;
    if (argc <= 1)
    {
        std::cout << "Usage - " << argv[0] << " <input xml filename>" << std::endl;
        // write to template in case failure
        std::vector<std::string> par_list = {"LambdaDir","ZA","ZAerror","scheme","nBoot","outputDir"};
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

    std::cout << "reader set up" << std::endl;
    //read the numerator and denominator infor
   
    std::string LambdaDir   =   parseParam<std::string>(reader,"LambdaDir");
    double      za          =   parseParam<double>(reader,"ZA");
    double      za_e        =   parseParam<double>(reader,"ZAerror");
    std::string scheme      =   parseParam<std::string>(reader,"scheme");
    std::string outputDir   =   parseParam<std::string>(reader,"outputDir");

    std::vector<std::string> vertices, out_vertices;
    if(scheme == "g")
    {
        vertices = {"Sg","Pg","Tg","Vqg"};
        out_vertices = vertices;
    }
    else if(scheme == "q")
    {
        vertices = {"Sg","Pg","Tg","Vq"};
        out_vertices = {"Sq","Pq","Tq","Vq"};
    }
    else
    {
        std::cout << "Error : scheme must be either g or q" << std::endl;
        return -1;
    }

    ////////////////////////////////////////////
    // Read in all the data
    ////////////////////////////////////////////
    Hdf5Reader readA(LambdaDir+"/LambdaA"+scheme+".h5"), readV(LambdaDir+"/LambdaV"+scheme+".h5");
    Hdf5Reader  readT(LambdaDir+"/LambdaTg.h5"), readP(LambdaDir+"/LambdaPg.h5"), readS(LambdaDir+"/LambdaSg.h5");
    std::vector<double> vecA, vecV, vecT, vecS, vecP;
    read(readA,"LambdaA"+scheme,vecA);
    read(readV,"LambdaV"+scheme,vecV);
    read(readS,"LambdaSg",vecS);
    read(readP,"LambdaPg",vecP);
    read(readT,"LambdaTg",vecT);

    ////////////////////////////////////////////
    // Create distributions
    ////////////////////////////////////////////
    Distribution<double> LambdaS(vecS,"bootstrap");
    Distribution<double> LambdaP(vecP,"bootstrap");
    Distribution<double> LambdaT(vecT,"bootstrap");
    Distribution<double> LambdaV(vecV,"bootstrap");
    Distribution<double> LambdaA(vecA,"bootstrap");
   
    int nBoot = LambdaS.get_Nmeas();

    ////////////////////////////////////////////
    // generate bootstraps for ZA
    ////////////////////////////////////////////
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(za,za_e);
    std::vector<double> za_vec;
    for(int i =0 ; i<nBoot; i++)
    {
        za_vec.push_back(distribution(generator));
    }
    za_vec.push_back(za); 
    Distribution<double> ZA(za_vec,"bootstrap");


    ////////////////////////////////////////////
    // Calculate Z from Lambda and ZA
    ////////////////////////////////////////////
    Distribution<double> ZS     =   LambdaA*ZA*(1.0/LambdaS);
    Distribution<double> ZP     =   LambdaA*ZA*(1.0/LambdaP);
    Distribution<double> ZT     =   LambdaA*ZA*(1.0/LambdaT);
    Distribution<double> ZV     =   LambdaA*ZA*(1.0/LambdaV);
    Distribution<double> Zm     =   1.0/ZS;

    ////////////////////////////////////////////
    // Print and write results
    ////////////////////////////////////////////
    std::cout << "ZS = " << ZS.get_central() << " ± " << ZS.get_std() << std::endl;
    std::cout << "ZP = " << ZP.get_central() << " ± " << ZP.get_std() << std::endl;
    std::cout << "ZT = " << ZT.get_central() << " ± " << ZT.get_std() << std::endl;
    std::cout << "ZV = " << ZV.get_central() << " ± " << ZV.get_std() << std::endl;
    std::cout << "Zm = " << Zm.get_central() << " ± " << Zm.get_std() << std::endl;

    save_result<std::vector<double>>(outputDir+"/ZS"+scheme+".h5","ZS"+scheme,ZS.get_values());
    save_result<std::vector<double>>(outputDir+"/Zm"+scheme+".h5","Zm"+scheme,Zm.get_values());
    save_result<std::vector<double>>(outputDir+"/ZP"+scheme+".h5","ZP"+scheme,ZP.get_values());
    save_result<std::vector<double>>(outputDir+"/ZT"+scheme+".h5","ZT"+scheme,ZT.get_values());
    save_result<std::vector<double>>(outputDir+"/ZV"+scheme+".h5","ZV"+scheme,ZV.get_values());
    save_result<std::vector<double>>(outputDir+"/ZA.h5","ZA",ZA.get_values());
    
    
    return 0;
}

