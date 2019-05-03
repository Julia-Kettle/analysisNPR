#include <iostream>
#include <string>
#include <math.h>
#include <sstream>
#include <tuple>
#include "Grid/Grid.h"
#include "AnalysisNPR.h"

int main(int argc, char *argv[])
{
    
    //////////////////////// Read parameter info from xml //////////////////////////////////
    std::string parameterFileName;

    std::cout << "starting xml read" << std::endl;
    if (argc <= 1)
    {
        std::cout << "Usage - " << argv[0] << " <input xml filename>" << std::endl;
        // write to template in case failure
        std::vector<std::string> par_list = {"latt_size","mass","momentum_list","twist_list","data_directory"};
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
    // Read data in numerical form
    std::vector<int>    momentum    = parseParam<std::vector<int>>(reader,"momentum_list");
    std::vector<double> twist       = parseParam<std::vector<double>>(reader,"twist_list");
    double              mass        = parseParam<double>(reader,"mass");  
    
    // Read data in label form
    std::vector<std::string>    momentum_label    = parseParam<std::vector<std::string>>(reader,"momentum_list");
    std::vector<std::string>    twist_label       = parseParam<std::vector<std::string>>(reader,"twist_list");
    std::string                 mass_label        = parseParam<std::string>(reader,"mass");  

    std::cout << "mom twist and mass read" << std::endl;
    //read the numerator and denominator infor
   
    std::string numName   =   parseParam<std::string>(reader,"numName");
    std::string numDir      =   parseParam<std::string>(reader,"numDir");
    std::string den1Name   =   parseParam<std::string>(reader,"den1Name");
    std::string den1Dir      =   parseParam<std::string>(reader,"den1Dir");
    std::string den2Name   =   parseParam<std::string>(reader,"den2Name");
    std::string den2Dir      =   parseParam<std::string>(reader,"den2Dir");
    // output
    std::string outputName   =   parseParam<std::string>(reader,"outputName");
    std::string outputDir      =   parseParam<std::string>(reader,"outputDir");
  

    std::cout << "numer and denom read" << std::endl;
    for(int i=0;i<momentum.size();i++)
    {
        std::string subdir = "/Lambda_m"+(mass_label)+"_m"+(mass_label)+"_p0"+(momentum_label[i])+(momentum_label[i])+"0_p";
                    subdir+=(momentum_label[i])+(momentum_label[i])+"00_tw"+(twist_label[i])+"/";
        std::cout << subdir << std::endl;

        std::vector<double> num, den1, den2;

        Hdf5Reader  numReader(numDir+subdir+numName+".h5");
        read(numReader,numName,num);
        
        std::cout << "num read" << std::endl;

        Hdf5Reader  den1Reader(den1Dir+subdir+den1Name+".h5");
        read(den1Reader,den1Name,den1);
        
        std::cout << "den1 read" << std::endl;
    
        if((den1Name == den2Name) and (den1Dir == den2Dir))
        {  
            den2 = den1; 
            std::cout << "den2 copied" << std::endl;
        }
        else
        {
            Hdf5Reader  den2Reader(den2Dir+subdir+den2Name+".h5");
            read(den2Reader,den2Name,den2);
            std::cout << "den2 read" << std::endl;
        }
        std::vector<double> result;
        std::cout << num << std::endl;
        std::cout << den1 << std::endl;
        std::cout << den2 << std::endl;
        std::cout << num/den2 << std::endl;
        result = num/(den1*den2);
        std::cout << result << std::endl;
        save_result<std::vector<double>>(outputDir+subdir+outputName+".h5",outputName,result);

    }

    return 0;
}


