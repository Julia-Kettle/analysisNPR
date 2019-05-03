#include <iostream>
#include <string>
#include <math.h>
#include <sstream>
#include <tuple>
#include "Grid/Grid.h"
#include "AnalysisNPR.h"


/*
std::string bilinearTable(std::string inputFileName)
{
    //std::vector<std::string>    vertices    =   {"Sg","Pg","Vg","Ag","Vq","Aq"};
    std::vector<double>         data;
    std::string table = "";

    for (auto vertex : vertices)
    {
        Hdf5Reader      h5reader(inputFileName+"/Lambda"+vertex+".h5");
        read(h5reader,"Lambda"+vertex,data);
        Distribution<double>    dist(data,"jackknife");
        table+=std::to_string(dist.get_mean())+"\t"+std::to_string(dist.get_std())+"\t";
    }
    return table;
}
*/

std::string getTable(std::string inputFileName,std::string prefix, std::vector<std::string> vertices, std::string suffix)
{
    std::vector<double>         data;
    std::string table = "";

    for (auto vertex : vertices)
    {
        Hdf5Reader      h5reader(inputFileName+"/"+prefix+vertex+suffix + ".h5");
        read(h5reader,prefix+vertex,data);
        Distribution<double>    dist(data,"jackknife");
        table+=std::to_string(dist.get_mean())+"\t"+std::to_string(dist.get_std())+"\t";
    }
    return table;
}

int main(int argc, char *argv[])
{
    //////////////////////// Read parameter info from xml //////////////////////////////////
    std::string parameterFileName;

    if (argc <= 1)
    {
        std::cout << "Usage - " << argv[0] << " <input xml filename>" << std::endl;
        // write to template in case failure
        std::vector<std::string> par_list = {"latt_size","ainv","mass","momentum_list","twist_list","data_directory"};
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
    std::vector<int>    latt_size   = parseParam<std::vector<int>>(reader,"latt_size");
    double              ainv        = parseParam<double>(reader,"ainv");  
    std::string         dirName     = parseParam<std::string>(reader,"data_directory");
    std::string         data_type   = parseParam<std::string>(reader,"data_type");
    std::string         data_suffix = parseParam<std::string>(reader,"data_suffix");
    std::string         data_prefix = parseParam<std::string>(reader,"data_prefix");

    // Read data in numerical form
    std::vector<int>    momentum    = parseParam<std::vector<int>>(reader,"momentum_list");
    std::vector<double> twist       = parseParam<std::vector<double>>(reader,"twist_list");
    double              mass        = parseParam<double>(reader,"mass");  
    
    // Read data in label form
    std::vector<std::string>    momentum_label    = parseParam<std::vector<std::string>>(reader,"momentum_list");
    std::vector<std::string>    twist_label       = parseParam<std::vector<std::string>>(reader,"twist_list");
    std::string                 mass_label        = parseParam<std::string>(reader,"mass");  
    
    std::vector<double>                 p(momentum.size());
    std::vector<std::vector<double>>    data(momentum.size());

    std::string     inputFileName;

    if(data_type=="bilinear" || data_type=="fourquark")
    {
        for(int i=0;i<momentum.size();i++)
        {
            // get the momentum
            p[i] = (Grid::sqrt(2)*2.0*M_PI/(static_cast<double>(latt_size[0])))*(static_cast<double>(momentum[i]) + twist[i])*ainv;

            //get the inputfile Name to be read
            inputFileName=dirName+"/Lambda_m"+(mass_label)+"_m"+(mass_label);
            inputFileName+="_p0"+(momentum_label[i])+(momentum_label[i])+"0_p";
            inputFileName+=(momentum_label[i])+(momentum_label[i])+"00_tw"+(twist_label[i]);
            std::vector<std::string> vertices;
            if (data_type == "bilinear")
            {
                vertices  = std::vector<std::string>({"Sg","Pg","Vg","Ag","Vq","Aq"});
            }
            else
            {
                for(int j=0;j<5;j++)
                for(int k=0;k<5;k++)
                {
                    vertices.push_back(std::to_string(j)+std::to_string(k));
                }
            }
            std::cout << p[i] << "\t" << getTable(inputFileName,data_prefix,vertices,data_suffix) << std::endl;
        }
    }
    else{ return -1; }
    
    return 0;

}
