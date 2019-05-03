#include "AnalysisNPR.h"
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <sstream>
#include <tuple>
#include <fstream>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
//              Code to fit data to a function and extrapolate to user provided values
//
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////



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
    std::vector<int>            latt_size       = parseParam<std::vector<int>>(reader,"latt_size");
    double                      ainv            = parseParam<double>(reader,"ainv");  
    std::string                 dirName         = parseParam<std::string>(reader,"data_directory");
    std::string                 fileName        = parseParam<std::string>(reader,"file_name");
    std::string                 vertex          = parseParam<std::string>(reader,"vertex");
    std::string                 resampling      = parseParam<std::string>(reader,"resamplingType");

    // Read data in numerical form
    std::vector<int>            momentum        = parseParam<std::vector<int>>(reader,"momentum_list");
    std::vector<double>         twist           = parseParam<std::vector<double>>(reader,"twist_list");
    double                      mass            = parseParam<double>(reader,"mass");  
    
    // Read data in label form
    std::vector<std::string>    momentum_label  = parseParam<std::vector<std::string>>(reader,"momentum_list");
    std::vector<std::string>    twist_label     = parseParam<std::vector<std::string>>(reader,"twist_list");
    std::string                 mass_label      = parseParam<std::string>(reader,"mass");  

    // read value of p to extrpolate to and output dir   
    std::string                 fitfunction     = parseParam<std::string>(reader,"fitfunction");                
    std::vector<double>         p_extrap        = parseParam<std::vector<double>>(reader,"p_extrap");
    std::string                 output_dir      = parseParam<std::string>(reader,"output_dir");                
    std::vector<double>         p_range         = parseParam<std::vector<double>>(reader,"p_range");
    ///////////////////////////////////////////////////////////////////////////////////

    
    ////////////////////////// Set up the Distribution of data sets for fitting ////////////////////////////
    std::vector<double>                 p, sigma;
    std::vector<std::vector<double>>    y_vector;
    std::string                         inputFileName;
    std::string                         outputTextFileName = output_dir+"/"+fileName+".txt";
    std::ofstream                            outputTextFile;

    outputTextFile.open (outputTextFileName);
    outputTextFile << "data" << std::endl;

    // calculate momenta and read all amputated vertex function amplitudes
    for(int i=0;i<momentum.size();i++)
    {
        // get the momentum
        double p_tmp = (Grid::sqrt(2)*2.0*M_PI/(static_cast<double>(latt_size[0])))*(static_cast<double>(momentum[i]) + twist[i])*ainv;
        
        std::cout << p_tmp << std::endl;
        //get the inputfile Name to be read
        inputFileName=dirName+"/Lambda_m"+(mass_label)+"_m"+(mass_label);
        inputFileName+="_p0"+(momentum_label[i])+(momentum_label[i])+"0_p";
        inputFileName+=(momentum_label[i])+(momentum_label[i])+"00_tw"+(twist_label[i]);
        
        // read vertex
        Hdf5Reader      h5reader(inputFileName+"/"+fileName+".h5");
        std::vector<double>         data_tmp;
        read(h5reader,vertex,data_tmp);

        // set up the vector of jackknifes 
        Distribution<double>    dist(data_tmp,resampling);
        // write the data to the text file for plotting
        outputTextFile << p_tmp << "\t" << dist.get_central() << "\t" << dist.get_std() << std::endl;


        // only use in interp if in right range
        if( p_tmp > p_range[0] && p_tmp < p_range[1] )
        {
            p.push_back(p_tmp);
            y_vector.push_back(dist.get_values());
        }
    }


    // convert vector of jackknifes to jackknife distribution of vectors
    y_vector = transpose(y_vector);
    Distribution<std::vector<double>> y(y_vector,resampling);


    ////////////////////////// Set up the jackknife fitter ////////////////////// 
    JackKnifeFitter fit(y,p);
    
    if(fitfunction == "inv_p2"){            fit.assignFitFunction(psq_inv_f, psq_inv_df, psq_inv, std::vector<double>(2,1.0)); }
    else if(fitfunction == "inv_p6"){       fit.assignFitFunction(psix_inv_f, psix_inv_df, psix_inv,std::vector<double>(2,1.0));}
    else if(fitfunction == "inv_p2_inv_p6"){fit.assignFitFunction(psq_psix_inv_f, psq_psix_inv_df, psq_psix_inv,std::vector<double>(3,1.0));}
    else if(fitfunction == "p2_inv_p2"){fit.assignFitFunction(psq_psq_inv_f, psq_psq_inv_df, psq_psq_inv,std::vector<double>(3,1.0));}
    else if(fitfunction == "p2"){           fit.assignFitFunction(psq_f,psq_df,psq,std::vector<double>(2,1.0));}
    else if(fitfunction == "p6"){           fit.assignFitFunction(psix_f,psix_df,psix,std::vector<double>(2,1.0));}
    else if(fitfunction == "p2_p6"){        fit.assignFitFunction(psq_psix_f,psq_psix_df,psq_psix,std::vector<double>(3,1.0));}
    else if(fitfunction == "polysq"){        fit.assignFitFunction(polysq_f,polysq_df,polysq,std::vector<double>(4,1.0));}
    else if(fitfunction == "inv_polysq"){        fit.assignFitFunction(inv_polysq_f,inv_polysq_df,inv_polysq,std::vector<double>(4,1.0));}
    else
    {
        std::cout << "Error: fit function must be one of: inv_p6, inv_p2, inv_p2_inv_p6, p2, p6, p2_p6, polysq, inv_polysq" << std::endl;
        return -1;
    }

    
    
    fit.fitAll();
    Distribution<std::vector<double>>    params  = Distribution<std::vector<double>>(fit.get_params(),resampling);
    Distribution<double>                 chi     = Distribution<double>(fit.get_chi(),resampling);
    
    std::cout << "Chi^2 / d.o.f. = " <<  chi.get_central() << " +/- " << chi.get_std() << std::endl;
    std::cout << "The parameters are : " << std::endl;
    for ( int i=0; i<params.get_central().size(); i++)
    {
        std::cout << "p" << i << " = " << params.get_central()[i] << " +/- " << params.get_std()[i] << std::endl;
    } 

    outputTextFile << "extrapolated" << std::endl;
    for(int i=0;i<p_extrap.size();i++)
    { 
        Distribution<double>    f_extrap(fit.extrapolate(p_extrap[i]),resampling);
        outputTextFile << p_extrap[i] << "\t" << f_extrap.get_central() << "\t" << f_extrap.get_std() << std::endl; 
        std::cout << "value at " << p_extrap[i] << "GeV : " << f_extrap.get_central() << " +/- " << f_extrap.get_std() << std::endl; 
        // hdf5 of extrapolated values
        save_result<std::vector<double>>(output_dir+"/"+fileName+"_p"+std::to_string(p_extrap[i])+"GeV.h5",vertex,f_extrap.get_values());
    }


    ///////////////////Writing///////////////////
    // Let's write the value 
    // Write an hdf5 of the jackknife values
    // Then we just want a txt file with
    
    //write the fit results to the text file
    outputTextFile << "fit" << std::endl;
    // find the min and max momentum
     
    std::vector<double>::iterator pmin_ptr = std::min_element(p.begin(), p.end());
    std::vector<double>::iterator pmax_ptr = std::max_element(p.begin(), p.end());
    double p_min = *pmin_ptr;
    double p_max = *pmax_ptr;
    std::cout << "p_min" << p_min << std::endl;
    // if desired range outside of p calculated update the values
    p_min = (p_min < p_range[0]) ? p_min : p_range[0];
    p_max = (p_max > p_range[1]) ? p_max: p_range[1];
    double p_plot;
    Distribution<double>    f_plot;
    for(int i=0; i<1001; i++)
    {
        p_plot = static_cast<float>(i)/1000*(p_max-p_min)+p_min;
        f_plot = Distribution<double>(fit.extrapolate(p_plot),resampling);
        outputTextFile << p_plot << "\t" << f_plot.get_central() << "\t" << f_plot.get_std() << std::endl;
    }


    
    return 0;
}
