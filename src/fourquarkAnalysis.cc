#include <iostream>
#include <string>
#include <math.h>
#include <sstream>
#include <tuple>
#include "Grid/Grid.h"
#include <Grid/Eigen/Core>
#include "AnalysisNPR.h"

using namespace Grid;
using namespace QCD;

int main(int argc, char *argv[])
{
    ////////////////////////////////////////////////////////////////////////////////////////
    ////                        Read parameter info from xml
    ////////////////////////////////////////////////////////////////////////////////////////
    std::string parameterFileName;

    if (argc <= 1)
    {
        std::cout << "Usage - " << argv[0] << " <input xml filename>" << std::endl;
        // write to template in case failure
        std::vector<std::string> par_list = {"latt_size","momentum1","twist1","momentum2","twist2","conf_start", "conf_inc","conf_end","prop1_file","prop2_file","vertex_file","output_dir"};
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
    std::vector<int> momentum1 = parseParam<std::vector<int>>(reader,"momentum1");
    std::vector<double> twist1 = parseParam<std::vector<double>>(reader,"twist1");
    std::vector<int> momentum2 = parseParam<std::vector<int>>(reader,"momentum2");
    std::vector<double> twist2 = parseParam<std::vector<double>>(reader,"twist2");
    bool qslash_4q              = static_cast<bool>(parseParam<int>(reader,"qslash_4q"));
    std::cout <<"Qslash bool is " << qslash_4q;
    int conf_start             = parseParam<int>(reader,"conf_start");
    int conf_inc               = parseParam<int>(reader,"conf_inc");
    int conf_end               = parseParam<int>(reader,"conf_end");
    int bootstraps              = parseParam<int>(reader,"bootstraps");
    std::string prop1_file     = parseParam<std::string>(reader,"prop1_file");
    std::string prop2_file     = parseParam<std::string>(reader,"prop2_file");
    std::string vertex_file    = parseParam<std::string>(reader,"vertex_file");
    std::string LambdaV_file   = parseParam<std::string>(reader,"LambdaV_file");
    std::string LambdaA_file   = parseParam<std::string>(reader,"LambdaA_file");
    std::string output_dir     = parseParam<std::string>(reader,"output_dir");
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    //      Calculate momenta
    ////////////////////////////////////////////////////////////////////////////////////////
    std::vector<double> p1(Nd),p2(Nd),q(Nd);

    for(int mu=0; mu<Nd; mu++)
    {
        p1[mu] = 2*M_PI*(momentum1[mu]+twist1[mu])/latt_size[mu];
        p2[mu] = 2*M_PI*(momentum2[mu]+twist2[mu])/latt_size[mu];
        q[mu]  = p1[mu]-p2[mu];
    }


    ////////////////////////////////////////////////////////////////////////////////////////
    //      get correct basis for scheme
    ////////////////////////////////////////////////////////////////////////////////////////
    // basis of vertex for tree is always gamma
    std::string                     scheme          = "_g";  
    std::string                     schemeZV;
    std::vector<DiracStructure>     vertex_basis    = gamma_basis();
    std::vector<DiracStructure>     basis           = vertex_basis;
    std::vector<bool>               colourMix(Nop,false);
    if(qslash_4q)
    {
        std::cout << "qslash" << std::endl;
        scheme      = "_q";
        basis       = qslash_basis(p1,p2); 
        colourMix   = {false,false,true,true,false};
    }
    if(LambdaV_file.find("g",0) != std::string::npos && LambdaV_file.find("g",0) != std::string::npos )
    {
        schemeZV = "g";
    }
    else if(LambdaV_file.find("q",0) != std::string::npos && LambdaV_file.find("q",0) != std::string::npos )
    {
        schemeZV = "q";
        std::cout << "V scheme set" << std::endl;
    }
    else
    {
        std::cout << "Error: both LambdaA and LambdaV must have same scheme, either g or q" << std::endl;
        return -1;
    }

    
    ////////////////////////////////////////////////////////////////////////////////////////
    //                  Read in the data for each config  
    ////////////////////////////////////////////////////////////////////////////////////////

    // get vector configs from start, stop, inc //
    std::vector<int>    configs;
    for(int ic=conf_start; ic < conf_end; ic += conf_inc){ configs.push_back(ic);   }
   
    // set up props and vertices
    std::vector<SpinColourMatrix>                           propin,propout;
    std::vector<std::vector<SpinColourSpinColourMatrix>>    fourQ;
    std::vector<double>                                     tmp_lambdaV,tmp_lambdaA;
    
    // Read in the data
    readDataByConfig(prop1_file, "SinAve", configs, propin);
    readDataByConfig(prop2_file, "SoutAve", configs, propout);
    readDataByConfig(vertex_file, "fourquark", configs, fourQ);
    // Read the LambdaV
    Grid::Hdf5Reader reader_V(LambdaV_file);
    read(reader_V,"LambdaV"+schemeZV, tmp_lambdaV);
    // Read the LambdaA
    Grid::Hdf5Reader reader_A(LambdaA_file);
    read(reader_A ,"LambdaA"+schemeZV, tmp_lambdaA);
    std::cout << "Data read" << std::endl;
    

    ////////////////////////////////////////////////////////////////////////////////////////
    //      Set up distributions
    ////////////////////////////////////////////////////////////////////////////////////////
    std::string resampling;
    (bootstraps > 0) ? resampling = "bootstrap" : resampling = "jackknife";
    Distribution<SpinColourMatrix>                          Sin(propin), Sout(propout);
    Distribution<std::vector<SpinColourSpinColourMatrix>>   vertex(fourQ);
    Distribution<double>                                    lambda_v(tmp_lambdaV,resampling);
    Distribution<double>                                    lambda_a(tmp_lambdaA,resampling);
    std::cout << "distributions set up" << std::endl;
       
    if(resampling == "bootstrap")
    {
        Sin     = Sin.bootstrap(bootstraps);
        Sout    = Sout.bootstrap(bootstraps);
        vertex  = vertex.bootstrap(bootstraps);
    }
    else if(resampling == "jackknife")
    {
        Sin     = Sin.jackknife();
        Sout    = Sout.jackknife();
        vertex  = vertex.jackknife();
    }
    std::cout << "dist lenght = " << vertex.get_values().size() <<  " " <<  Sin.get_values().size() << " " << Sout.get_values().size() << std::endl;
    std::cout << trace(vertex.get_values().back()[0]) << std::endl;
    std::cout << trace(Sin.get_values().back()) << std::endl;
    std::cout << trace(Sout.get_values().back()) << std::endl;
    //std::vector<Distribution<SpinColourSpinColourMatrix>> 
    //auto vertex = get_vector_resample(get_vector_distributions(fourQ),resampling,bootstraps);
    

    //////////////////////////////////////////////////////////
    // Perform the projections on tree to get F and invert 
    //////////////////////////////////////////////////////////
    Eigen::MatrixXd tree(Nop,Nop);
    for(int i=0;i<Nop;i++)
    {
        for(int j=0;j<Nop;j++)
        {
            tree(i,j) = projectTree(vertex_basis[i],basis[j],colourMix[j]);
            std::cout << tree(i,j) << "\t";
        }
        std::cout << std::endl;
    }
    Eigen::MatrixXd treeInv =   tree.inverse();
    std::cout << "tree inverted" << std::endl;

    //////////////////////////////////////////////////////////
    // Perform the projections on spin matrix 1 as a debugging check 
    //////////////////////////////////////////////////////////
    SpinColourSpinColourMatrix vertex_one;
    vertex_one = vertex_one + Complex(1,0);
    
    
    
    std::vector<SpinColourSpinColourMatrix> vertex_rho(Grid::QCD::Gamma::nGamma,vertex_one);
    for(int i=0; i<Grid::QCD::Gamma::nGamma; i++)
    {
        SpinColourMatrix  rho;
        rho = rho + Complex(1,0);
        rho = rho*Gamma(i);
        for(int si=0; si < Ns; ++si){
        for(int sj=0; sj < Ns; ++sj){
            for (int ci=0; ci < Nc; ++ci){
            for (int cj=0; cj < Nc; ++cj){
              vertex_rho[i]()(si,sj)(ci,cj)=rho()(si,sj)(ci,cj)*rho();
            }}
        }}
    }

    SpinColourMatrix  rho;
    rho = rho + Complex(1,0);
    auto debug_lambda = projectFourQuark(rho,rho,vertex_rho,vertex_basis,basis,colourMix);
    std::cout << debug_lambda << std::endl;

    
    //////////////////////////////////////////////////////////
    // Perform the projections on vertex data 
    //////////////////////////////////////////////////////////
    //Eigen::MatrixXd                 tmp(Nop,Nop);
    //Distribution<Eigen::MatrixXd>   lambda(std::vector<Eigen::MatrixXd>(configs.size(),tmp));
    Distribution<Eigen::MatrixXd>   lambda  = projectFourQuark(Sin,Sout,vertex,vertex_basis,basis,colourMix);
    std::cout << "vertex projected" << std::endl;
    std::cout << lambda.get_value(0) << std::endl;
    
    std::cout << "length lambda = " << lambda.get_values().size() << std::endl;

    //////////////////////////////////////////////////////////
    // Normalise and jackknife the projected vertices
    //////////////////////////////////////////////////////////
    Distribution<Eigen::MatrixXd> lambda_norm = lambda*treeInv;
    std::cout << lambda_norm.get_mean() << std::endl;
    std::cout << "normalisation and jk done" << std::endl;
    std::cout << "length lambda_norm = " << lambda_norm.get_values().size() << std::endl;

    //////////////////////////////////////////////////////////
    // Divide by Lambda_(A/V) 
    //////////////////////////////////////////////////////////
    Distribution<Eigen::MatrixXd> lambda_ij_vsq; 
    Distribution<Eigen::MatrixXd> lambda_ij_asq;
    std::vector<Eigen::MatrixXd> tmp_v, tmp_a;
    for(int i=0; i<lambda_norm.get_values().size();i++)
    {
        tmp_v.push_back(lambda_norm.get_value(i)*(1/pow(lambda_v.get_value(i),2)));
        tmp_a.push_back(lambda_norm.get_value(i)*(1/pow(lambda_a.get_value(i),2)));
    }
    lambda_ij_vsq = Distribution<Eigen::MatrixXd>(tmp_v,resampling);
    lambda_ij_asq = Distribution<Eigen::MatrixXd>(tmp_a,resampling);
   
    std::cout << lambda_v.get_mean() <<  "    " << lambda_a.get_mean() << std::endl;

    //////////////////////////////////////////////////////////
    // inverte Lambda_ij/Lambda_(A/V) to get Zij/Z(v/a) 
    //////////////////////////////////////////////////////////
    Distribution<Eigen::MatrixXd> Zij_Zvsq = invert(lambda_ij_vsq);
    Distribution<Eigen::MatrixXd> Zij_Zasq = invert(lambda_ij_asq);
    std::cout << Zij_Zasq.get_mean() << std::endl;
    std::cout << Zij_Zvsq.get_mean() << std::endl;
    std::cout << 0.5*(Zij_Zasq.get_mean() + Zij_Zvsq.get_mean()) << std::endl;

    //////////////////////////////////////////////////////////
    // restructure data from dist<matrix> -> matrix<dist> for writing 
    //////////////////////////////////////////////////////////
    auto lambdaNorm_matrix = get_matrix_distributions(lambda_norm);
    auto lambdaNorm_v_matrix = get_matrix_distributions(lambda_ij_vsq);
    auto lambdaNorm_a_matrix = get_matrix_distributions(lambda_ij_asq);
    auto Zij_a_matrix = get_matrix_distributions(Zij_Zvsq);
    auto Zij_v_matrix = get_matrix_distributions(Zij_Zasq);
    
    //////////////////////////////////////////////////////////
    // average of Za and Zv results 
    //////////////////////////////////////////////////////////
    auto lambdaNorm_av_matrix = 0.5*(lambdaNorm_a_matrix + lambdaNorm_v_matrix);
    auto Zij_av_matrix = 0.5*(Zij_a_matrix+Zij_v_matrix);
    
    //////////////////////////////////////////////////////////
    // write the results to file
    /////////////////////////////////////////////////////////
   
    for(int i=0;i<lambdaNorm_matrix.size();i++)
    for(int j=0;j<lambdaNorm_matrix[0].size();j++)
    {
        save_result<std::vector<double>>(output_dir+"/Lambda"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+".h5","Lambda"+std::to_string(i)+std::to_string(j)+scheme+schemeZV,lambdaNorm_matrix[i][j].get_values());
        save_result<std::vector<double>>(output_dir+"/Lambda"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_Vsq.h5","Lambda"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_Vsq",lambdaNorm_v_matrix[i][j].get_values());
        save_result<std::vector<double>>(output_dir+"/Lambda"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_Asq.h5","Lambda"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_Asq",lambdaNorm_a_matrix[i][j].get_values());
        save_result<std::vector<double>>(output_dir+"/Lambda"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_aveVAsq.h5","Lambda"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_aveVAsq",lambdaNorm_av_matrix[i][j].get_values());
        save_result<std::vector<double>>(output_dir+"/Z"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_Vsq.h5","Z"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_Vsq",Zij_v_matrix[i][j].get_values());
        save_result<std::vector<double>>(output_dir+"/Z"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_Asq.h5","Z"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_Asq",Zij_a_matrix[i][j].get_values());
        save_result<std::vector<double>>(output_dir+"/Z"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_aveVAsq.h5","Z"+std::to_string(i)+std::to_string(j)+scheme+schemeZV+"_aveVAsq",Zij_av_matrix[i][j].get_values());


    }

     


    return 0;
}
