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

int main()
{

    std::vector<double> q  = {4,0,-4,0};
    std::vector<double> p1 = {4,4,0,0};
    std::vector<double> p2 = {0,4,4,0};

    // set gamma indices for projection - S,P,V,A,T
    std::vector<Gamma::Algebra> I           = {Gamma::Algebra::Identity};
    std::vector<Gamma::Algebra> g5          = {Gamma::Algebra::Gamma5};
    std::vector<Gamma::Algebra> gmu         = {Gamma::Algebra::GammaT,Gamma::Algebra::GammaX,Gamma::Algebra::GammaY,Gamma::Algebra::GammaZ};
    std::vector<Gamma::Algebra> gmug5       = {Gamma::Algebra::GammaTGamma5,Gamma::Algebra::GammaXGamma5,Gamma::Algebra::GammaYGamma5,Gamma::Algebra::GammaZGamma5};
    std::vector<Gamma::Algebra> sigma_mu_nu = {Gamma::Algebra::SigmaXT,Gamma::Algebra::SigmaXY,Gamma::Algebra::SigmaXZ,
                                                    Gamma::Algebra::SigmaYT,Gamma::Algebra::SigmaYZ,Gamma::Algebra::SigmaZT};
    // dirac structures for full basis
    DiracStructure  VVpAA;
    DiracStructure  VVmAA;
    DiracStructure  SSmPP;
    DiracStructure  SSpPP;
    DiracStructure  TT;
    
    
    DiracStructure  VVpAAq;
    DiracStructure  VVmAAq;
    DiracStructure  TTq;
        
    SpinMatrix rho;
    rho = rho+Complex(1.0,0);
    std::cout << rho << std::endl;

    double qsq     = 0;
    double p1sq    = 0;
    double p2sq    = 0;
    double p1dotp2 = 0;
    for(int i=0;i<4;i++)
    {
        qsq+=q[i]*q[i];
        p1sq+=p1[i]*p1[i];
        p2sq+=p2[i]*p2[i];
        p1dotp2+=p1[i]*p2[i];
    }

    std::cout << q << qsq << Grid::sqrt(qsq) << std::endl;

    /////////////// gamma basis ////////////////
    //VVpAA
    for(int i=0;i<gmu.size();i++){  VVpAA.gammas.push_back(rho*Gamma(gmu[i])); VVpAA.indices.push_back(gmu[i]);   }
    for(int i=0;i<gmug5.size();i++) {  VVpAA.gammas.push_back(rho*Gamma(gmug5[i])); VVpAA.indices.push_back(gmug5[i]);  }
    VVpAA.signs     = {1,1,1,1,1,1,1,1};
    
    //VVmAA
    VVmAA.gammas    = VVpAA.gammas;
    VVmAA.indices    = VVpAA.indices;
    VVmAA.signs     = {1,1,1,1,-1,-1,-1,-1};
    
    //SSpPP
    SSpPP.signs     = {1,1};
    SSpPP.gammas.push_back(rho*Gamma(I[0]));
    SSpPP.gammas.push_back(rho*Gamma(g5[0]));
    SSpPP.indices.push_back(I[0]);
    SSpPP.indices.push_back(g5[0]);
    
    //SSmPP
    SSmPP.signs     = {1,-1};
    SSmPP.gammas    = SSpPP.gammas;
    SSmPP.indices   = SSpPP.indices;
    
    //TT
    TT.signs        =   {1,1,1,1,1,1};
    for(int i=0;i<sigma_mu_nu.size();i++){  TT.gammas.push_back(rho*Gamma(sigma_mu_nu[i])); TT.indices.push_back(sigma_mu_nu[i]); }
    
    std::vector<DiracStructure> basis = {VVpAA,VVmAA,SSmPP,SSpPP,TT};

    for(int i=0;i<basis.size();i++)
    {
        for(int j=0;j<basis.size();j++)
        {
            std::cout << projectTree(basis[i],basis[j],false) << "\t";
        }
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    

    /////////////////// qslash basis //////////////////
    // VVpAA
    VVpAAq.signs    = VVpAA.signs;
    for(int i=0;i<VVpAA.gammas.size();i++)
    { 
        VVpAAq.gammas.push_back(rho*VVpAA.gammas[i]*(q[i%4]/Grid::sqrt(qsq)));
    }

    // VVmAA
    VVmAAq = VVpAAq;
    VVmAAq.signs = VVmAA.signs;


    //TT
    TTq=TT;
    int count=0;
    for(int i=0;i<Nd;i++)
    for(int j=i+1;j<Nd;j++)
    {
        TTq.gammas[count] = TTq.gammas[count]*0.5*(p1[i]*p2[j]/Grid::sqrt(p1sq*p2sq-pow(p1dotp2,2)));
        TTq.gammas[count] = TTq.gammas[count] - TTq.gammas[count]*Grid::QCD::Gamma(g5[0]);
        count++;
    }

    // vector of dirac structures
    std::vector<DiracStructure> basis_q = {VVpAAq,VVmAAq,VVmAAq,TTq,TTq};
    std::vector<bool>               basis_q_mix = {false,false,true,true,false};

    for(int i=0;i<basis.size();i++)
    {
        for(int j=0;j<basis.size();j++)
        {
            std::cout << projectTree(basis[i],basis_q[j],basis_q_mix[j]) << "\t";
        }
        std::cout << std::endl;
    }


}
