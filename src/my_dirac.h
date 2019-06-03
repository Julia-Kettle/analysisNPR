#ifndef MY_DIRAC_H
#define MY_DIRAC_H

#include <AnalysisNPR.h>
    
////////////////////////////////////////////////////////////////////////////////////////
//      Set up gamma structures and bases  
////////////////////////////////////////////////////////////////////////////////////////

std::vector<DiracStructure> gamma_basis()
{
    /////////////// gamma basis ////////////////
    //VVpAA

    // dirac structures for full basis
    DiracStructure  VVpAA;
    DiracStructure  VVmAA;
    DiracStructure  SSmPP;
    DiracStructure  SSpPP;
    DiracStructure  TT;
    
    SpinMatrix rho;
    rho = rho+Complex(1.0,0);
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
    return basis;
}

std::vector<DiracStructure> qslash_basis(std::vector<double> p1, std::vector<double> p2)
{
    
    std::vector<DiracStructure>     gg_basis = gamma_basis();

    DiracStructure VVpAA = gg_basis[0];
    DiracStructure VVmAA = gg_basis[1];
    DiracStructure TT    = gg_basis[4];
    
    //qq
    DiracStructure  VVpAAq;
    DiracStructure  VVmAAq;
    DiracStructure  TTq;

    SpinMatrix rho;
    rho = rho+Complex(1.0,0);

    /////////////////// qslash basis //////////////////
    // VVpAA
    std::vector<double> q  = p1-p2;

    std::cout << p1 << std::endl;
    std::cout << p2 << std::endl;
    std::cout << q << std::endl;

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

    //qslash
    SpinMatrix qslash;
    for(int i=0;i<4;i++)
    {
        SpinMatrix Grho;
        Grho = Grho+ Complex(1,0);
        qslash = qslash + (Grho*Grid::QCD::Gamma(gmu[i]))*q[i];
    }
    
    int count=0;
    SpinMatrix psigp;
    for(int i=0;i<Nd;i++)
    {
        for(int j=i+1;j<Nd;j++)
        {
        SpinMatrix Grho;
        Grho = Grho+ Complex(1,0);
        psigp = psigp + TT.gammas[count]*0.5*(p1[j]*p2[i]/Grid::sqrt(p1sq*p2sq-pow(p1dotp2,2)));
        psigp = psigp - TT.gammas[count]*Grid::QCD::Gamma(g5[0])*0.5*(p1[j]*p2[i]/Grid::sqrt(p1sq*p2sq-pow(p1dotp2,2)));
        count++;
        }
    }

    VVpAAq.signs    = VVpAA.signs;
    for(int i=0;i<VVpAA.gammas.size();i++)
    { 
        if(i>=4)
        {
            VVpAAq.gammas.push_back(qslash*(1.0/Grid::sqrt(qsq)));
        }
        else
        {
            VVpAAq.gammas.push_back(qslash*Grid::QCD::Gamma(g5[0])*(1.0/Grid::sqrt(qsq)));
        }
    }

    // VVmAA
    VVmAAq = VVpAAq;
    VVmAAq.signs = VVmAA.signs;


    //TT
    TTq=TT;
    // j < i
    count = 0;
    for(int i=0;i<Nd;i++)
        for(int j=i+1;j<Nd;j++)
    {
        //TTq.gammas[count] = TTq.gammas[count]*0.5*(p1[j]*p2[i]/Grid::sqrt(p1sq*p2sq-pow(p1dotp2,2)));
        //TTq.gammas[count] = TTq.gammas[count] - TTq.gammas[count]*Grid::QCD::Gamma(g5[0]);
        TTq.gammas[count] = psigp;
        count++;
    }


    std::vector<DiracStructure> basis = {VVpAAq,VVmAAq,VVmAAq,TTq,TTq};
    return basis;
}
#endif
