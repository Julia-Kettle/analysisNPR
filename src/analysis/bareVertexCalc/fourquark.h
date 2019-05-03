#ifndef FOURQUARK_H
#define FOURQUARK_H

#include "Grid/Grid.h"
#include <string>
#include <iostream>
#include <Grid/Eigen/Core>
#include <Grid/Eigen/Dense>
#include <Grid/Eigen/SVD>
#include <complex>
#include "amputation.h"
#include "distribution/distribution.h"


using namespace Grid;
using namespace QCD;


struct DiracStructure{
    std::vector<Gamma::Algebra> indices;
    std::vector<SpinMatrix>     gammas;
    std::vector<int>            signs;
};

////////////////////////////////////
//  Four Quark Projections and Amputation
//  parts:
//      g5*adj(Sinv(s1,c1))*g5 * g5*adj(Sinv(s3,c3))*g5 * V((s1,s2),(c1,c2),(s3,s4),(c3,c4)) * Sinv(s2,c2) * Sinv(s4*c4)
//
//      The two Wick contractions are
//      tr [ g5 adj(s(p2;x)) * g5 * Gamma * d(p1;x) ] x tr[ g5 adj(s(p2;x)) * g5 * Gamma * d(p1;x) ] (1)
//      - tr [ g5 adj(s(p2;x)) * g5 * Gamma *  d(p1;x) * g5 adj(s(p2;x)) * g5 * Gamma * d(p1;x) ] (2) 
//
//      notice that the bilinear parts are the same
//
//      V_4q_{ijkl;abcd}    =   V_bil_{ij;ab}*V_bil_{kl;cd}
//                          =   g5*adj(Sout_i,a)*g5*Gamma*Sin_j,b * g5*adj(Sout_k,c)*g5*Gamma*Sin_l,d 
//
//      // figure of eight
//      P[Pi_4q_comp1]      =  [ inv(Sin) * Gamma_proj * g5*adj(inv(Sout))*g5 ]_ji,ba * [ V_4q_{ijkl;abcd} ] * [ inv(Sin) * Gamma_proj * g5*adj(Sout)*g5 ]_lk,dc
//      // circle
//      P[Pi_4q_comp1]      =  [ inv(Sin) * Gamma_proj * g5*adj(inv(Sout))*g5 ]_li,da * [ V_4q_{ijkl;abcd} ] * [ inv(Sin) * Gamma_proj * g5*adj(Sout)*g5 ]_jk,bc
//
//
//     
//
//      Generalised problem : mat_sc((si,sj),(ci,cj)) * mat_sc_sc((sk,sl),(ck,cl),(sm,sn),(cm,cn)) contract over some index eg delta_jk delta_il  or delta_jk delta_in 
//      Could do outer product then trace? 
//
////////////////////////////////////


/*
struct DiracStructure{
    std::vector<Gamma::Algebra> gammas;
    std::vector<int>            signs;
};
*/


/*
Real projectTree(DiracStructure vertex_structure, DiracStructure projector)
{
    ComplexD    figure8 = 0;
    ComplexD    circle  = 0;
    ComplexD    tr      = 0; 
    
    //vertex indices
    for(int nu=0;nu<vertex_structure.signs.size();nu++)
    {
        SpinColourMatrix vertex = Complex(1,0);
        vertex = vertex*Gamma(vertex_structure.gammas[nu]);
        //projector indices
        for(int mu=0;mu<projector.signs.size();mu++)
        {
            SpinColourMatrix leg = Complex(1,0);
            leg = leg*Gamma(projector.gammas[mu]);

            int ns(Grid::QCD::Ns);
            int nc(Grid::QCD::Nc);

            figure8 = 0;
            circle  = 0;

            figure8 = trace(leg*vertex)*trace(leg*vertex);
            circle  = trace(leg*vertex*leg*vertex);


            tr += vertex_structure.signs[nu]*projector.signs[mu]*2*(figure8-circle);
        }
    }
    return real(tr);
}
*/
Real projectFourQuark(Grid::QCD::SpinColourMatrix prop1, Grid::QCD::SpinColourMatrix prop2, std::vector<Grid::QCD::SpinColourSpinColourMatrix> vertices, DiracStructure vertex_structure,DiracStructure projector, bool colourMix)
{
    ComplexD    figure8 = 0;
    ComplexD    circle  = 0;
    ComplexD    tr      = 0; 
    
    //invert propagators
    auto propInv1 = invert(prop1);
    auto propInv2 = invert(prop2);

    // let's just do SS for now. Will then change to VV and AA. 
    // adjoint and g5 already saved in the props on disc.
    
    //vertex indices
    for(int nu=0;nu<vertex_structure.signs.size();nu++)
    {
        auto vertex = vertices[vertex_structure.indices[nu]];
        //projector indices
        for(int mu=0;mu<projector.signs.size();mu++)
        {

            auto sc_mat_1 = propInv1*projector.gammas[mu]*propInv2;
            auto sc_mat_2 = sc_mat_1;

            int ns(Grid::QCD::Ns);
            int nc(Grid::QCD::Nc);

            figure8 = 0;
            circle  = 0;

            // trace over spin and colour
            for(int cc=0;cc<nc;cc++)
            for(int cd=0;cd<nc;cd++)
            for(int sc=0;sc<ns;sc++)
            for(int sd=0;sd<ns;sd++)
            {
            for(int ca=0;ca<nc;ca++)
            for(int cb=0;cb<nc;cb++)
            for(int sa=0;sa<ns;sa++)
            for(int sb=0;sb<ns;sb++)
            {
                if(colourMix)
                {
                    figure8   += sc_mat_1()(sb,sa)(cd,ca)*vertex()(sa,sb)(ca,cb)(sc,sd)(cc,cd)*sc_mat_2()(sd,sc)(cb,cc);
                    circle    += sc_mat_1()(sd,sa)(cb,ca)*vertex()(sa,sb)(ca,cb)(sc,sd)(cc,cd)*sc_mat_2()(sb,sc)(cd,cc);
                }
                else
                {
                    figure8   += sc_mat_1()(sb,sa)(cb,ca)*vertex()(sa,sb)(ca,cb)(sc,sd)(cc,cd)*sc_mat_2()(sd,sc)(cd,cc);
                    circle    += sc_mat_1()(sd,sa)(cd,ca)*vertex()(sa,sb)(ca,cb)(sc,sd)(cc,cd)*sc_mat_2()(sb,sc)(cb,cc);
                }
            }
            }
            tr += vertex_structure.signs[nu]*projector.signs[mu]*2*(figure8-circle);
            //std::cout << nu << " " << mu << " " << 2*(figure8-circle) << std::endl; 
        }
    }
    return real(tr);
}

Eigen::MatrixXd projectFourQuark(Grid::QCD::SpinColourMatrix prop1, Grid::QCD::SpinColourMatrix prop2, std::vector<Grid::QCD::SpinColourSpinColourMatrix> vertices, std::vector<DiracStructure> vertex_structure,std::vector<DiracStructure> projector, std::vector<bool> colourMix)
{
    Eigen::MatrixXd trace(vertex_structure.size(),projector.size());
    for(int i=0;i<vertex_structure.size();i++)
    for(int j=0;j<projector.size();j++)
    {
        trace(i,j) = projectFourQuark(prop1,prop2,vertices,vertex_structure[i],projector[j],colourMix[j]);
    }
    return trace;
}

Distribution<Eigen::MatrixXd> projectFourQuark(Distribution<Grid::QCD::SpinColourMatrix> prop1, Distribution<Grid::QCD::SpinColourMatrix> prop2, Distribution<std::vector<Grid::QCD::SpinColourSpinColourMatrix>> vertices, std::vector<DiracStructure> vertex_structure,std::vector<DiracStructure> projector, std::vector<bool> colourMix)
{
    int nSamples = prop1.get_Nmeas();
    std::vector<Eigen::MatrixXd> trace(nSamples,Eigen::MatrixXd(vertex_structure.size(),projector.size()));
    trace.reserve(nSamples);

    for(int i=0;i<nSamples;i++)
    {
        trace[i] = projectFourQuark(prop1.get_value(i),prop2.get_value(i),vertices.get_value(i),vertex_structure,projector,colourMix);
    }
    return Distribution<Eigen::MatrixXd>(trace);

}

Distribution<Real> projectFourQuark(Distribution<Grid::QCD::SpinColourMatrix> prop1, Distribution<Grid::QCD::SpinColourMatrix> prop2, Distribution<std::vector<Grid::QCD::SpinColourSpinColourMatrix>> vertices, DiracStructure vertex_structure,DiracStructure projector, bool colourMix)
{
    int nSamples = prop1.get_Nmeas();

    std::vector<Real> trace;
    trace.reserve(nSamples);

    for(int i=0;i<nSamples;i++)
    {
        trace[i] = projectFourQuark(prop1.get_value(i),prop2.get_value(i),vertices.get_value(i),vertex_structure,projector,colourMix);
    }
    return Distribution<Real>(trace);
}


ComplexD projMix(SpinColourMatrix vertex, SpinColourMatrix leg)
{
    ComplexD    figure8;
    ComplexD    circle;

    SpinMatrix      tr_c = traceIndex<ColourIndex>(leg*vertex);
    ColourMatrix    tr_s = traceIndex<SpinIndex>(leg*vertex);

    figure8     = traceIndex<ColourIndex>(tr_s*tr_s);
    circle      = traceIndex<SpinIndex>(tr_c*tr_c);
            

    return figure8 - circle;
}

ComplexD projUnmix(SpinColourMatrix vertex, SpinColourMatrix leg)
{
    ComplexD    figure8;
    ComplexD    circle;

    figure8 = trace(leg*vertex)*trace(leg*vertex);
    circle  = trace(leg*vertex*leg*vertex);

    return figure8 - circle;
}

Real projectTree(DiracStructure vertex_structure, DiracStructure projector, bool colourMix)
{
    ComplexD    tr      = 0; 
    
    //vertex indices
    for(int nu=0;nu<vertex_structure.signs.size();nu++)
    {
        SpinColourMatrix vertex = Complex(1,0);
        vertex = vertex*vertex_structure.gammas[nu];
        //projector indices
        for(int mu=0;mu<projector.signs.size();mu++)
        {
            SpinColourMatrix leg = Complex(1,0);
            leg = leg*projector.gammas[mu];

            if(colourMix)   {   tr += vertex_structure.signs[nu]*projector.signs[mu]*2*projMix(vertex,leg);      }
            else            {   tr += vertex_structure.signs[nu]*projector.signs[mu]*2*projUnmix(vertex,leg);    }
        }
    }
    return round(real(tr));
}



#endif
