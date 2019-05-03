#ifndef GLOBAL_H
#define GLOBAL_H

#include "Grid/Grid.h"

// set gamma indices for projection - S,P,V,A,T
const std::vector<Gamma::Algebra> I           = {Gamma::Algebra::Identity};
const std::vector<Gamma::Algebra> g5          = {Gamma::Algebra::Gamma5};
const std::vector<Gamma::Algebra> gmu         = {Gamma::Algebra::GammaX,Gamma::Algebra::GammaY,Gamma::Algebra::GammaZ,Gamma::Algebra::GammaT};
const std::vector<Gamma::Algebra> gmug5       = {Gamma::Algebra::GammaXGamma5,Gamma::Algebra::GammaYGamma5,Gamma::Algebra::GammaZGamma5,Gamma::Algebra::GammaTGamma5};
const std::vector<Gamma::Algebra> sigma_mu_nu = {Gamma::Algebra::SigmaXY,Gamma::Algebra::SigmaXZ,Gamma::Algebra::SigmaXT,
                                                Gamma::Algebra::SigmaYZ,Gamma::Algebra::SigmaYT,Gamma::Algebra::SigmaZT};
#define ColourIndex     2
#define SpinIndex       1
#define Nc              3
#define Ns              4
#define Nd              4
#define Nop             5


#endif
