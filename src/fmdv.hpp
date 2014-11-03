#ifndef _FMDV_HPP_
#define _FMDV_HPP_ 1

#include <vector>
#include "parameter.hpp"

template<typename ParamClass>
void FMDV_Mardones_Exponential(
    std::map<ParamClass,TypedParameter<ParamClass>>& parameters) {
  using MyParm=TypedParameter<ParamClass>;
  parameters[ParamClass::Latent]=MyParm{ParamClass::Latent, "latent", 1/3.59,
    "latent period"};
  parameters[ParamClass::Gamma]=MyParm{ParamClass::Gamma, "gamma", 1/4.39,
    "recovery rate"};
}

template<typename ParamClass>
void FMDV_Mardones_Nonexponential(
    std::map<ParamClass,TypedParameter<ParamClass>>& parameters) {
  using MyParm=TypedParameter<ParamClass>;
  // Shape
  parameters[ParamClass::LatentAlpha]=MyParm{ParamClass::LatentAlpha, "latentalpha", 1.782,
    "latent period alpha"};
  // Rate (not scale=1/rate)
  parameters[ParamClass::LatentBeta]=MyParm{ParamClass::LatentBeta, "latentbeta", 1/3.974,
    "latent period beta"};
  // Shape
  parameters[ParamClass::GammaAlpha]=MyParm{ParamClass::GammaAlpha, "gammaalpha", 3.969,
    "recovery rate alpha"};
  // Rate (not scale=1/rate)
  parameters[ParamClass::GammaBeta]=MyParm{ParamClass::GammaBeta, "gammabeta", 1/1.107,
    "recovery rate beta"};
  parameters[ParamClass::SubClinicalAlpha]=MyParm{ParamClass::SubClinicalAlpha, "scalpha",
    1.22, "Subclinical period shape"};
  parameters[ParamClass::SubClinicalBeta]=MyParm{ParamClass::SubClinicalBeta, "scbeta",
    1.672, "Subclinical period scale"};
}


#endif
