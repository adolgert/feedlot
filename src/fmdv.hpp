#ifndef _FMDV_HPP_
#define _FMDV_HPP_ 1

#include <vector>
#include "parameter.hpp"

template<typename ParamClass>
void FMDV_Mardones_Exponential(
    std::vector<TypedParameter<ParamClass>>& parameters) {
  using MyParm=TypedParameter<ParamClass>;
  parameters.emplace_back(MyParm{ParamClass::Latent, "latent", 1/3.59,
    "latent period"});
  parameters.emplace_back(MyParm{ParamClass::Gamma, "gamma", 1/4.39,
    "recovery rate"});
}

template<typename ParamClass>
void FMDV_Mardones_Nonexponential(
    std::vector<TypedParameter<ParamClass>>& parameters) {
  using MyParm=TypedParameter<ParamClass>;
  // Shape
  parameters.emplace_back(MyParm{ParamClass::LatentAlpha, "latentalpha", 1.782,
    "latent period alpha"});
  // Rate (not scale=1/rate)
  parameters.emplace_back(MyParm{ParamClass::LatentBeta, "latentbeta", 1/3.974,
    "latent period beta"});
  // Shape
  parameters.emplace_back(MyParm{ParamClass::GammaAlpha, "gammaalpha", 3.969,
    "recovery rate alpha"});
  // Rate (not scale=1/rate)
  parameters.emplace_back(MyParm{ParamClass::GammaBeta, "gammabeta", 1/1.107,
    "recovery rate beta"});
}


#endif
