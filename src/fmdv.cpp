#include "fmdv.hpp"
#include "rider_enums.hpp"


template<typename ParamClass>
void FMDV_Mardones_Exponential(std::vector<ParamClass>& parameters) {
  parameters.emplace_back(MyParm{ParamClass::Latent, "latent", 1/3.59,
    "latent period"})
  parameters.emplace_back(MyParm{ParamClass::Gamma, "gamma", 1/4.39,
    "recovery rate"});
}


void FMDV_Mardones_Nonexponential(std::vector<ParamClass>& parameters) {
  parameters.emplace_back(MyParm{ParamClass::LatentAlpha, "latentalpha", 1/1.782,
    "latent period alpha"})
  parameters.emplace_back(MyParm{ParamClass::LatentBeta, "latentbeta", 1/3.974,
    "latent period beta"})
  parameters.emplace_back(MyParm{ParamClass::GammaAlpha, "gammaalpha", 1/3.969,
    "recovery rate alpha"});
  parameters.emplace_back(MyParm{ParamClass::GammaBeta, "gammabeta", 1/1.107,
    "recovery rate beta"});
}
