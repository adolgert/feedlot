#ifndef _RIDER_ENUMS_HPP_
#define _RIDER_ENUMS_HPP_ 1
#include <iostream>

enum class SIRParam { Beta0, Beta1, Beta2, Latent, Gamma, RiderMove, RiderRecover,
    RiderInfect, RiderGetInfected, LatentAlpha, LatentBeta,
    GammaAlpha, GammaBeta, SubClinicalAlpha, SubClinicalBeta };

enum class TransitionType : int64_t { none, infect0, infect1, infect2,
  infectious, recover, movers, moveri, recoverr, infectr, infectbyr,
  reset, subclinical };

std::ostream& operator<<(std::ostream&, TransitionType t);

#endif
