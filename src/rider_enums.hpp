#ifndef _RIDER_ENUMS_HPP_
#define _RIDER_ENUMS_HPP_ 1

enum class SIRParam { Beta0, Beta1, Beta2, Gamma, RiderMove, RiderRecover,
    RiderInfect, RiderGetInfected };


enum class TransitionType : int64_t { none, infect0, infect1, infect2,
  infectious, recover, movers, moveri, recoverr, infectr, infectbyr };

std::ostream& operator<<(std::ostream&, TransitionType t);

#endif
