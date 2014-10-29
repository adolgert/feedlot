#include "rider_enums.hpp"


std::ostream& operator<<(std::ostream& os, TransitionType t) {
  switch (t) {
    case TransitionType::infect0 :
      os << "i0";
      break;
    case TransitionType::infect1 :
      os << "i1";
      break;
    case TransitionType::infect2 :
      os << "i2";
      break;
    case TransitionType::infectious :
      os << "is";
      break;
    case TransitionType::recover :
      os << "r";
      break;
    case TransitionType::movers :
      os << "ms";
      break;
    case TransitionType::moveri :
      os << "mi";
      break;
    case TransitionType::recoverr :
      os << "rr";
      break;
    case TransitionType::infectr :
      os << "ir";
      break;
    case TransitionType::infectbyr :
      os << "ib";
      break;
    case TransitionType::subclinical :
      os << "sc";
      break;
    default:
      os << "unknown transition";
      break;
  }
  return os;
}
