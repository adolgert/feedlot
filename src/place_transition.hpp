#ifndef _PLACE_TRANSITION_HPP_
#define _PLACE_TRANSITION_HPP_ 1

#include <cstdint>
#include <cstddef>
#include "rider_enums.hpp"

struct SIRPlace
{
  int64_t individual;
  int64_t location;
  int64_t disease;
  SIRPlace()=default;
  SIRPlace(int64_t i, int64_t l, int64_t d);
  friend bool operator<(const SIRPlace& a, const SIRPlace& b);

  friend bool operator==(const SIRPlace& a, const SIRPlace& b);

  friend std::ostream& operator<<(std::ostream& os, const SIRPlace& cp);
};


struct SIRTKey
{
  int64_t i;
  int64_t j;
  TransitionType kind;
  SIRTKey()=default;
  SIRTKey(int64_t i, int64_t j, TransitionType k);
  friend bool operator<(const SIRTKey& a, const SIRTKey& b);
  friend bool operator==(const SIRTKey& a, const SIRTKey& b);
  friend std::ostream& operator<<(std::ostream& os, const SIRTKey& cp);
};


#endif
