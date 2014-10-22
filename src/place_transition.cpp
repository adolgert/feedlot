#include "smv_algorithm.hpp"
#include "place_transition.hpp"


namespace smv=afidd::smv;
using namespace smv;

SIRPlace::SIRPlace(int64_t i, int64_t l, int64_t d)
  : individual(i), location(l), disease(d) {}

bool operator<(const SIRPlace& a, const SIRPlace& b) {
  return smv::LazyLess(a.individual, b.individual, a.location, b.location,
      a.disease, b.disease);
}

bool operator==(const SIRPlace& a, const SIRPlace& b) {
  return (a.individual==b.individual) && (a.location==b.location) &&
      (a.disease==b.disease);
}

std::ostream& operator<<(std::ostream& os, const SIRPlace& cp) {
  return os << '(' << cp.individual << ',' << cp.location << ',' <<
      cp.disease << ')';
}


SIRTKey::SIRTKey(int64_t i, int64_t j, TransitionType k) : i(i), j(j), kind(k) {}

bool operator<(const SIRTKey& a, const SIRTKey& b) {
  return smv::LazyLess(a.i, b.i, a.j, b.j, a.kind, b.kind);
}

bool operator==(const SIRTKey& a, const SIRTKey& b) {
  return (a.i==b.i) && (a.j==b.j) && (a.kind==b.kind);
}

std::ostream& operator<<(std::ostream& os, const SIRTKey& cp) {
  return os << '(' << cp.i << ',' << cp.j << ',' << cp.kind << ')';
}
