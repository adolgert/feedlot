#ifndef _FEEDLOT_PLACES_HPP_
#define _FEEDLOT_PLACES_HPP_

struct IndividualToken
{
  IndividualToken()=default;

  inline friend
  std::ostream& operator<<(std::ostream& os, const IndividualToken& it){
    return os << "T";
  }
};


struct SIRPlace
{
  int64_t individual;
  int64_t location;
  int64_t disease;
  SIRPlace()=default;
  SIRPlace(int64_t i, int64_t l, int64_t d)
    : individual(i), location(l), disease(d) {}
  friend inline
  bool operator<(const SIRPlace& a, const SIRPlace& b) {
    return afidd::smv::LazyLess(a.individual, b.individual,
        a.location, b.location, a.disease, b.disease);
  }

  friend inline
  bool operator==(const SIRPlace& a, const SIRPlace& b) {
    return (a.individual==b.individual) && (a.location==b.location) &&
        (a.disease==b.disease);
  }

  friend inline
  std::ostream& operator<<(std::ostream& os, const SIRPlace& cp) {
    return os << '(' << cp.individual << ',' << cp.location << ',' <<
        cp.disease << ')';
  }
};


struct SIRTKey
{
  int64_t i;
  int64_t j;
  TransitionType kind;
  SIRTKey()=default;
  SIRTKey(int64_t i, int64_t j, TransitionType k) : i(i), j(j), kind(k) {}

  friend inline
  bool operator<(const SIRTKey& a, const SIRTKey& b) {
    return afidd::smv::LazyLess(a.i, b.i, a.j, b.j, a.kind, b.kind);
  }

  friend inline
  bool operator==(const SIRTKey& a, const SIRTKey& b) {
    return (a.i==b.i) && (a.j==b.j) && (a.kind==b.kind);
  }

  friend inline
  std::ostream& operator<<(std::ostream& os, const SIRTKey& cp) {
    return os << '(' << cp.i << ',' << cp.j << ',' << cp.kind << ')';
  }
};

#endif
