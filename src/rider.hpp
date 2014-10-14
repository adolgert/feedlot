#ifndef _SIR_EXP_HPP_
#define _SIR_EXP_HPP_ 1
#include <map>
#include <random>
#include "boost/random/mersenne_twister.hpp"
#include "mt19937.hpp"

using RandGen=afidd::rng::mt19937;
//using RandGen=boost::mt19937;

enum class SIRParam { Beta0, Beta1, Beta2, Gamma, RiderMove, RiderRecover,
    RiderInfect, RiderGetInfected };
struct Parameter {
  SIRParam kind;
  std::string name;
  double value;
  std::string description;
  Parameter(SIRParam k, std::string n, double v, std::string desc)
  : kind(k), name(n), value(v), description(desc) {}
  Parameter()=default;
};

struct TrajectoryEntry {
  int64_t s;
  int64_t e;
  int64_t i;
  int64_t r;
  double t;
  TrajectoryEntry(int64_t s, int64_t e, int64_t i, int64_t r, double t)
  : s(s), e(e), i(i), r(r), t(t) {}
  TrajectoryEntry()=default;
};

class TrajectoryObserver
{
public:
  virtual void Step(TrajectoryEntry sirt)=0;
  virtual const std::vector<TrajectoryEntry>& Trajectory() const =0;
};


int64_t SEIR_run(double time_limit, const std::vector<int64_t>& seir_cnt,
    const std::vector<Parameter>& parameters, TrajectoryObserver& observer,
    RandGen& rng, int block_cnt, int row_cnt);

#endif
