
#include <tuple>
#include <map>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <memory>
#include <set>
#include <functional>
#include "stochnet.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/log/core.hpp"
#include "boost/math/constants/constants.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/property_map/dynamic_property_map.hpp"
#include "boost/graph/graphml.hpp"
#include "boost/mpl/vector.hpp"
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "gspn_random.hpp"
#include "rider_enums.hpp"
#include "parameter.hpp"
#include "rider.hpp"

namespace smv=afidd::smv;
using namespace smv;

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
    default:
      os << "unknown transition";
      break;
  }
  return os;
}

struct IndividualToken
{
  int64_t id;
  IndividualToken()=default;
  IndividualToken(int64_t id) : id(id) {}

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
    return LazyLess(a.individual, b.individual, a.location, b.location,
        a.disease, b.disease);
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
    return LazyLess(a.i, b.i, a.j, b.j, a.kind, b.kind);
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


// This is as much of the marking as the transition will see.
using Local=LocalMarking<Uncolored<IndividualToken>>;
// Extra state to add to the system state. Will be passed to transitions.
struct WithParams {
  // Put our parameters here.
  std::map<SIRParam,double> params;
};


// The transition needs to know the local marking and any extra state.
using SIRTransition=ExplicitTransition<Local,RandGen,WithParams>;

using Dist=TransitionDistribution<RandGen>;
using ExpDist=ExponentialDistribution<RandGen>;



// For one infected to infect many susceptibles, the hazard
// is just multiplied by the number of susceptibles, so we can
// coalesce this transition. That means that, if it fires, it
// picks, randomly, which susceptible will become latent.
class InfectPen : public SIRTransition
{
  int64_t susceptible_cnt_;
public:
  InfectPen(int64_t susceptible_cnt) : susceptible_cnt_(susceptible_cnt) {}

  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    // If these are just size_t, then the rate calculation overflows.
    int64_t I=lm.template Length<0>(0);
    int64_t S=lm.template Length<0>(1);
    double rate=S*s.params.at(SIRParam::Beta0);
    if (S>0 && I>0 && rate>0.0) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"IP enable "<<S<< " I "<<I);
      return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"infection disable");
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    //SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm);
    // s0 i1 r2 i3 r4
    // If all Susceptibles are together, and we are taking one, which
    // exposed gets that susceptible? Look at its ID within the pen.
    auto ta=lm.template GetToken<0>(1, [](const IndividualToken& t)->int64_t {
      return t.id;
    });
    assert(ta.second==true);
    int64_t individual_idx=ta.first;
    const auto checker=[individual_idx](IndividualToken& t) {
      assert(t.id==individual_idx);
    };
    lm.template Move<0,0,decltype(checker)>(1, 2, 1, checker);
    lm.template Add<0>(3+individual_idx, IndividualToken{individual_idx});
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"IP Fire "<<individual_idx);
  }
};


class InfectFence : public SIRTransition
{
  int64_t susceptible_cnt_;
public:
  InfectFence(int64_t susceptible_cnt) : susceptible_cnt_(susceptible_cnt) {}

  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    // If these are just size_t, then the rate calculation overflows.
    int64_t I=lm.template Length<0>(0);
    int64_t S=lm.template Length<0>(1);
    double rate=S*s.params.at(SIRParam::Beta1);
    if (S>0 && I>0 && rate>0.0) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"IF enable "<<S<< " I "<<I);
      return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"infection disable");
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    //SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm);
    // s0 i1 r2 i3 r4
    auto ta=lm.template GetToken<0>(1, [](const IndividualToken& t)->int64_t {
      return t.id;
    });
    assert(ta.second==true);
    int64_t individual_idx=ta.first;
    const auto checker=[individual_idx](IndividualToken& t) {
      assert(t.id==individual_idx);
    };
    lm.template Move<0,0,decltype(checker)>(1, 2, 1, checker);
    lm.template Add<0>(3+individual_idx, IndividualToken{individual_idx});
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"IF Fire "<<individual_idx);
  }
};


class InfectOther : public SIRTransition
{
  int64_t susceptible_cnt_;
public:
  InfectOther(int64_t susceptible_cnt) : susceptible_cnt_(susceptible_cnt) {}

  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    // If these are just size_t, then the rate calculation overflows.
    int64_t I=lm.template Length<0>(0);
    int64_t S=lm.template Length<0>(1);
    double rate=S*s.params.at(SIRParam::Beta2);
    if (S>0 && I>0 && rate>0.0) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"IO enable "<<S<< " I "<<I);
      return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"infection disable");
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    //SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm);
    // s0 i1 r2 i3 r4
    auto ta=lm.template GetToken<0>(1, [](const IndividualToken& t)->int64_t {
      return t.id;
    });
    assert(ta.second==true);
    int64_t individual_idx=ta.first;
    const auto checker=[individual_idx](IndividualToken& t) {
      assert(t.id==individual_idx);
    };
    lm.template Move<0,0,decltype(checker)>(1, 2, 1, checker);
    lm.template Add<0>(3+individual_idx, IndividualToken{individual_idx});
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"IO Fire "<<individual_idx);
  }
};


// Now make specific transitions.
class InfectiousExponential : public SIRTransition
{
  int64_t individual_;
public:
  InfectiousExponential(int64_t individual) : individual_(individual) {}
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      double rate=I*s.params.at(SIRParam::Gamma);
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<ExpDist>(
        new ExpDist(rate, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    SMVLOG(BOOST_LOG_TRIVIAL(debug) << "Fire infectious " <<individual_
        << " marking " << lm);
    lm.template Move<0, 0>(0, 2, 1); // Change the individual
    lm.template Move<0, 0>(1, 3, 1); // Change the summary count
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Infectious Fire ");
  }
};

template<typename BaseTransition>
class Infectious : public BaseTransition
{
public:
  virtual std::pair<bool, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>>
  Enabled(const typename BaseTransition::UserState& s,
    const typename BaseTransition::LocalMarking& lm,
    double te, double t0, typename BaseTransition::RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Infectious enable I "<<I);
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<
        afidd::smv::WeibullDistribution<typename BaseTransition::RandGen>>(
        new afidd::smv::WeibullDistribution<typename BaseTransition::RandGen>(
          1/s.params.at(SIRParam::LatentBeta),
          s.params.at(SIRParam::LatentAlpha), te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>(nullptr)};
    }
  }

  virtual void Fire(typename BaseTransition::UserState& s,
    typename BaseTransition::LocalMarking& lm, double t0,
      typename BaseTransition::RandGen& rng) override {
    lm.template Move<0, 0>(0, 2, 1); // Change the individual
    lm.template Move<0, 0>(1, 3, 1); // Change the summary count
  }
};

// Now make specific transitions.
class RecoverExponential : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      double rate=I*s.params.at(SIRParam::Gamma);
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<ExpDist>(
        new ExpDist(rate, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    SMVLOG(BOOST_LOG_TRIVIAL(debug) << "Fire recover " << lm);
    lm.template Remove<0>(0, 1, rng);
    lm.template Move<0, 0>(1, 2, 1);
  }
};


// Now make specific transitions.
template<typename BaseTransition>
class Recover : public BaseTransition
{
  virtual std::pair<bool, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>>
  Enabled(const typename BaseTransition::UserState& s,
    const typename BaseTransition::LocalMarking& lm,
    double te, double t0, typename BaseTransition::RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"recover enable "<< I<<" te "<<te
          <<" t0 "<<t0);
      return {true, std::unique_ptr<afidd::smv::GammaDistribution<
        typename BaseTransition::RandGen>>(
        new afidd::smv::GammaDistribution<
        typename BaseTransition::RandGen>(
          s.params.at(SIRParam::GammaAlpha),
          s.params.at(SIRParam::GammaBeta), te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>(nullptr)};    }
  }

  virtual void Fire(typename BaseTransition::UserState& s,
    typename BaseTransition::LocalMarking& lm, double t0,
      typename BaseTransition::RandGen& rng) override {
    SMVLOG(BOOST_LOG_TRIVIAL(debug) << "Fire recover " << lm);
    lm.template Remove<0>(0, 1, rng);
    lm.template Move<0, 0>(1, 2, 1);
  }
};


// Now make specific transitions.
class StartRider : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"Start rider enable "<<I);
    if (I>0) {
      // ta is the next occurrence of 6am.
      double ta=std::ceil(te)-te;
      double tb=ta + 15.0/(60.0*24.0);
      SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"Start rider enable "<<I<<
          " ta "<<ta<<" tb "<<tb<<" te "<<te<<" t0 "<<t0);
      return {true, std::unique_ptr<Dist>(
        new afidd::smv::UniformDistribution<RandGen>(ta, tb, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Start rider fire " << lm);
    lm.template Move<0, 0>(0, 1, 1);
  }
};

// Now make specific transitions.
class MoveRider : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      double rate=I*s.params.at(SIRParam::RiderMove);
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Move rider enable");
      return {true, std::unique_ptr<ExpDist>(
        new ExpDist(rate, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Move rider fire " << lm);
    lm.template Move<0, 0>(0, 1, 1);
  }
};


// Now make specific transitions.
class RecoverRider : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      double rate=I*s.params.at(SIRParam::RiderRecover);
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"recover rider rate "<< rate);
      return {true, std::unique_ptr<ExpDist>(
        new ExpDist(rate, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    //SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Fire recover " << lm);
    lm.template Move<0, 0>(0, 1, 1);
  }
};


// Now make specific transitions.
class InfectRider : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t S=lm.template Length<0>(0);
    int64_t I=lm.template Length<0>(1);
    if (S>0 && I>0) {
      double rate=s.params.at(SIRParam::RiderGetInfected);
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"infect rider "<< rate
        <<" at "<<t0);
      return {true, std::unique_ptr<ExpDist>(
        new ExpDist(rate, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    //SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Fire recover " << lm);
    lm.template Move<0, 0>(0, 2, 1);
  }
};


// Now make specific transitions.
class RiderInfects : public SIRTransition
{
  int64_t susceptible_cnt_;
public:
  RiderInfects(int64_t susceptible_cnt) : susceptible_cnt_(susceptible_cnt) {}

  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    // If these are just size_t, then the rate calculation overflows.
    int64_t I=lm.template Length<0>(0);
    int64_t S=lm.template Length<0>(1);
    double rate=S*s.params.at(SIRParam::RiderInfect);
    if (S>0 && I>0 && rate>0.0) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Rider infects enable"<<S<< " I "<<I);
      return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"infection disable");
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    //SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm);
    // s0 i1 r2 i3 r4
    auto ta=lm.template GetToken<0>(1, [](const IndividualToken& t)->int64_t {
      return t.id;
    });
    assert(ta.second==true);
    int64_t individual_idx=ta.first;
    const auto checker=[individual_idx](IndividualToken& t) {
      assert(t.id==individual_idx);
    };
    lm.template Move<0,0,decltype(checker)>(1, 2, 1, checker);
    lm.template Add<0>(3+individual_idx, IndividualToken{individual_idx});
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Rider infects Fire "<<individual_idx);
  }
};

// This kind of adjacency list uses integers for the vertex id, which
// is special to those using vecS and vecS.
using PenContactGraph=boost::adjacency_list<boost::vecS,
    boost::vecS,boost::undirectedS>;
/*! Feedlots look like suburbs. There are long blocks with streets
 *  between. Each block has row_cnt pens, sitting back-to-back.
 *  Returns block_cnt*row_cnt*2 pens.
 */
PenContactGraph BlockStructure(int block_cnt, int row_cnt) {
  PenContactGraph g(block_cnt*row_cnt*2);
  for (int bidx=0; bidx<2*block_cnt; ++bidx) {
    int base=bidx*row_cnt;
    for (int right_idx=0; right_idx<row_cnt-1; ++right_idx) {
      add_edge(base+right_idx, base+right_idx+1, g);
    }
  }
  return g;
}

bool AdjacentPens(int i, int j, const PenContactGraph& g) {
  using AdjIter=boost::graph_traits<PenContactGraph>::adjacency_iterator;
  AdjIter start, end;
  assert(i<num_vertices(g));
  std::tie(start, end)=adjacent_vertices(i, g);
  for ( ; start!=end; ++start) {
    if (*start==j) {
      return true;
    }
  }
  return false;
}

int64_t pen_of(int64_t individual, int64_t per_pen) {
  return individual/per_pen;
}

// The GSPN itself.
using SIRGSPN=
    ExplicitTransitions<SIRPlace, SIRTKey, Local, RandGen, WithParams>;

/*! SIR infection on an all-to-all graph of uncolored tokens.
 */
void BuildSystem(SIRGSPN& bg, int64_t individual_cnt,
    const PenContactGraph& g)
{
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

  int64_t pen_cnt=num_vertices(g);
  int64_t per_pen=individual_cnt/pen_cnt;
  assert((individual_cnt % per_pen)==0);
  BOOST_LOG_TRIVIAL(info) << "individuals "<<individual_cnt<<" pens "<<pen_cnt
      <<" per pen "<<per_pen;
  enum : int64_t { s, e, i, r };

  // 2N
  const int64_t location=0;
  for (int64_t ap_idx=0; ap_idx<individual_cnt; ++ap_idx) {
    for (int64_t place : std::vector<int64_t>{e, i}) {
      bg.AddPlace({ap_idx, location, place}, 0);
    }
  }
  // 4P
  // Give each pen a summary count of disease states within.
  const int64_t pen_summary=1;
  for (int64_t sp_idx=0; sp_idx<pen_cnt; ++sp_idx) {
    for (int64_t place : std::vector<int64_t>{s, e, i, r}) {
      bg.AddPlace({sp_idx, pen_summary, place}, 0);
    }
  }

  const int64_t rider=2;
  for (int64_t pp_idx=0; pp_idx<pen_cnt+1; ++pp_idx) {
    bg.AddPlace({pp_idx, rider, s}, 0);
    bg.AddPlace({pp_idx, rider, i}, 0);
  }

  // 2N
  for (int64_t ind_idx=0; ind_idx<individual_cnt; ++ind_idx) {
    int64_t ind_pen=pen_of(ind_idx, per_pen);
    bg.AddTransition({ind_idx, ind_idx, TransitionType::infectious},
      {Edge{{ind_idx, location, e}, -1},
       Edge{{ind_pen, pen_summary, e},-1}, Edge{{ind_idx, location, i}, 1},
       Edge{{ind_pen, pen_summary, i}, 1}},
      std::unique_ptr<SIRTransition>(new Infectious<SIRTransition>())
      );
    bg.AddTransition({ind_idx, ind_idx, TransitionType::recover},
      {Edge{{ind_idx, location, i}, -1},
       Edge{{ind_pen, pen_summary, i}, -1},
       Edge{{ind_pen, pen_summary, r}, 1}},
      std::unique_ptr<SIRTransition>(new Recover<SIRTransition>())
      );
  }

  // 2*CE*(N/P) + N + P (for the rider)
  std::vector<Edge> infect_vec(3+per_pen);
  for (int64_t d_idx=0; d_idx<pen_cnt; ++d_idx) {
    for (int64_t s_idx=0; s_idx<pen_cnt; ++s_idx) {
      int64_t s_base=s_idx*per_pen;
      infect_vec[1]=Edge{{s_idx, pen_summary, s}, -1};
      infect_vec[2]=Edge{{s_idx, pen_summary, e}, 1};
      for (int64_t targ_idx=0; targ_idx<per_pen; ++targ_idx) {
        infect_vec[3+targ_idx]=Edge{{s_base+targ_idx, location, e},1};
      }
      if (s_idx==d_idx) {
        int64_t src_base=d_idx*per_pen;
        for (int64_t src_idx=0; src_idx<per_pen; ++src_idx) {
          int64_t src=src_base+src_idx;
          infect_vec[0]=Edge{{src, location, i}, -1};
          bg.AddTransition({src, s_idx, TransitionType::infect0}, infect_vec,
            std::unique_ptr<SIRTransition>(new InfectPen(per_pen))
            );
        }
        // Also let the rider infect this pen.
        infect_vec[0]=Edge{{d_idx, rider, i}, -1};
        bg.AddTransition({d_idx, d_idx, TransitionType::infectbyr},
            infect_vec,
            std::unique_ptr<SIRTransition>(new RiderInfects(per_pen)));
      } else if (AdjacentPens(d_idx, s_idx, g)) {
        int64_t src_base=d_idx*per_pen;
        for (int64_t src_idx=0; src_idx<per_pen; ++src_idx) {
          int64_t src=src_base+src_idx;
          infect_vec[0]=Edge{{src, location, i}, -1};
          bg.AddTransition({src, s_idx, TransitionType::infect1}, infect_vec,
            std::unique_ptr<SIRTransition>(new InfectFence(per_pen))
            );
        }
      } else {
        ;
      }
    }
  }

  auto rider_out_s=SIRPlace{pen_cnt, rider, s};
  auto rider_out_i=SIRPlace{pen_cnt, rider, i};
  auto first_pen_s=SIRPlace{0, rider, s};
  auto first_pen_i=SIRPlace{0, rider, i};
  bg.AddTransition({pen_cnt, 0, TransitionType::movers},
    {Edge{rider_out_s, -1}, Edge{first_pen_s, 1}},
    std::unique_ptr<SIRTransition>(new StartRider()));
  bg.AddTransition({pen_cnt, 0, TransitionType::moveri},
    {Edge{rider_out_i, -1}, Edge{first_pen_i, 1}},
    std::unique_ptr<SIRTransition>(new StartRider()));
  //P (recover) + N (infect) + (CG+1)*2 (move)
  for (int64_t rp_idx=0; rp_idx<pen_cnt; ++rp_idx) {
    // Move
    auto pen_s=SIRPlace{rp_idx, rider, s};
    auto pen_i=SIRPlace{rp_idx, rider, i};
    int64_t next_pen=static_cast<int64_t>(rp_idx+1);
    BOOST_LOG_TRIVIAL(debug)<<"rider connects "<<rp_idx<<" to "<<next_pen;
    auto pen_s_n=SIRPlace{next_pen, rider, s};
    auto pen_i_n=SIRPlace{next_pen, rider, i};
    bg.AddTransition({rp_idx, next_pen, TransitionType::movers},
      {Edge{pen_s, -1}, Edge{pen_s_n, 1}},
      std::unique_ptr<SIRTransition>(new MoveRider()));
    bg.AddTransition({rp_idx, next_pen, TransitionType::moveri},
      {Edge{pen_i, -1}, Edge{pen_i_n, 1}},
      std::unique_ptr<SIRTransition>(new MoveRider()));
    bg.AddTransition({rp_idx, rp_idx, TransitionType::recoverr},
      {Edge{pen_i, -1}, Edge{pen_s, 1}},
      std::unique_ptr<SIRTransition>(new RecoverRider()));
    for (int64_t c_idx=rp_idx*per_pen; c_idx<(rp_idx+1)*per_pen; ++c_idx) {
      bg.AddTransition({rp_idx, c_idx, TransitionType::infectr},
        {Edge{pen_s, -1}, Edge{{c_idx, location, i},-1},
         Edge{pen_i, 1}},
         std::unique_ptr<SIRTransition>(new InfectRider()));
    }
  }
}


template<typename GSPN, typename Marking>
bool CheckMarking(const GSPN& gspn, const Marking& marking,
    int64_t individual_cnt, int64_t pen_cnt, int64_t per_pen) {
  BOOST_LOG_TRIVIAL(debug)<<"Checking for exactly one state "<<individual_cnt
      <<" "<<pen_cnt<<" "<<per_pen;
  // Check that each individual is in exactly one state.
  enum : int64_t { s, e, i, r };
  std::vector<int64_t> seir{0,0,0,0};
  const int64_t single=0;
  const int64_t pen_summary=1;
  const int64_t rider=2;

  int64_t ind_mark_total=0;
  for (int64_t cp_idx=0; cp_idx<pen_cnt; ++cp_idx) {
    for (int64_t ps_idx=0; ps_idx<4; ++ps_idx) {
      auto pen_place=gspn.PlaceVertex({cp_idx, pen_summary, ps_idx});
      int64_t pen_kind=Length<0>(marking, pen_place);
      ind_mark_total+=pen_kind;
      seir[ps_idx]+=pen_kind;
      if (ps_idx==1 || ps_idx==2) {
        int64_t kind_in_pen=0;
        for (int64_t c_idx=cp_idx*per_pen; c_idx<(cp_idx+1)*per_pen; ++c_idx) {
          auto c_place=gspn.PlaceVertex({c_idx, single, ps_idx});
          kind_in_pen+=Length<0>(marking, c_place);
        }
        if (kind_in_pen!=pen_kind) {
          BOOST_LOG_TRIVIAL(error)<<"CheckMarking expected "<<pen_kind<<
            " but found "<<kind_in_pen<<" tokens in pen "<<cp_idx<<" of type "
            <<ps_idx;
          assert(kind_in_pen==pen_kind);
        }
      }
    }
  }
  int64_t rider_cnt=0;
  for (int64_t ridx=0; ridx<pen_cnt+1; ++ridx) {
    auto rsplace=gspn.PlaceVertex({ridx, rider, s});
    rider_cnt+=Length<0>(marking, rsplace);
    auto riplace=gspn.PlaceVertex({ridx, rider, i});
    rider_cnt+=Length<0>(marking, riplace);
  }
  BOOST_LOG_TRIVIAL(debug)<<"Total of "<<ind_mark_total<<" tokens "
    <<"and "<<rider_cnt<<" riders";
  BOOST_LOG_TRIVIAL(debug)<<"Marking has "<<seir[0]<<", "<<seir[1]
      <<", "<<seir[2]<<", "<<seir[3];
  assert(ind_mark_total==individual_cnt);
  return true;
}

/*!
 * Given a vector of checkpoint times, for the state of the system
 * at each of those times, count the number of
 * susceptibles and infecteds at that time. Form into two vectors,
 * one for susceptibles, one for infecteds. Then repeat the whole
 * simulation 10^4 times. Return these vectors.
 */
template<typename GSPN, typename SIRState>
struct SEIROutput
{
  using StateArray=std::array<int64_t,4>;
  std::shared_ptr<PenTrajectoryObserver> observer_;
  const GSPN& gspn_;
  NonHomogeneousPoissonProcesses<int64_t,RandGen>& propagator_;
  StateArray seir_;
  int64_t step_cnt{0};
  int64_t individual_cnt_;
  int64_t per_pen_;
  int64_t pen_cnt_;

  SEIROutput(const GSPN& gspn,
      NonHomogeneousPoissonProcesses<int64_t,RandGen>& propagator,
      std::shared_ptr<PenTrajectoryObserver> observer,
      const std::vector<int64_t>& initial, int64_t per_pen)
  : gspn_(gspn), propagator_(propagator), observer_(observer), per_pen_(per_pen)
  {
    std::get<0>(seir_)=initial[0];
    std::get<1>(seir_)=initial[1];
    std::get<2>(seir_)=initial[2];
    std::get<3>(seir_)=initial[3];
    individual_cnt_=std::accumulate(initial.begin(), initial.end(), int64_t{0});
    pen_cnt_=individual_cnt_/per_pen_;
  };

  enum : int64_t { s, e, i, r };
  enum : int64_t { none, infect0, infect1, infect2, infectious, recover };
  bool operator()(const SIRState& state) {
    auto transition=gspn_.VertexTransition(state.last_transition);
    int64_t individual=transition.i;
    int64_t compartment=0;
    bool affected=true;
    switch (transition.kind) {
      case TransitionType::infect0: // infect
      case TransitionType::infect1:
      case TransitionType::infect2:
        compartment=0;
        get<0>(seir_)-=1;
        get<1>(seir_)+=1;
        break;
      case TransitionType::infectbyr:
        compartment=0;
        individual=transition.j;
        get<0>(seir_)-=1;
        get<1>(seir_)+=1;
        break;
      case TransitionType::infectious:
        compartment=1;
        get<1>(seir_)-=1;
        get<2>(seir_)+=1;
        break;
      case TransitionType::recover: // recover
        compartment=2;
        get<2>(seir_)-=1;
        get<3>(seir_)+=1;
        break;
      default:
        affected=false;
        if (transition.kind!=TransitionType::movers) {
          SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"observer transition kind: "
              <<transition.kind);
        }
        break;
    }

    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"observer: transition "<<transition);
    if (affected) {
      int64_t pen=pen_of(individual, per_pen_);
      observer_->Step({individual, pen, compartment, state.CurrentTime()});
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"step cnt "<<step_cnt<<" e "
          <<std::get<1>(seir_)<<" i "<<std::get<2>(seir_)<<" r "
          <<std::get<3>(seir_));
    }
    //CheckMarking(gspn_, state.marking, individual_cnt_, pen_cnt_, per_pen_);

    // Check that every infectious individual has a recovery transition enabled.
    // We need allowance b/c a newly-infectious individual won't have its
    // recovery yet.
    // int64_t allowance=0;
    // if (transition.kind==TransitionType::infectious) {
    //   allowance=1;
    // }
    // for (int64_t iidx=0; iidx<individual_cnt_; ++iidx) {
    //   auto sp=SIRPlace{iidx,0,i};
    //   int64_t pv=gspn_.PlaceVertex(sp);
    //   if (Length<0>(state.marking, pv)==1) {
    //     int64_t tv=gspn_.TransitionVertex({iidx, iidx, TransitionType::recover});
    //     auto enabled=propagator_.Enabled(tv);
    //     double firing_time=-1;
    //     if (!std::get<0>(enabled)) {
    //       allowance-=1;
    //     } else {
    //       firing_time=propagator_.FiringTime(tv);
    //     }
    //     BOOST_LOG_TRIVIAL(debug)<<"individual still infectious "<<sp
    //         <<" with recovery enabling time "<<std::get<1>(enabled)
    //         <<" firing time "<<firing_time
    //         <<" at current time "<<state.CurrentTime();
    //   }
    // }
    // if (allowance<0) {
    //   BOOST_LOG_TRIVIAL(error)<<"infectious without recovery "<<-allowance;
    //   assert(allowance>=0);
    // }

    ++step_cnt;
    return (std::get<1>(seir_)+std::get<2>(seir_)>0);
  }

  void final(const SIRState& state) {
    BOOST_LOG_TRIVIAL(info) << "Took "<< step_cnt << " transitions.";
  }
};



int64_t SEIR_run(double end_time, const std::vector<int64_t>& seir_cnt,
    const std::vector<TypedParameter<SIRParam>>& parameters,
    std::shared_ptr<PenTrajectoryObserver> observer,
    RandGen& rng, int block_cnt, int row_cnt)
{
  int64_t individual_cnt=std::accumulate(seir_cnt.begin(), seir_cnt.end(),
    int64_t{0});
  auto pen_contact=BlockStructure(block_cnt, row_cnt);
  // The goal is to put all latent and infecteds in the same pen.
  int64_t pen_cnt=num_vertices(pen_contact);
  int64_t animals_per_pen=individual_cnt/pen_cnt;

  int64_t N=individual_cnt;
  int64_t P=num_vertices(pen_contact);
  int64_t CE=num_edges(pen_contact);
  int64_t animal_places=2*N + P*4;
  int64_t rider_places=(P+1)*2;
  int64_t individual_transitions=2*N; // infectious and recover.
  int64_t animal_infect=(P+CE*2)*animals_per_pen;
  int64_t rider_infect=P+N; // infect and infected by.
  int64_t rider_move=(P+1)*2;
  int64_t rider_recover=P;
  int64_t guess_cnt=animal_places+rider_places+individual_transitions
      +animal_infect+rider_infect+rider_move+rider_recover;
  SIRGSPN gspn(guess_cnt);
  BuildSystem(gspn, individual_cnt, pen_contact);
  BOOST_LOG_TRIVIAL(debug)<<"GSPN vertex count "<<gspn.VerticesUsed()
      << " guess " << guess_cnt;

  // Marking of the net.
  static_assert(std::is_same<int64_t,SIRGSPN::PlaceKey>::value,
    "The GSPN's internal place type is int64_t.");
  using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<IndividualToken>>;
  using SIRState=GSPNState<Mark,SIRGSPN::TransitionKey,WithParams>;

  SIRState state;
  std::set<SIRParam> scale_param{ SIRParam::Beta0, SIRParam::Beta1,
      SIRParam::Beta2, SIRParam::RiderInfect };
  for (auto& cp : parameters) {
    double value=cp.value;
    if (scale_param.find(cp.kind)!=scale_param.end()) {
      value=cp.value/animals_per_pen;
    }
    state.user.params[cp.kind]=value;
  }

  BOOST_LOG_TRIVIAL(debug)<<"Creating susceptibles. "<<individual_cnt;
  const int64_t location=0;
  const int64_t sumloc=1;
  for (int64_t sus_idx=0; sus_idx<individual_cnt; ++sus_idx) {
    auto pen_idx=pen_of(sus_idx, animals_per_pen);
    auto summary_id=gspn.PlaceVertex({pen_idx, sumloc, 0});
    Add<0>(state.marking, summary_id, IndividualToken{sus_idx%animals_per_pen});
  }

  int64_t infected_pen=smv::uniform_index(rng, pen_cnt);
  BOOST_LOG_TRIVIAL(debug)<<"Moving susceptibles to other states. First "
      <<"infected pen is "<<infected_pen;
  int64_t first_in_pen=animals_per_pen*infected_pen;
  for (int reinit=1; reinit<4; ++reinit) {
    for (int64_t mv_idx=0; mv_idx<seir_cnt[reinit]; ++mv_idx) {
      if (reinit>0 && reinit<3) {
        auto to_place=SIRPlace{first_in_pen, location, reinit};
        auto to_id=gspn.PlaceVertex(to_place);
        Add<0>(state.marking, to_id, IndividualToken{});
      }
      auto move_pen=pen_of(first_in_pen, animals_per_pen);
      auto sus_pen=gspn.PlaceVertex({move_pen, sumloc, 0});
      auto to_pen_place=SIRPlace{move_pen, sumloc, reinit};
      auto to_pen=gspn.PlaceVertex(to_pen_place);
      Move<0,0>(state.marking, sus_pen, to_pen, 1);
      ++first_in_pen;
    }
  }

  // The rider starts outside the pens.
  auto rider_id=gspn.PlaceVertex({0, 2, 0});
  Add<0>(state.marking, rider_id, IndividualToken{});
  // BOOST_LOG_TRIVIAL(debug)<<"Checking marking just after building.";
  // CheckMarking(gspn, state.marking, individual_cnt, pen_cnt, animals_per_pen);

  //using Propagator=PropagateCompetingProcesses<int64_t,RandGen>;
  using Propagator=NonHomogeneousPoissonProcesses<int64_t,RandGen>;
  Propagator competing;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RandGen>;
  Dynamics dynamics(gspn, {&competing});

  SEIROutput<SIRGSPN,SIRState> output_function(gspn, competing,
    observer, seir_cnt, animals_per_pen);

  dynamics.Initialize(&state, &rng);

  BOOST_LOG_TRIVIAL(info)<<"Starting main loop";
  bool running=true;
  auto nothing=[](SIRState&)->void {};
  double last_time=state.CurrentTime();
  while (running && state.CurrentTime()<end_time) {
    running=dynamics(state);
    if (running) {
      double new_time=state.CurrentTime();
      if (new_time-last_time<-1e-12) {
        BOOST_LOG_TRIVIAL(warning) << "last time "<<last_time <<" "
          << " new_time "<<new_time;
      }
      last_time=new_time;
      running=output_function(state);
    } else {
      BOOST_LOG_TRIVIAL(info)<<"No transitions left to fire "
          <<state.CurrentTime();
    }
    // auto v=competing.content_size();
    // SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Competing Processes: size "
    //     <<v.first<<" infinities "<<v.second);
  }
  BOOST_LOG_TRIVIAL(info)<<"Reached end time "<<state.CurrentTime();
  output_function.final(state);
  return 0;
}

