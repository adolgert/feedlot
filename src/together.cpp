
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
#include "boost/mpl/vector.hpp"
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "gspn_random.hpp"
#include "seir_exp.hpp"

namespace smv=afidd::smv;
using namespace smv;


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
  int64_t kind;
  SIRTKey()=default;
  SIRTKey(int64_t i, int64_t j, int64_t k) : i(i), j(j), kind(k) {}

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
    int64_t m=0;
    for (m=3; m<susceptible_cnt_+3; ++m) {
      // Doesn't matter which one we choose. Individuals are identical,
      // so choose first.
      if (lm.template Length<0>(m)==0) {
        lm.template Add<0>(m, IndividualToken{});
        break;
      }
    }
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"IP Fire "<<m);
    // Update the pen summary.
    lm.template Move<0,0>(1, 2, 1);
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
    int64_t m=0;
    for (m=3; m<susceptible_cnt_+3; ++m) {
      if (lm.template Length<0>(m)==0) {
        lm.template Add<0>(m, IndividualToken{});
        break;
      }
    }
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"IF Fire "<<m);
    // Update the pen summary.
    lm.template Move<0,0>(1, 2, 1);
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
    int64_t m=0;
    for (m=3; m<susceptible_cnt_+3; ++m) {
      if (lm.template Length<0>(m)==0) {
        lm.template Add<0>(m, IndividualToken{});
        break;
      }
    }
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"IO Fire "<<m);
    // Update the pen summary.
    lm.template Move<0,0>(1, 2, 1);
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
      double rate=I*s.params.at(SIRParam::Latent);
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
    lm.template Move<0, 0>(0, 2, 1); // Move individual to summary count.
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
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
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
    lm.template Move<0, 0>(0, 2, 1); // Move individual to summary count.
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
  for (int bidx=0; bidx<block_cnt; ++bidx) {
    int base=bidx*row_cnt*2;
    add_edge(base, base+1, g);
    for (int row_idx=1; row_idx<row_cnt; ++row_idx) {
      add_edge(base+row_idx*2, base+row_idx*2+1, g);
      add_edge(base+row_idx*2, base+(row_idx-1)*2, g);
      add_edge(base+row_idx*2+1,  base+(row_idx-1)*2+1, g);
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
void BuildSystem(SIRGSPN& bg, int64_t individual_cnt, int block_cnt, int row_cnt)
{
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

  auto g=BlockStructure(block_cnt, row_cnt);
  auto pen_cnt=num_vertices(g);
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

  enum : int64_t { none, infect0, infect1, infect2, infectious, recover };

  // 2N
  for (int64_t ind_idx=0; ind_idx<individual_cnt; ++ind_idx) {
    int64_t ind_pen=pen_of(ind_idx, per_pen);
    bg.AddTransition({ind_idx, ind_idx, infectious},
      {Edge{{ind_idx, location, e}, -1},
       Edge{{ind_pen, pen_summary, e},-1}, Edge{{ind_idx, location, i}, 1},
       Edge{{ind_pen, pen_summary, i}, 1}},
      std::unique_ptr<SIRTransition>(new Infectious<SIRTransition>())
      );
    bg.AddTransition({ind_idx, ind_idx, recover},
      {Edge{{ind_idx, location, i}, -1},
       Edge{{ind_pen, pen_summary, i}, -1},
       Edge{{ind_pen, pen_summary, r}, 1}},
      std::unique_ptr<SIRTransition>(new Recover<SIRTransition>())
      );
  }

  // NP
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
          bg.AddTransition({src, s_idx, infect0}, infect_vec,
            std::unique_ptr<SIRTransition>(new InfectPen(per_pen))
            );
        }
      } else if (AdjacentPens(d_idx, s_idx, g)) {
        int64_t src_base=d_idx*per_pen;
        for (int64_t src_idx=0; src_idx<per_pen; ++src_idx) {
          int64_t src=src_base+src_idx;
          infect_vec[0]=Edge{{src, location, i}, -1};
          bg.AddTransition({src, s_idx, infect1}, infect_vec,
            std::unique_ptr<SIRTransition>(new InfectFence(per_pen))
            );
        }
      } else {
        int64_t src_base=d_idx*per_pen;
        for (int64_t src_idx=0; src_idx<per_pen; ++src_idx) {
          int64_t src=src_base+src_idx;
          infect_vec[0]=Edge{{src, location, i}, -1};
          bg.AddTransition({src, s_idx, infect2}, infect_vec,
            std::unique_ptr<SIRTransition>(new InfectOther(per_pen))
            );
        }
      }
    }
  }
}


template<typename GSPN, typename Marking>
bool CheckMarking(const GSPN& gspn, const Marking& marking,
    int64_t individual_cnt, int64_t pen_cnt, int64_t per_pen) {
  BOOST_LOG_TRIVIAL(debug)<<"Checking for exactly one state "<<individual_cnt
      <<" "<<pen_cnt<<" "<<per_pen;
  // Check that each individual is in exactly one state.
  std::vector<int64_t> seir{0,0,0,0};
  int ind_mark_total=0;
  for (int64_t ind_idx=0; ind_idx<individual_cnt; ++ind_idx) {
    for (int64_t state_idx=1; state_idx<3; ++state_idx) {
      auto pid=gspn.PlaceVertex({ind_idx, 0, state_idx});
      auto l=Length<0>(marking, pid);
      seir[state_idx]+=l;
      ind_mark_total+=l;
    }
  }
  for (int64_t cp_idx=0; cp_idx<pen_cnt; ++cp_idx) {
    for (int64_t ps_idx=0; ps_idx<4; ps_idx+=3) {
      auto pen_place=gspn.PlaceVertex({cp_idx, 1, ps_idx});
      ind_mark_total+=Length<0>(marking, pen_place);
      seir[ps_idx]=Length<0>(marking, pen_place);
    }
  }
  BOOST_LOG_TRIVIAL(debug)<<"Total of "<<ind_mark_total<<" tokens";
  BOOST_LOG_TRIVIAL(debug)<<"Marking has "<<seir[0]<<", "<<seir[1]
      <<", "<<seir[2]<<", "<<seir[3];
  assert(ind_mark_total==individual_cnt);

  BOOST_LOG_TRIVIAL(debug)<<"Checking pen counts are correct.";
  // Check that each pen count is correct.
  for (int64_t pen_idx=0; pen_idx<pen_cnt; ++pen_idx) {
    for (int64_t state_idx=1; state_idx<3; ++state_idx) {
      int64_t state_pen_id=gspn.PlaceVertex({pen_idx, 1, state_idx});
      int64_t state_pen_cnt=Length<0>(marking, state_pen_id);

      int64_t state_pen_total=0;
      for (int64_t in_pen=0; in_pen<per_pen; ++in_pen) {
        auto ind_in_pen=gspn.PlaceVertex({per_pen*pen_idx+in_pen, 0, state_idx});
        state_pen_total+=Length<0>(marking, ind_in_pen);
      }
      assert(state_pen_total==state_pen_cnt);
      if (state_pen_total!=state_pen_cnt) {
        return false;
      }
    }
  }
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
  TrajectoryObserver& observer_;
  const GSPN& gspn_;
  StateArray seir_;
  int64_t step_cnt{0};

  SEIROutput(const GSPN& gspn, TrajectoryObserver& observer,
      const std::vector<int64_t>& initial)
  : gspn_(gspn), observer_(observer)
  {
    std::get<0>(seir_)=initial[0];
    std::get<1>(seir_)=initial[1];
    std::get<2>(seir_)=initial[2];
    std::get<3>(seir_)=initial[3];
  };

  enum : int64_t { none, infect0, infect1, infect2, infectious, recover };
  void operator()(const SIRState& state) {
    auto transition=gspn_.VertexTransition(state.last_transition);
    switch (transition.kind) {
      case infect0: // infect
      case infect1:
      case infect2:
        get<0>(seir_)-=1;
        get<1>(seir_)+=1;
        break;
      case infectious:
        get<1>(seir_)-=1;
        get<2>(seir_)+=1;
        break;
      case recover: // recover
        get<2>(seir_)-=1;
        get<3>(seir_)+=1;
        break;
      default:
        assert("unknown transition kind");
        break;
    }

    ++step_cnt;
    observer_.Step({get<0>(seir_), get<1>(seir_), get<2>(seir_), get<3>(seir_),
        state.CurrentTime()});
  }

  void final(const SIRState& state) {
    BOOST_LOG_TRIVIAL(info) << "Took "<< step_cnt << " transitions.";
  }
};


int64_t SEIR_run(double end_time, const std::vector<int64_t>& seir_cnt,
    const std::vector<TypedParameter<SIRParam>>& parameters,
    TrajectoryObserver& observer,
    RandGen& rng, int block_cnt, int row_cnt)
{
  int64_t individual_cnt=std::accumulate(seir_cnt.begin(), seir_cnt.end(),
    int64_t{0});

  // The goal is to put all latent and infecteds in the same pen.
  int64_t pen_cnt=2*block_cnt*row_cnt;
  int64_t animals_per_pen=individual_cnt/pen_cnt;
    // 2N 4P 2N NP
  int64_t N=individual_cnt;
  int64_t P=2*block_cnt*row_cnt;
  int64_t guess_cnt=2*N + 4*P + 2*N + N*P;
  SIRGSPN gspn(guess_cnt);
  BuildSystem(gspn, individual_cnt, block_cnt, row_cnt);
  BOOST_LOG_TRIVIAL(debug)<<"GSPN vertex count "<<gspn.VerticesUsed()
      << " guess " << guess_cnt;
  BOOST_LOG_TRIVIAL(info)<<"Created GSPN";
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
    Add<0>(state.marking, summary_id, IndividualToken{});
  }

  //CheckMarking(gspn, state.marking, individual_cnt, pen_cnt, animals_per_pen);

  BOOST_LOG_TRIVIAL(debug)<<"Moving susceptibles to other states.";
  int64_t infected_pen=smv::uniform_index(rng, pen_cnt);
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

  //CheckMarking(gspn, state.marking, individual_cnt, pen_cnt, animals_per_pen);

  //using Propagator=PropagateCompetingProcesses<int64_t,RandGen>;
  using Propagator=NonHomogeneousPoissonProcesses<int64_t,RandGen>;
  Propagator competing;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RandGen>;
  Dynamics dynamics(gspn, {&competing});

  SEIROutput<SIRGSPN,SIRState> output_function(gspn, observer, seir_cnt);

  dynamics.Initialize(&state, &rng);

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
      output_function(state);
    }
  }
  if (running) {
    BOOST_LOG_TRIVIAL(info)<<"Reached end time "<<state.CurrentTime();
  } else {
    BOOST_LOG_TRIVIAL(info)<<"No transitions left to fire at time "<<last_time;
  }
  output_function.final(state);
  return 0;
}

