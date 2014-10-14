
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
    int64_t S=0;
    for (int64_t s=1; s<susceptible_cnt_+1; ++s) {
      S+=lm.template Length<0>(s);
    }
    double rate=S*s.params.at(SIRParam::Beta0);
    if (S>0 && I>0 && rate>0.0) {
      BOOST_LOG_TRIVIAL(debug)<<"IP Enabled S "<<S<<" I "<<I<<" r "<<rate;
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
    int64_t S=0;
    for (int64_t s=1; s<susceptible_cnt_+1; ++s) {
      S+=lm.template Length<0>(s);
    }
    int64_t move_it=smv::uniform_index(rng, S);
    int64_t M=0;
    for (int64_t m=1; m<susceptible_cnt_+1; ++m) {
      auto len=lm.template Length<0>(m);
      if (M>=move_it && len>0) {
        lm.template Move<0,0>(m, m+susceptible_cnt_, 1);
      }
      M+=len;
    }
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
    int64_t S=0;
    for (int64_t s=1; s<susceptible_cnt_+1; ++s) {
      S+=lm.template Length<0>(s);
    }
    double rate=S*s.params.at(SIRParam::Beta1);
    if (S>0 && I>0 && rate>0.0) {
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
    int64_t S=0;
    for (int64_t s=1; s<susceptible_cnt_+1; ++s) {
      S+=lm.template Length<0>(s);
    }
    int64_t move_it=smv::uniform_index(rng, S);
    int64_t M=0;
    for (int64_t m=1; m<susceptible_cnt_+1; ++m) {
      auto len=lm.template Length<0>(m);
      if (M>=move_it && len>0) {
        lm.template Move<0,0>(m, m+susceptible_cnt_, 1);
      }
      M+=len;
    }
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
    int64_t S=0;
    for (int64_t s=1; s<susceptible_cnt_+1; ++s) {
      S+=lm.template Length<0>(s);
    }
    double rate=S*s.params.at(SIRParam::Beta2);
    if (S>0 && I>0 && rate>0.0) {
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
    int64_t S=0;
    for (int64_t s=1; s<susceptible_cnt_+1; ++s) {
      S+=lm.template Length<0>(s);
    }
    int64_t move_it=smv::uniform_index(rng, S);
    int64_t M=0;
    for (int64_t m=1; m<susceptible_cnt_+1; ++m) {
      auto len=lm.template Length<0>(m);
      if (M>=move_it && len>0) {
        lm.template Move<0,0>(m, m+susceptible_cnt_, 1);
      }
      M+=len;
    }
  }
};


// Now make specific transitions.
class Infectious : public SIRTransition
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
    //SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Fire recover " << lm);
    lm.template Move<0, 0>(0, 1, 1);
  }
};

// Now make specific transitions.
class Recover : public SIRTransition
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
    //SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Fire recover " << lm);
    lm.template Move<0, 0>(0, 1, 1);
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


// The GSPN itself.
using SIRGSPN=
    ExplicitTransitions<SIRPlace, SIRTKey, Local, RandGen, WithParams>;

/*! SIR infection on an all-to-all graph of uncolored tokens.
 */
SIRGSPN
BuildSystem(int64_t individual_cnt, int block_cnt, int row_cnt)
{
  BuildGraph<SIRGSPN> bg;
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

  auto g=BlockStructure(block_cnt, row_cnt);
  auto pen_cnt=num_vertices(g);
  int64_t per_pen=individual_cnt/pen_cnt;
  assert((individual_cnt % per_pen)==0);
  BOOST_LOG_TRIVIAL(info) << "individuals "<<individual_cnt<<" pens "<<pen_cnt
      <<" per pen "<<per_pen;
  enum : int64_t { s, e, i, r };

  const int64_t location=0;
  for (int64_t ap_idx=0; ap_idx<individual_cnt; ++ap_idx) {
    for (int64_t place : std::vector<int64_t>{s, e, i, r}) {
      bg.AddPlace({ap_idx, location, place}, 0);
    }
  }

  enum : int64_t { none, infect0, infect1, infect2, infectious, recover };

  for (int64_t ind_idx=0; ind_idx<individual_cnt; ++ind_idx) {
    bg.AddTransition({ind_idx, ind_idx, infectious},
      {Edge{{ind_idx, location, e}, -1}, Edge{{ind_idx, location, i}, 1}},
      std::unique_ptr<SIRTransition>(new Infectious())
      );
    bg.AddTransition({ind_idx, ind_idx, recover},
      {Edge{{ind_idx, location, i}, -1}, Edge{{ind_idx, location, r}, 1}},
      std::unique_ptr<SIRTransition>(new Recover())
      );
  }

  std::vector<Edge> infect_vec(1+2*per_pen);
  for (int64_t d_idx=0; d_idx<pen_cnt; ++d_idx) {
    for (int64_t s_idx=0; s_idx<pen_cnt; ++s_idx) {
      int64_t s_base=s_idx*per_pen;
      for (int64_t targ_idx=0; targ_idx<per_pen; ++targ_idx) {
        infect_vec[1+targ_idx]=Edge{{s_base+targ_idx, location, s},-1};
        infect_vec[1+per_pen+targ_idx]=Edge{{s_base+targ_idx, location, e},1};
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
  // std::move the transitions because they contain unique_ptr.
  return std::move(bg.Build());
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
    const std::vector<Parameter>& parameters, TrajectoryObserver& observer,
    RandGen& rng, int block_cnt, int row_cnt)
{
  int64_t individual_cnt=std::accumulate(seir_cnt.begin(), seir_cnt.end(),
    int64_t{0});
  auto gspn=BuildSystem(individual_cnt, block_cnt, row_cnt);

  // Marking of the net.
  static_assert(std::is_same<int64_t,SIRGSPN::PlaceKey>::value,
    "The GSPN's internal place type is int64_t.");
  using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<IndividualToken>>;
  using SIRState=GSPNState<Mark,SIRGSPN::TransitionKey,WithParams>;

  SIRState state;
  for (auto& cp : parameters) {
    state.user.params[cp.kind]=cp.value;
  }

  const int64_t location=0;
  for (int64_t sir_idx=0; sir_idx<4; ++sir_idx) {
    for (int64_t sus_idx=0; sus_idx<seir_cnt[sir_idx]; ++sus_idx) {
      auto place_id=gspn.PlaceVertex({sus_idx, location, sir_idx});
      Add<0>(state.marking, place_id, IndividualToken{});
    }
  }
  // The goal is to put all latent and infecteds in the same pen.
  int64_t pen_cnt=2*block_cnt*row_cnt;
  int64_t animals_per_pen=individual_cnt/pen_cnt;
  int64_t infected_pen=smv::uniform_index(rng, pen_cnt);
  int64_t first_in_pen=animals_per_pen*infected_pen;
  for (int reinit=1; reinit<4; ++reinit) {
    for (int64_t mv_idx=0; mv_idx<seir_cnt[reinit]; ++mv_idx) {
      auto sus_id=gspn.PlaceVertex({first_in_pen, location, 0});
      auto to_id=gspn.PlaceVertex({first_in_pen, location, reinit});
      Move<0,0>(state.marking, sus_id, to_id, 1);
      ++first_in_pen;
    }
  }

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

