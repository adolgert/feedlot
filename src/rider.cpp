
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



// Now make specific transitions.
class InfectPen : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    // If these are just size_t, then the rate calculation overflows.
    int64_t S=lm.template Length<0>(0);
    int64_t I=lm.template Length<0>(1);
    double rate=s.params.at(SIRParam::Beta0);
    if (S>0 && I>0 && rate>0.0) {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"infection rate "<<rate<<" beta0 "<<
      //  s.params.at(SIRParam::Beta0) << " beta1 " <<
       // s.params.at(SIRParam::Beta1) << " t0 " << t0 << " N "<<(S+I+R)
       // << " S "<<S <<" I " << I);
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
    lm.template Move<0,0>(0, 2, 1);
  }
};


// Now make specific transitions.
class InfectFence : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    // If these are just size_t, then the rate calculation overflows.
    int64_t S=lm.template Length<0>(0);
    int64_t I=lm.template Length<0>(1);
    double rate=s.params.at(SIRParam::Beta1);
    if (S>0 && I>0 && rate>0.0) {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"infection rate "<<rate<<" beta0 "<<
      //  s.params.at(SIRParam::Beta0) << " beta1 " <<
       // s.params.at(SIRParam::Beta1) << " t0 " << t0 << " N "<<(S+I+R)
       // << " S "<<S <<" I " << I);
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
    lm.template Move<0,0>(0, 2, 1);
  }
};


// Now make specific transitions.
class Infect : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    // If these are just size_t, then the rate calculation overflows.
    int64_t S=lm.template Length<0>(0);
    int64_t I=lm.template Length<0>(1);
    double rate=s.params.at(SIRParam::Beta2);
    if (S>0 && I>0 && rate>0.0) {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"infection rate "<<rate<<" beta0 "<<
      //  s.params.at(SIRParam::Beta0) << " beta1 " <<
       // s.params.at(SIRParam::Beta1) << " t0 " << t0 << " N "<<(S+I+R)
       // << " S "<<S <<" I " << I);
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
    lm.template Move<0,0>(0, 2, 1);
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



// Now make specific transitions.
class MoveRider : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      double rate=I*s.params.at(SIRParam::RiderMove);
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
class RecoverRider : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      double rate=I*s.params.at(SIRParam::RiderRecover);
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
class InfectRider : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t S=lm.template Length<0>(0);
    int64_t I=lm.template Length<0>(1);
    if (S>0 && I>0) {
      double rate=I*s.params.at(SIRParam::RiderInfect);
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
    lm.template Move<0, 0>(0, 2, 1);
  }
};


// Now make specific transitions.
class RiderInfects : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t S=lm.template Length<0>(0);
    int64_t I=lm.template Length<0>(1);
    if (S>0 && I>0) {
      double rate=I*s.params.at(SIRParam::RiderGetInfected);
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
    lm.template Move<0, 0>(0, 2, 1);
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
SIRGSPN
BuildSystem(int64_t individual_cnt, int block_cnt, int row_cnt)
{
  BuildGraph<SIRGSPN> bg;
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

  auto g=BlockStructure(block_cnt, row_cnt);
  if (false) {
    boost::dynamic_properties dp;
    std::ofstream contact_file("blocks.graphml");
    boost::write_graphml(contact_file, g, dp, true);
  }
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
  const int64_t rider=2;
  for (int64_t pp_idx=0; pp_idx<pen_cnt; ++pp_idx) {
    bg.AddPlace({pp_idx, rider, s}, 0);
    bg.AddPlace({pp_idx, rider, i}, 0);
  }

  for (int64_t ind_idx=0; ind_idx<individual_cnt; ++ind_idx) {
    bg.AddTransition({ind_idx, ind_idx, TransitionType::infectious},
      {Edge{{ind_idx, location, e}, -1}, Edge{{ind_idx, location, i}, 1}},
      std::unique_ptr<SIRTransition>(new Infectious())
      );
    bg.AddTransition({ind_idx, ind_idx, TransitionType::recover},
      {Edge{{ind_idx, location, i}, -1}, Edge{{ind_idx, location, r}, 1}},
      std::unique_ptr<SIRTransition>(new Recover())
      );
  }

  for (int64_t d_idx=0; d_idx<individual_cnt; ++d_idx) {
    int pen_d=(d_idx-1)/per_pen;
    for (int64_t s_idx=0; s_idx<individual_cnt; ++s_idx) {
      int pen_s=(s_idx-1)/per_pen;
      if (s_idx!=d_idx) {
        if (pen_d==pen_s) {
          bg.AddTransition({d_idx, s_idx, TransitionType::infect0},
            {Edge{{d_idx, location, s}, -1}, Edge{{s_idx, location, i}, -1}, 
                Edge{{d_idx, location, e}, 1}, Edge{{s_idx, location, i}, 1}},
            std::unique_ptr<SIRTransition>(new InfectPen())
            );
        } else if (AdjacentPens(pen_s, pen_d, g)) {
          bg.AddTransition({d_idx, s_idx, TransitionType::infect1},
            {Edge{{d_idx, location, s}, -1}, Edge{{s_idx, location, i}, -1}, 
                Edge{{d_idx, location, e}, 1}, Edge{{s_idx, location, i}, 1}},
            std::unique_ptr<SIRTransition>(new InfectFence())
            );
        } else {
          ; //No long-range transitions.
        }
      }
    }
  }

  for (int64_t rp_idx=0; rp_idx<pen_cnt; ++rp_idx) {
    // Move
    auto pen_s=SIRPlace{rp_idx, rider, s};
    auto pen_i=SIRPlace{rp_idx, rider, i};
    int64_t next_pen=static_cast<int64_t>((rp_idx+1)%pen_cnt);
    auto pen_s_n=SIRPlace{next_pen, rider, s};
    auto pen_i_n=SIRPlace{next_pen, rider, i};
    bg.AddTransition({rp_idx, rp_idx, TransitionType::movers},
      {Edge{pen_s, -1}, Edge{pen_s_n, 1}},
      std::unique_ptr<SIRTransition>(new MoveRider()));
    bg.AddTransition({rp_idx, rp_idx, TransitionType::moveri},
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
      bg.AddTransition({rp_idx, c_idx, TransitionType::infectbyr},
        {Edge{{c_idx, location, s},-1}, Edge{pen_i, -1}, 
         Edge{{c_idx, location, e}, 1}},
         std::unique_ptr<SIRTransition>(new RiderInfects()));
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
  std::shared_ptr<PenTrajectoryObserver> observer_;
  const GSPN& gspn_;
  StateArray seir_;
  int64_t step_cnt{0};
  int64_t individual_cnt_;
  int64_t per_pen_;

  SEIROutput(const GSPN& gspn, std::shared_ptr<PenTrajectoryObserver> observer,
      const std::vector<int64_t>& initial, int64_t per_pen)
  : gspn_(gspn), observer_(observer), per_pen_(per_pen)
  {
    std::get<0>(seir_)=initial[0];
    std::get<1>(seir_)=initial[1];
    std::get<2>(seir_)=initial[2];
    std::get<3>(seir_)=initial[3];
    individual_cnt_=std::accumulate(initial.begin(), initial.end(), int64_t{0});
  };

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
        break;
    }

    ++step_cnt;
    if (affected) {
      int64_t pen=pen_of(individual, per_pen_);
      observer_->Step({individual, pen, compartment, state.CurrentTime()});
    }
    return (std::get<3>(seir_)<individual_cnt_);
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
  // The rider
  auto rider_id=gspn.PlaceVertex({0, 2, 0});
  Add<0>(state.marking, rider_id, IndividualToken{});

  //using Propagator=PropagateCompetingProcesses<int64_t,RandGen>;
  using Propagator=NonHomogeneousPoissonProcesses<int64_t,RandGen>;
  Propagator competing;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RandGen>;
  Dynamics dynamics(gspn, {&competing});

  SEIROutput<SIRGSPN,SIRState> output_function(gspn, observer, seir_cnt,
    animals_per_pen);

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

