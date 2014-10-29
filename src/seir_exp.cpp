
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
#include "rider_enums.hpp"
#include "place_transition.hpp"
#include "pen.hpp"

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
using GammaDist=GammaDistribution<RandGen>;
using WeibullDist=WeibullDistribution<RandGen>;


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
      double th=1.0/s.params.at(SIRParam::LatentBeta);
      double a=s.params.at(SIRParam::LatentAlpha);
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<WeibullDistribution<RandGen>>(
        new WeibullDistribution<RandGen>(th, a, te))};
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
class InfectiousGamma : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      double a=s.params.at(SIRParam::LatentAlpha);
      double b=s.params.at(SIRParam::LatentBeta);
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<GammaDistribution<RandGen>>(
        new GammaDistribution<RandGen>(a, b, te))};
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
class InfectiousExponential : public SIRTransition
{
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
      double a=s.params.at(SIRParam::GammaAlpha);
      double b=s.params.at(SIRParam::GammaBeta);
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<GammaDistribution<RandGen>>(
        new GammaDistribution<RandGen>(a, b, te))};
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
    //SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Fire recover " << lm);
    lm.template Move<0, 0>(0, 1, 1);
  }
};


// Now make specific transitions.
class SubClinical : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    int64_t N=lm.template Length<0>(0);
    int64_t I=lm.template Length<0>(1);
    if (N>0 && I>0) {
      double a=s.params.at(SIRParam::SubClinicalAlpha);
      double b=s.params.at(SIRParam::SubClinicalBeta);
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<GammaDistribution<RandGen>>(
        new GammaDistribution<RandGen>(a, b, te))};
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

// The GSPN itself.
using SIRGSPN=
    ExplicitTransitions<SIRPlace, SIRTKey, Local, RandGen, WithParams>;

/*! SIR infection on an all-to-all graph of uncolored tokens.
 */
SIRGSPN
BuildSystem(int64_t individual_cnt, const std::map<ModelOptions,bool>& opts,
  const PenContactGraph& pen_graph)
{
  BuildGraph<SIRGSPN> bg;
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

  auto pen_cnt=num_vertices(pen_graph);
  int64_t per_pen=individual_cnt/pen_cnt;
  assert((individual_cnt % per_pen)==0);
  BOOST_LOG_TRIVIAL(info) << "individuals "<<individual_cnt<<" pens "<<pen_cnt
      <<" per pen "<<per_pen;
  enum : int64_t { s, e, i, r, n, c };

  const int64_t location=0;
  for (int64_t ap_idx=0; ap_idx<individual_cnt; ++ap_idx) {
    for (int64_t place : std::vector<int64_t>{s, e, i, r, n, c}) {
      bg.AddPlace({ap_idx, location, place}, 0);
    }
  }

  for (int64_t ind_idx=0; ind_idx<individual_cnt; ++ind_idx) {
    std::unique_ptr<SIRTransition> infectious;
    std::unique_ptr<SIRTransition> recover;
    if (opts.at(ModelOptions::ExponentialTransitions)) {
      infectious.reset(new InfectiousExponential());
      recover.reset(new RecoverExponential());
    } else if (opts.at(ModelOptions::DoubleGamma)) {
      infectious.reset(new InfectiousGamma());
      recover.reset(new Recover());
    } else {
      infectious.reset(new Infectious());
      recover.reset(new Recover());
    }
    bg.AddTransition({ind_idx, ind_idx, TransitionType::infectious},
      {Edge{{ind_idx, location, e}, -1}, Edge{{ind_idx, location, i}, 1}},
      std::move(infectious)
      );
    bg.AddTransition({ind_idx, ind_idx, TransitionType::recover},
      {Edge{{ind_idx, location, i}, -1}, Edge{{ind_idx, location, r}, 1}},
      std::move(recover)
      );
    bg.AddTransition({ind_idx, ind_idx, TransitionType::subclinical},
      {Edge{{ind_idx, location, n}, -1}, Edge{{ind_idx, location, i}, -1},
       Edge{{ind_idx, location, c}, 1}},
      std::unique_ptr<SIRTransition>(new SubClinical())
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
        } else if (AdjacentPens(pen_s, pen_d, pen_graph)) {
          bg.AddTransition({d_idx, s_idx, TransitionType::infect1},
            {Edge{{d_idx, location, s}, -1}, Edge{{s_idx, location, i}, -1}, 
                Edge{{d_idx, location, e}, 1}, Edge{{s_idx, location, i}, 1}},
            std::unique_ptr<SIRTransition>(new InfectFence())
            );
        } else {
          bg.AddTransition({d_idx, s_idx, TransitionType::infect2},
            {Edge{{d_idx, location, s}, -1}, Edge{{s_idx, location, i}, -1}, 
                Edge{{d_idx, location, e}, 1}, Edge{{s_idx, location, i}, 1}},
            std::unique_ptr<SIRTransition>(new Infect())
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
  std::shared_ptr<PenTrajectoryObserver> observer_;
  const GSPN& gspn_;
  StateArray seir_;
  int64_t step_cnt{0};
  int64_t pen_cnt_;
  int64_t per_pen_;

  SEIROutput(const GSPN& gspn,
      std::shared_ptr<PenTrajectoryObserver> observer,
      const std::vector<int64_t>& initial, int64_t pen_cnt,
      int64_t per_pen)
  : gspn_(gspn), observer_(observer), pen_cnt_(pen_cnt), per_pen_(per_pen)
  {
    std::get<0>(seir_)=initial[0];
    std::get<1>(seir_)=initial[1];
    std::get<2>(seir_)=initial[2];
    std::get<3>(seir_)=initial[3];
  };

  void operator()(const SIRState& state) {
    if (step_cnt==0) {
      this->initial(state);
    }

    auto transition=gspn_.VertexTransition(state.last_transition);
    auto individual=transition.i;
    auto pen=pen_of(individual, per_pen_);
    int64_t compartment=-1;
    switch (transition.kind) {
      case TransitionType::infect0: // infect
      case TransitionType::infect1:
      case TransitionType::infect2:
        get<0>(seir_)-=1;
        get<1>(seir_)+=1;
        compartment=0;
        break;
      case TransitionType::infectious:
        get<1>(seir_)-=1;
        get<2>(seir_)+=1;
        compartment=1;
        break;
      case TransitionType::recover: // recover
        get<2>(seir_)-=1;
        get<3>(seir_)+=1;
        compartment=2;
        break;
      case TransitionType::subclinical:
        compartment=-1;
        break;
      default:
        assert("unknown transition kind");
        break;
    }

    ++step_cnt;
    if (compartment>=0) {
      observer_->Step({individual, pen, compartment, state.CurrentTime()});
    }
  }

  void initial(const SIRState& state) {
    const int64_t pen_summary=1;
    std::vector<TrajectoryEntry> seir_init_(pen_cnt_);
    for (int64_t pen_idx=0; pen_idx<pen_cnt_; ++pen_idx) {
      int64_t splace=gspn_.PlaceVertex({pen_idx, pen_summary, 0});
      seir_init_[pen_idx].s=Length<0>(state.marking, splace);
      splace=gspn_.PlaceVertex({pen_idx, pen_summary, 1});
      seir_init_[pen_idx].e=Length<0>(state.marking, splace);
      splace=gspn_.PlaceVertex({pen_idx, pen_summary, 2});
      seir_init_[pen_idx].i=Length<0>(state.marking, splace);
      splace=gspn_.PlaceVertex({pen_idx, pen_summary, 3});
      seir_init_[pen_idx].r=Length<0>(state.marking, splace);
      seir_init_[pen_idx].t=0;
    }
    observer_->SetInitial(seir_init_);
  }

  void final(const SIRState& state) {
    BOOST_LOG_TRIVIAL(info) << "Took "<< step_cnt << " transitions.";
  }
};



int64_t SEIR_run(double end_time, const std::vector<int64_t>& seir_cnt,
    const std::vector<TypedParameter<SIRParam>>& parameters,
    std::map<ModelOptions,bool> opts,
    const PenContactGraph& pen_graph,
    std::shared_ptr<PenTrajectoryObserver> observer,
    RandGen& rng)
{
  int64_t individual_cnt=std::accumulate(seir_cnt.begin(), seir_cnt.end(),
    int64_t{0});
  auto gspn=BuildSystem(individual_cnt, opts, pen_graph);

  // Marking of the net.
  static_assert(std::is_same<int64_t,SIRGSPN::PlaceKey>::value,
    "The GSPN's internal place type is int64_t.");
  using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<IndividualToken>>;
  using SIRState=GSPNState<Mark,SIRGSPN::TransitionKey,WithParams>;

  int64_t pen_cnt=num_vertices(pen_graph);
  int64_t animals_per_pen=individual_cnt/pen_cnt;

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

  enum : int64_t { s, e, i, r, n, c };
  const int64_t location=0;
  for (int64_t sus_idx=0; sus_idx<individual_cnt; ++sus_idx) {
    auto place_id=gspn.PlaceVertex({sus_idx, location, s});
    Add<0>(state.marking, place_id, IndividualToken{});
  }
  for (int64_t nc_idx=0; nc_idx<individual_cnt; ++nc_idx) {
    auto nonclinical=gspn.PlaceVertex({nc_idx, location, n});
    Add<0>(state.marking, nonclinical, IndividualToken{});
  }
  // The goal is to put all latent and infecteds in the same pen.
  int64_t infected_pen=smv::uniform_index(rng, pen_cnt);
  int64_t first_in_pen=animals_per_pen*infected_pen;
  for (int reinit=1; reinit<4; ++reinit) {
    for (int64_t mv_idx=0; mv_idx<seir_cnt[reinit]; ++mv_idx) {
      auto sus_id=gspn.PlaceVertex({first_in_pen, location, s});
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

  SEIROutput<SIRGSPN,SIRState> output_function(gspn, observer, seir_cnt,
      pen_cnt, animals_per_pen);

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

