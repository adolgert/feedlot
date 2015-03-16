#include <string>
#include <sstream>
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "boost/random/mersenne_twister.hpp"
//#include "mt19937.hpp"
#include "rider_enums.hpp"
#include "feedlot_version.hpp"
#include "feedlot_places.hpp"
#include "feedlot_individual_transitions.hpp"
#include "fmdv.hpp"


namespace smv=afidd::smv;
using namespace smv;

//using RandGen=afidd::rng::mt19937;
using RandGen=boost::mt19937;

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
// The GSPN itself.
using SIRGSPN=
    ExplicitTransitions<SIRPlace, SIRTKey, Local, RandGen, WithParams>;

using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<IndividualToken>>;
using SIRState=GSPNState<Mark,SIRGSPN::TransitionKey,WithParams>;

// Now make specific transitions.
template<typename BaseTransition>
class Reset : public BaseTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const typename BaseTransition::UserState& s,
    const typename BaseTransition::LocalMarking& lm,
    double te, double t0, typename BaseTransition::RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<ExponentialDistribution<RandGen>>(
        new ExponentialDistribution<RandGen>(1, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<TransitionDistribution< 
        typename BaseTransition::RandGen>>(nullptr)};
    }
  }

  virtual void Fire(typename BaseTransition::UserState& s,
    typename BaseTransition::LocalMarking& lm, double t0,
      typename BaseTransition::RandGen& rng) override {
    SMVLOG(BOOST_LOG_TRIVIAL(debug) << "Fire recover " << lm);
    lm.template Move<0, 0>(0, 2, 1);
    lm.template Move<0, 0>(1, 3, 1);
  }
};


struct IndividualOutput {
  const SIRGSPN& gspn_;
  double time_;
  std::vector<std::array<double,2>> trajectory_;
  int64_t step_;

  IndividualOutput(const SIRGSPN& gspn, size_t step_cnt)
      : gspn_(gspn), time_{0.}, trajectory_{step_cnt}, step_{0} {}

  bool operator()(const SIRState& state) {
    double time_delta=state.CurrentTime()-time_;
    auto transition=gspn_.VertexTransition(state.last_transition);
    if (transition.kind==TransitionType::infectious) {
      std::get<0>(trajectory_[step_])=time_delta;
    } else if (transition.kind==TransitionType::recover) {
      std::get<1>(trajectory_[step_])=time_delta;
      ++step_;
    }
    time_=state.CurrentTime();
    bool running=(step_<trajectory_.size());
    BOOST_LOG_TRIVIAL(debug)<<"step "<<step_<<" trans "<<
        static_cast<int>(transition.kind)<<" run " << running;
    return running;
  }

  void final(const SIRState& state) {}

  friend
  std::ostream& operator<<(std::ostream& os, const IndividualOutput& io) {
    for (int64_t ind=0; ind<io.trajectory_.size(); ++ind) {
      os << std::get<0>(io.trajectory_[ind]) << '\t'
        << std::get<1>(io.trajectory_[ind]) << std::endl;
    }
    return os;
  }
};


void BuildExponentialGSPN(SIRGSPN& gspn) {
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;
  enum : int64_t { s, e, i, r };
  gspn.AddPlace({0, 0, e});
  gspn.AddPlace({0, 0, i});
  gspn.AddPlace({0, 0, r});
  gspn.AddPlace({0, 1, e});
  gspn.AddPlace({0, 1, i});
  gspn.AddPlace({0, 1, r});

  int64_t ind_idx=0;
  int64_t ind_pen=0;
  const int64_t location=0;
  const int64_t pen_summary=1;
  gspn.AddTransition({ind_idx, location, TransitionType::infectious},
    {Edge{{ind_idx, location, e}, -1},
    Edge{{ind_pen, pen_summary, e},-1},
    Edge{{ind_idx, location, i}, 1},
    Edge{{ind_pen, pen_summary, i}, 1}},
    std::unique_ptr<SIRTransition>(new InfectiousExponential<SIRTransition>())
    );
  gspn.AddTransition({ind_idx, ind_idx, TransitionType::recover},
    {Edge{{ind_idx, location, i}, -1},
    Edge{{ind_pen, pen_summary, i}, -1},
       Edge{{ind_idx, location, r}, 1},
     Edge{{ind_pen, pen_summary, r}, 1}},
    std::unique_ptr<SIRTransition>(new RecoverExponential<SIRTransition>())
    );
  gspn.AddTransition({ind_idx, ind_idx, TransitionType::reset},
    {Edge{{ind_idx, location, r}, -1},
    Edge{{ind_pen, pen_summary, r}, -1},
    Edge{{ind_idx, location, e}, 1},
    Edge{{ind_pen, pen_summary, e}, 1}},
    std::unique_ptr<SIRTransition>(new Reset<SIRTransition>())
    );
}

void BuildGSPN(SIRGSPN& gspn) {
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;
  enum : int64_t { s, e, i, r };
  gspn.AddPlace({0, 0, e});
  gspn.AddPlace({0, 0, i});
  gspn.AddPlace({0, 0, r});
  gspn.AddPlace({0, 1, e});
  gspn.AddPlace({0, 1, i});
  gspn.AddPlace({0, 1, r});

  int64_t ind_idx=0;
  int64_t ind_pen=0;
  const int64_t location=0;
  const int64_t pen_summary=1;
  gspn.AddTransition({ind_idx, location, TransitionType::infectious},
    {Edge{{ind_idx, location, e}, -1},
    Edge{{ind_pen, pen_summary, e},-1},
    Edge{{ind_idx, location, i}, 1},
    Edge{{ind_pen, pen_summary, i}, 1}},
    std::unique_ptr<SIRTransition>(new Infectious<SIRTransition>())
    );
  gspn.AddTransition({ind_idx, ind_idx, TransitionType::recover},
    {Edge{{ind_idx, location, i}, -1},
    Edge{{ind_pen, pen_summary, i}, -1},
       Edge{{ind_idx, location, r}, 1},
     Edge{{ind_pen, pen_summary, r}, 1}},
    std::unique_ptr<SIRTransition>(new Recover<SIRTransition>())
    );
  gspn.AddTransition({ind_idx, ind_idx, TransitionType::reset},
    {Edge{{ind_idx, location, r}, -1},
    Edge{{ind_pen, pen_summary, r}, -1},
    Edge{{ind_idx, location, e}, 1},
    Edge{{ind_pen, pen_summary, e}, 1}},
    std::unique_ptr<SIRTransition>(new Reset<SIRTransition>())
    );
}

int main(int argc, char *argv[]) {
  namespace po=boost::program_options;
  po::options_description desc("Test of transitions for an individual.");
  int64_t run_cnt=10;
  bool is_exponential=false;
  std::string log_level;
  desc.add_options()
    ("help", "show help message")
    ("loglevel", po::value<std::string>(&log_level)->default_value("info"),
      "Set the logging level to trace, debug, info, warning, error, or fatal.")
    ("runcnt",
      po::value<int64_t>(&run_cnt)->default_value(run_cnt),
      "number of runs")
    ("exponential",
      po::value<bool>(&is_exponential)->default_value(is_exponential),
      "whether to use all-exponential distributions");

  po::variables_map vm;
  auto parsed_options=po::parse_command_line(argc, argv, desc);
  po::store(parsed_options, vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  afidd::LogInit(log_level);
  RandGen rng{7432};

  SIRGSPN gspn{5};
  std::map<SIRParam,TypedParameter<SIRParam>> parameters;
  if (is_exponential) {
    BuildExponentialGSPN(gspn);
    FMDV_Mardones_Exponential(parameters);
  } else {
    BuildGSPN(gspn);
    FMDV_Mardones_Nonexponential(parameters);
  }

  SIRState state;
  for (auto& cp : parameters) {
    state.user.params[cp.second.kind]=cp.second.value;
  }
  auto e_cattlebeast=gspn.PlaceVertex({0, 0, 1});
  Add<0>(state.marking, e_cattlebeast, IndividualToken{});
  auto e_summary=gspn.PlaceVertex({0, 1, 1});
  Add<0>(state.marking, e_summary, IndividualToken{});

  IndividualOutput output_function{gspn, static_cast<size_t>(run_cnt)};

  //using Propagator=PropagateCompetingProcesses<int64_t,RandGen>;
  using Propagator=NonHomogeneousPoissonProcesses<int64_t,RandGen>;
  Propagator competing;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RandGen>;
  Dynamics dynamics(gspn, {&competing});

  dynamics.Initialize(&state, &rng);

  BOOST_LOG_TRIVIAL(info)<<"Starting main loop";
  bool running=true;
  auto nothing=[](SIRState&)->void {};
  double last_time=state.CurrentTime();
  while (running) {
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
      BOOST_LOG_TRIVIAL(info)<<"no enabled transitions";
    }
  }
  BOOST_LOG_TRIVIAL(info)<<"time "<<last_time;
  output_function.final(state);

  std::cout << output_function;

  return 0;
}
