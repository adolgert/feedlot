#include <string>
#include <sstream>
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "together.hpp"
#include "hdf_file.hpp"
#include "ensemble.hpp"
#include "feedlot_version.hpp"




/*! Save the whole trajectory.
 */
class TrajectorySave : public TrajectoryObserver
{
  std::vector<TrajectoryEntry> trajectory_;
 public:
  virtual void Step(TrajectoryEntry seirt) override {
    trajectory_.emplace_back(seirt);
  }
  virtual const std::vector<TrajectoryEntry>& Trajectory() const {
    return trajectory_; }
};



/*! Save the trajectory every time any of SIR change by a percent.
 */
class PercentTrajectorySave : public TrajectoryObserver
{
  int64_t step_{0};
  int64_t threshhold_{0};
  double percent_{0.0001};

  TrajectoryEntry last_{0,0,0,0,0.0};
  std::vector<TrajectoryEntry> trajectory_;
 public:
  PercentTrajectorySave() {}

  virtual void Step(TrajectoryEntry seirt) override {
    if (0==step_) {
      last_=seirt;
      threshhold_=std::floor(percent_*(seirt.s+seirt.e+seirt.i+seirt.r));
      trajectory_.emplace_back(seirt);
    } else {
      bool ps=std::abs(seirt.s-last_.s)>threshhold_;
      bool pe=std::abs(seirt.e-last_.e)>threshhold_;
      bool pi=std::abs(seirt.i-last_.i)>threshhold_;
      bool pr=std::abs(seirt.r-last_.r)>threshhold_;
      if (ps||pe||pi||pr) {
        trajectory_.emplace_back(seirt);
        last_=seirt;
      }
    }
    ++step_;
  }
  virtual const std::vector<TrajectoryEntry>& Trajectory() const {
    return trajectory_; }
};




int main(int argc, char *argv[]) {
  namespace po=boost::program_options;
  po::options_description desc("Well-mixed SIR with demographics.");
  int64_t individual_cnt=1024;
  int64_t exposed_cnt=1;
  int64_t infected_cnt=0;
  int64_t recovered_cnt=0;

  int run_cnt=1;
  size_t rand_seed=1;
  // Time is in years.
  std::vector<Parameter> parameters;
  parameters.emplace_back(Parameter{SIRParam::Beta0, "beta0", 1/0.26,
    "density-dependent infection rate within a pen"});
  parameters.emplace_back(Parameter{SIRParam::Beta1, "beta1", 0.1/0.26,
    "density-dependent infection rate across a fence"});
  parameters.emplace_back(Parameter{SIRParam::Beta2, "beta2", 0.001/0.26,
    "density-dependent infection rate to any other animal"});
  parameters.emplace_back(Parameter{SIRParam::Gamma, "gamma", 1/8.0,
    "recovery rate"});
  double end_time=30.0;
  bool exacttraj=true;
  bool exactinfect=false;
  int thread_cnt=1;
  std::string log_level;
  std::string data_file("wellmixed.h5");
  bool save_file=false;
  std::string translation_file;
  bool test=false;

  desc.add_options()
    ("help", "show help message")
    ("threadcnt,j",
      po::value<int>(&thread_cnt)->default_value(thread_cnt),
      "number of threads")
    ("runcnt",
      po::value<int>(&run_cnt)->default_value(run_cnt),
      "number of runs")
    ("size,s",
      po::value<int64_t>(&individual_cnt)->default_value(individual_cnt),
      "size of the population")
    ("exposed,e",
      po::value<int64_t>(&exposed_cnt),
      "number of exposed")
    ("infected,i",
      po::value<int64_t>(&infected_cnt),
      "number of infected")
    ("recovered,r",
      po::value<int64_t>(&recovered_cnt),
      "number of recovered")
    ("seed",
      po::value<size_t>(&rand_seed)->default_value(rand_seed),
      "seed for random number generator")
    ("endtime",
      po::value<double>(&end_time)->default_value(end_time),
      "how many years to run")
    ("exacttraj",
      po::value<bool>(&exacttraj)->default_value(exacttraj),
      "save trajectory only when it changes by a certain amount")
    ("exactinfect",
      po::value<bool>(&exactinfect)->default_value(exactinfect),
      "set true to use exact distribution for seasonal infection")
    ("datafile",
      po::value<std::string>(&data_file)->default_value(data_file),
      "Write to this data file.")
    ("save",
      po::value<bool>(&save_file)->default_value(save_file),
      "Add data to file instead of erasing it with new data.")
    ("loglevel", po::value<std::string>(&log_level)->default_value("info"),
      "Set the logging level to trace, debug, info, warning, error, or fatal.")
    ("translate",
      po::value<std::string>(&translation_file)->default_value(""),
      "write file relating place ids to internal ids")
    ("info", "show provenance of program")
    ;

  for (auto& p : parameters) {
    desc.add_options()(p.name.c_str(),
      po::value<double>(&p.value)->default_value(p.value),
      p.description.c_str());
  }

  po::variables_map vm;
  auto parsed_options=po::parse_command_line(argc, argv, desc);
  po::store(parsed_options, vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  std::map<std::string,std::string> compile_info {
    {"VERSION", VERSION}, {"COMPILETIME", COMPILETIME},
    {"CONFIG", CFG}
  };
  if (vm.count("info")) {
    for (auto& kv : compile_info) {
      std::cout << kv.second << "\n\n";
    }
    return 0;
  }

  afidd::LogInit(log_level);

  if (test) {
    ;
  }

  std::map<SIRParam,double*> params;
  for (auto& pm : parameters) {
    assert(params.find(pm.kind)==params.end());
    params[pm.kind]=&pm.value;
  }

  int64_t susceptible_cnt=individual_cnt-(exposed_cnt+infected_cnt+recovered_cnt);

  assert(susceptible_cnt>0);
  if (susceptible_cnt<0) {
    BOOST_LOG_TRIVIAL(error)<<"Number of susceptibles is "<<susceptible_cnt;
    return -2;
  }
  std::vector<int64_t> seir_init{susceptible_cnt, exposed_cnt, infected_cnt,
        recovered_cnt};
  BOOST_LOG_TRIVIAL(info)<<"Starting with sir="<<seir_init[0]<<" "<<seir_init[1]
    <<" "<<seir_init[2] << " "<<seir_init[3];

  for (auto& showp : parameters) {
    BOOST_LOG_TRIVIAL(info)<<showp.name<<" "<<showp.value;
  }

  HDFFile file(data_file);
  if (!file.Open(!save_file)) {
    BOOST_LOG_TRIVIAL(error)<<"could not open output file: "<<data_file;
    return -1;
  }
  file.WriteExecutableData(compile_info, parsed_options, seir_init);

  auto runnable=[=](RandGen& rng, size_t single_seed, size_t idx)->void {
    std::shared_ptr<TrajectoryObserver> observer=0;
    if (exacttraj) {
      observer=std::make_shared<TrajectorySave>();
    } else {
      observer=std::make_shared<PercentTrajectorySave>();
    }

    SEIR_run(end_time, seir_init, parameters, *observer, rng, 2, 4);
    file.SaveTrajectory(parameters, single_seed, idx, observer->Trajectory());
  };

  afidd::smv::Ensemble<decltype(runnable),RandGen> ensemble(runnable, thread_cnt,
      run_cnt, rand_seed);
  ensemble.Run();
  BOOST_LOG_TRIVIAL(debug)<<"Finished running ensemble.";

  file.Close();

  return 0;
}
