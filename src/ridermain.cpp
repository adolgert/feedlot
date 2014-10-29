#include <string>
#include <sstream>
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "rider_enums.hpp"
#include "parameter.hpp"
#include "trajectory.hpp"
#include "rider.hpp"
#include "hdf_file.hpp"
#include "ensemble.hpp"
#include "fmdv.hpp"
#include "model_options.hpp"
#include "feedlot_version.hpp"

int main(int argc, char *argv[]) {
  namespace po=boost::program_options;
  po::options_description desc("Feedlot with rider.");
  int64_t individual_cnt=1024;
  int64_t exposed_cnt=1;
  int64_t infected_cnt=0;
  int64_t recovered_cnt=0;

  int64_t block_cnt=2;
  int64_t row_cnt=8;
  int64_t disconnected_pens=0;

  int run_cnt=1;
  size_t rand_seed=1;
  // Time is in years.
  using MyParm=TypedParameter<SIRParam>;
  std::vector<MyParm> parameters;
  parameters.emplace_back(MyParm{SIRParam::Beta0, "beta0", 1/0.26,
    "density-dependent infection rate within a pen"});
  parameters.emplace_back(MyParm{SIRParam::Beta1, "beta1", 0.1/0.26,
    "density-dependent infection rate across a fence"});
  parameters.emplace_back(MyParm{SIRParam::Beta2, "beta2", 0.001/0.26,
    "density-dependent infection rate to any other animal"});
  parameters.emplace_back(MyParm{SIRParam::RiderMove,
    "ridermove", 1.0/(15/(60*24.0)), "rate for rider to move pens"});
  parameters.emplace_back(MyParm{SIRParam::RiderRecover,
    "riderrecover", 24, "rider recover from carrying disease"});
  parameters.emplace_back(MyParm{SIRParam::RiderInfect,
    "riderinfect", 1.0/0.26, "rate for rider to infect in pen"});
  parameters.emplace_back(MyParm{SIRParam::RiderGetInfected,
    "ridergetinfected", 1.0/0.26, "rate for rider to pick up infection"});
  FMDV_Mardones_Nonexponential(parameters);
  FMDV_Mardones_Exponential(parameters);
  auto model_opts=model_options();
  model_opts[ModelOptions::AllToAllInfection]=false;
  model_opts[ModelOptions::Rider]=true;

  double end_time=std::numeric_limits<double>::infinity();
  int thread_cnt=1;
  std::string log_level;
  std::string data_file("rider.h5");
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
    ("penblocks",
      po::value<int64_t>(&block_cnt),
      "number of city blocks of pens")
    ("penrows",
      po::value<int64_t>(&row_cnt),
      "number of rows within a block of pens")
    ("disconnected",
      po::value<int64_t>(&disconnected_pens),
      "How many pens to make with no adjacency.")
    ("seed",
      po::value<size_t>(&rand_seed)->default_value(rand_seed),
      "seed for random number generator")
    ("endtime",
      po::value<double>(&end_time)->default_value(end_time),
      "how many years to run")
    ("exponential",
      po::value<bool>(&model_opts[ModelOptions::ExponentialTransitions])->
        default_value(model_opts[ModelOptions::ExponentialTransitions]),
      "Use exponentially-distributed latent and infectious periods")
    ("doublegamma",
      po::value<bool>(&model_opts[ModelOptions::DoubleGamma])->
        default_value(model_opts[ModelOptions::DoubleGamma]),
      "Use gamma-distributed latent and infectious periods")
    ("rider",
      po::value<bool>(&model_opts[ModelOptions::Rider])->default_value(
        model_opts[ModelOptions::Rider]),
      "Whether a rider carries infection to next pens")
    ("otherpens",
      po::value<bool>(&model_opts[ModelOptions::AllToAllInfection])->
      default_value(model_opts[ModelOptions::AllToAllInfection]),
      "Whether there is a background infection between any two pens.")
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

  PenContactGraph pen_graph;
  if (disconnected_pens>0) {
    pen_graph=DisconnectedPens(disconnected_pens);
  } else {
    pen_graph=BlockStructure(block_cnt, row_cnt);
  }

  HDFFile file(data_file);
  if (!file.Open(!save_file)) {
    BOOST_LOG_TRIVIAL(error)<<"could not open output file: "<<data_file;
    return -1;
  }
  file.WriteExecutableData(compile_info, parsed_options, seir_init);

  auto runnable=[=, &file](RandGen& rng, size_t single_seed, size_t idx)->void {
    std::shared_ptr<PenTrajectorySave> observer=
      std::make_shared<PenTrajectorySave>(static_cast<size_t>(individual_cnt));
    auto trajectory_save=std::make_shared<TrajectorySave>(4*individual_cnt);

    SEIR_run(end_time, seir_init, parameters, model_opts, pen_graph,
        observer, trajectory_save, rng);
    file.SavePenTrajectory(parameters, single_seed, idx, observer->Trajectory(),
      observer->PenInitial());
    file.SaveTrajectory(parameters, single_seed, idx,
        trajectory_save->Trajectory());
  };

  afidd::smv::Ensemble<decltype(runnable),RandGen> ensemble(runnable, thread_cnt,
      run_cnt, rand_seed);
  ensemble.Run();
  BOOST_LOG_TRIVIAL(debug)<<"Finished running ensemble.";

  file.Close();

  return 0;
}
