#define BOOST_LOG_DYN_LINK 1
#include "boost/log/core.hpp"
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "figtree.h"
#include "hdf_file.hpp"

namespace smv=afidd::smv;
using namespace smv;

int main(int argc, char* argv[]) {
  namespace po=boost::program_options;
  po::options_description desc("Generate ensemble plot from datasets.");

  std::string log_level{"debug"};


  afidd::LogInit(log_level);


  HDFFile file("rider.h5");
  file.OpenRead();

  auto traj_names=file.Trajectories();
  for (auto traj_name : traj_names) {
    BOOST_LOG_TRIVIAL(debug)<<"Loading file "<<traj_name;
    auto trajectory=file.LoadTrajectoryFromPens(traj_name);
  }

  file.Close();
  return 0;
}
