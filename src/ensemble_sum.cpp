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
  auto initial=file.InitialValues();
  int64_t total=std::accumulate(initial.begin(), initial.end(), int64_t{0});
  BOOST_LOG_TRIVIAL(debug)<<"Total animals "<<total;
  int64_t trajectory_cnt{0};
  int64_t event_cnt{0};
  auto evt=file.EventsInFile();
  trajectory_cnt=std::get<0>(evt);
  event_cnt=std::get<1>(evt);
  BOOST_LOG_TRIVIAL(debug)<<"Trajectories "<<trajectory_cnt
      <<" events "<<event_cnt;

  // Because the state at t=0.0 is added to the event count.
  int64_t sample_cnt=trajectory_cnt+event_cnt;
  std::vector<std::vector<double>> sei;
  for (int64_t in_idx=0; in_idx<3; ++in_idx) {
    sei.push_back(std::vector<double>(2*sample_cnt, 0));
  }
  // Pre-scaling before single-variate kernel density estimation.
  double days_per_individual=100.0/total;

  double max_time{0.0};
  int64_t entry_idx=0;
  auto traj_names=file.Trajectories();
  for (auto traj_name : traj_names) {
    BOOST_LOG_TRIVIAL(debug)<<"Loading file "<<traj_name;
    auto trajectory=file.LoadTrajectoryFromPens(traj_name);
    for (const auto& entry : trajectory) {
      sei[0][2*entry_idx]=days_per_individual*entry.s;
      sei[1][2*entry_idx]=days_per_individual*entry.e;
      sei[2][2*entry_idx]=days_per_individual*entry.i;
      for (int64_t st_idx=0; st_idx<3; ++st_idx) {
        sei[st_idx][2*entry_idx+1]=entry.t;
      }
      if (entry.t>max_time) max_time=entry.t;
      ++entry_idx;
    }
  }

  int hres=100;
  int vres=100;
  int target_cnt=hres*vres;
  std::vector<double> target(2*target_cnt, 0.0);
  for (int vidx=0; vidx<vres; ++vidx) {
    for (int hidx=0; hidx<hres; ++hidx) {
      target[2*(vidx*hres+hidx)]=hidx*(max_time/hres);
      target[2*(vidx*hres+hidx)+1]=vidx*(total/static_cast<double>(vres));
    }
  }
  std::vector<double> weight(target_cnt, 1.0);
  int weight_kind_cnt=1;

  int sample_dimension=2;
  int source_cnt=sample_cnt;
  double h=4.0; // About 100x100 is domain.
  double epsilon=1e-2;

  std::vector<double> interpolant(weight_kind_cnt*source_cnt, 0.0);
  figtree(sample_dimension, source_cnt, target_cnt, weight_kind_cnt,
      &sei[0][0], h, &weight[0], &target[0], epsilon, &interpolant[0]);

  file.Close();
  return 0;
}
