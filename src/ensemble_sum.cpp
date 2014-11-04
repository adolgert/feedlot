#define BOOST_LOG_DYN_LINK 1
#include "boost/log/core.hpp"
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "figtree.h"
#include "hdf_file.hpp"

namespace smv=afidd::smv;
using namespace smv;

void watermark_for_orientation(std::vector<double>& interpolant,
    int hres, int vres) {
  for (int vwidx=0; vwidx<vres/10; ++vwidx) {
    for (int hwidx=0; hwidx<hres/10; ++hwidx) {
      interpolant[vwidx*hres+hwidx]=1.0;
    }
  }
  for (int vsidx=vres/10; vsidx<vres; ++vsidx) {
    for (int hsidx=0; hsidx<hres/10; ++hsidx) {
      interpolant[vsidx*hres+hsidx]=0.5;
    }
  }
}

class ClinicalGreaterThan {
 public:
  ClinicalGreaterThan(double fraction, int64_t total, int64_t trajectory_cnt)
  : found_{false}, when_(trajectory_cnt), idx_{0} {
    val_=static_cast<int64_t>(std::ceil(total*fraction));
  }
  void observe(int64_t const* seirc, double time) {
    if (!found_ && seirc[4]>=val_) {
      when_[idx_]=time;
    }
  }
  void done_trajectory() {
    idx_++;
  }
 private:
  int64_t val_;
  bool found_;
  size_t idx_;
  std::vector<double> when_;
};


class TrajectoryDensity {
 public:
  TrajectoryDensity(int64_t total_individuals, double max_time,
      int64_t vres, int64_t hres, int64_t total_event_cnt)
  : total_individuals_{total_individuals+1}, max_time_{max_time+0.1},
    vres_{vres}, hres_{hres}, data_(hres*vres, 0.0),
    scale_(1.0/total_event_cnt) {
  }
  void observe(int64_t const* seirc, double time) {
    int64_t infecteds=seirc[1]+seirc[2];
    int64_t x=static_cast<int64_t>(time*hres_/max_time_);
    int64_t y=static_cast<int64_t>(
        infecteds*static_cast<double>(vres_)/total_individuals_);
    data_[y*hres_+x]+=scale_;
  }
  void done_trajectory() {
  }
  void report() {
    // The x value runs faster. time is the x.
    std::vector<double> x(hres_);
    for (int xidx=0; xidx<hres_; ++xidx) {
      x[xidx]=xidx*(max_time_/hres_);
    }
    std::vector<double> y(vres_);
    for (int yidx=0; yidx<vres_; ++yidx) {
      // Don't rescale our reported axes.
      y[yidx]=yidx*(total_individuals_/static_cast<double>(vres_));
    }
  }
 private:
  std::vector<double> data_;
  int64_t total_individuals_;
  double max_time_;
  int64_t vres_;
  int64_t hres_;
  double scale_;
};


class TotalInfected {
 public:
  TotalInfected(int64_t trajectory_cnt) : count(trajectory_cnt), idx_{0} {}
  void observe(const std::vector<int64_t>& seirc,
      const std::vector<double>& time, int64_t entry_cnt) {
    count[idx_]=seirc[5*(entry_cnt-1)+3];
    ++idx_;
  }
 private:
  std::vector<int64_t> count;
  size_t idx_;
};


int main(int argc, char* argv[]) {
  std::string filename("rider.h5");
  std::string outfilename("image.h5");
  int64_t hres=500;
  int64_t vres=500;

  namespace po=boost::program_options;
  po::options_description desc("Generate ensemble plot from datasets.");

  std::string log_level{"debug"};


  afidd::LogInit(log_level);


  HDFFile file(filename);
  file.OpenRead();
  auto initial=file.InitialValues();
  int64_t total_individuals=
      std::accumulate(initial.begin(), initial.end(), int64_t{0});
  BOOST_LOG_TRIVIAL(debug)<<"Total animals "<<total_individuals;
  int64_t trajectory_cnt{0};
  int64_t event_cnt{0};
  auto evt=file.EventsInFile();
  trajectory_cnt=std::get<0>(evt);
  event_cnt=std::get<1>(evt);
  BOOST_LOG_TRIVIAL(debug)<<"Trajectories "<<trajectory_cnt
      <<" events "<<event_cnt;

  auto end=file.EndTimes();
  const auto& times=std::get<0>(end);
  const auto& events=std::get<1>(end);
  BOOST_LOG_TRIVIAL(info)<<"Found end times "<<std::get<0>(end).size();
  for (int et_idx=0; et_idx<std::get<0>(end).size(); ++et_idx) {
    std::cout << std::get<0>(end)[et_idx] << '\t' <<
        std::get<1>(end)[et_idx] << std::endl;
  }
  double max_tr_time=*std::max_element(times.begin(), times.end());
  int64_t max_event=*std::max_element(events.begin(), events.end());
  int64_t total_event_cnt=std::accumulate(events.begin(), events.end(),
    int64_t{0});
  BOOST_LOG_TRIVIAL(info)<<"max time "<<max_tr_time<<" max event "<<max_event
      <<" total events "<<total_event_cnt;

  int compartment_cnt=5;
  std::vector<int64_t> seirc(max_event*compartment_cnt);
  std::vector<double> time(max_event);

  ClinicalGreaterThan one_percent_clinical(0.01, total_individuals,
      trajectory_cnt);
  ClinicalGreaterThan five_percent_clinical(0.01, total_individuals,
      trajectory_cnt);
  TrajectoryDensity density(total_individuals, max_tr_time, vres, hres,
      total_event_cnt);
  TotalInfected total_infected(trajectory_cnt);

  std::string max_trajectory_name;
  int64_t max_trajectory_len{0};

  auto traj_names=file.Trajectories();
  for (auto traj_name : traj_names) {
    BOOST_LOG_TRIVIAL(debug)<<"Loading file "<<traj_name;
    int64_t count_cnt;
    int64_t time_cnt;
    file.LoadTrajectoryCounts(traj_name, seirc, count_cnt);
    file.LoadTrajectoryTimes(traj_name, time, time_cnt);
    assert(count_cnt==time_cnt);
    if (time_cnt>max_trajectory_len) {
      max_trajectory_len=time_cnt;
      max_trajectory_name=traj_name;
    }

    for (int64_t i=0; i<time_cnt; ++i) {
      int64_t* entry=&seirc[i*compartment_cnt];
      one_percent_clinical.observe(entry, time[i]);
      five_percent_clinical.observe(entry, time[i]);
      density.observe(entry, time[i]);
    }
    one_percent_clinical.done_trajectory();
    five_percent_clinical.done_trajectory();
    density.done_trajectory();

    total_infected.observe(seirc, time, time_cnt);
  }
  
  file.Close();

  HDFFile out(outfilename);
  out.Open();
  //out.Save2DPDF(interpolant, x, y, "ensemble2d");
  out.Close();
  return 0;
}
