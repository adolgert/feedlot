#define BOOST_LOG_DYN_LINK 1
#include <algorithm>
#include <vector>
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
  : found_{false}, when_(trajectory_cnt,0), idx_{0} {
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
  const std::vector<double> data() { return when_; }
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
  void done_trajectory() {}
  std::pair<std::vector<double>,std::vector<double>> xy() {
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
    return std::make_pair(x, y);
  }
  const std::vector<double>& data() { return data_; }
 private:
  std::vector<double> data_;
  int64_t total_individuals_;
  double max_time_;
  int64_t vres_;
  int64_t hres_;
  double scale_;
};


class InfectionTimes {
 public:
  InfectionTimes(int64_t total_individuals, int64_t trajectory_cnt)
  : compartment_cnt_(5), when_(trajectory_cnt*total_individuals,0), idx_{0},
    have_init_(true), last_(compartment_cnt_, 0),
    trajectory_cnt_(trajectory_cnt), trajectory_idx_(0) { }
  void observe(int64_t const* seirc, double time) {
    if (have_init_ && trajectory_idx_<trajectory_cnt_) {
      if (seirc[1]>last_[1]) {
        when_[idx_]=time;
        ++idx_;
      }
    } else {
      have_init_=true;
    }
    std::copy(seirc, seirc+compartment_cnt_, last_.begin());
  }
  void done_trajectory() {
    have_init_=false;
    ++trajectory_idx_;
  }
  const std::vector<double> data() {
    BOOST_LOG_TRIVIAL(debug)<<"InfectionTimes size "<<idx_;
    when_.resize(idx_);
    std::sort(when_.begin(), when_.end());
    return when_;
  }
 private:
  bool have_init_;
  size_t idx_;
  int64_t trajectory_idx_;
  int64_t trajectory_cnt_;
  std::vector<double> when_;
  int compartment_cnt_;
  std::vector<int64_t> last_;
};


class PrevalenceIncidence {
 public:
  PrevalenceIncidence(int day_cnt, const std::vector<int64_t> initial)
  : compartment_cnt_(5), prevalence_day_start_(compartment_cnt_*(day_cnt+1), 0),
    incidence_on_day_(compartment_cnt_*(day_cnt+1), 0), initial_(initial),
    last_prevalence_(initial), last_time_(0), max_time_(day_cnt+1),
    new_since_last_prevalence_(compartment_cnt_,0)
  {
    assert(initial_.size()==compartment_cnt_);
    assert(day_cnt>0);
    assert(prevalence_day_start_.size()>0);
    assert(prevalence_day_start_.size()==(day_cnt+1)*compartment_cnt_);
  }
  void observe(int64_t const* seirc, double time) {
    while (time-last_time_>1) {
      for (int cidx=0; cidx<compartment_cnt_; ++cidx) {
        auto lp=last_prevalence_.at(cidx);
        prevalence_day_start_.at(last_time_*compartment_cnt_+cidx)+=lp;
      }
      for (int aidx=0; aidx<compartment_cnt_; ++aidx) {
        incidence_on_day_.at(last_time_*compartment_cnt_+aidx)
            +=new_since_last_prevalence_.at(aidx);
      }
      new_since_last_prevalence_.assign(compartment_cnt_, 0);
      last_time_+=1;
    }
    for (int didx=0; didx<compartment_cnt_; ++didx) {
      auto diff=seirc[didx]-last_prevalence_.at(didx);
      if (diff>0) {
        new_since_last_prevalence_.at(didx)+=diff;
      }
    }
    std::copy(seirc, seirc+compartment_cnt_, last_prevalence_.begin());
  }
  void done_trajectory() {
    for (int t=last_time_+1; t<max_time_; ++t) {
       for (int cidx=0; cidx<compartment_cnt_; ++cidx) {
        auto lp=last_prevalence_.at(cidx);
        prevalence_day_start_.at(last_time_*compartment_cnt_+cidx)=lp;
      }
    }
    last_prevalence_=initial_;
    last_time_=0;
  }
  const std::vector<int64_t>& prevalence() { return prevalence_day_start_; }
  const std::vector<int64_t>& incidence() { return incidence_on_day_; }
 private:
  int compartment_cnt_;
  // day_cnt_ x 5
  std::vector<int64_t> prevalence_day_start_;
  // day_cnt_ x 5
  std::vector<int64_t> incidence_on_day_;
  // initial states, 5
  std::vector<int64_t> initial_;
  // 5
  std::vector<int64_t> last_prevalence_;
  // 5
  std::vector<int64_t> new_since_last_prevalence_;
  int last_time_;
  int max_time_;
};


class TotalInfected {
 public:
  TotalInfected(int64_t trajectory_cnt) : count_(trajectory_cnt), idx_{0} {}
  void observe(const std::vector<int64_t>& seirc,
      const std::vector<double>& time, int64_t entry_cnt) {
    count_[idx_]=seirc[5*(entry_cnt-1)+3];
    ++idx_;
  }
  const std::vector<int64_t> data() { return count_; }
 private:
  std::vector<int64_t> count_;
  size_t idx_;
};


int main(int argc, char* argv[]) {
  std::string filename("rider.h5");
  std::string outfilename("image.h5");
  std::string log_level("info");
  int64_t hres=500;
  int64_t vres=500;
  int64_t infection_times_cnt=100;

  namespace po=boost::program_options;
  po::options_description desc("Generate ensemble plot from datasets.");
  desc.add_options()
      ("file",
      po::value<std::string>(&filename)->default_value(filename),
      "Read data from this file.")
      ("outfile",
      po::value<std::string>(&outfilename)->default_value(outfilename),
      "Write to this data file.")
      ("res",
        po::value<int64_t>(&hres)->default_value(hres),
        "Resolution of grid for aggregation")
      ("infectiontimes",
        po::value<int64_t>(&infection_times_cnt)->default_value(
          infection_times_cnt),
        "Resolution of grid for aggregation")
      ("loglevel",
      po::value<std::string>(&log_level)->default_value(log_level),
      "trace, debug, info, warn, error")
      ;

  po::variables_map vm;
  auto parsed_options=po::parse_command_line(argc, argv, desc);
  po::store(parsed_options, vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  afidd::LogInit(log_level);
  vres=hres;

  BOOST_LOG_TRIVIAL(info)<<"Reading file "<<filename;
  HDFFile file(filename);
  file.OpenRead(false);
  auto initial=file.InitialValues();
  int64_t total_individuals=
      std::accumulate(initial.begin(), initial.end(), int64_t{0});
  BOOST_LOG_TRIVIAL(debug)<<"Total animals "<<total_individuals;
  auto end=file.EndTimes();
  const auto& times=std::get<0>(end);
  const auto& events=std::get<1>(end);
  BOOST_LOG_TRIVIAL(info)<<"Found end times "<<std::get<0>(end).size();
  int sample_cnt=(std::min)(times.size(), 10ul);
  for (int et_idx=0; et_idx<sample_cnt; ++et_idx) {
    std::cout << std::get<0>(end)[et_idx] << '\t' <<
        std::get<1>(end)[et_idx] << std::endl;
  }
  int64_t trajectory_cnt=times.size();
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
  PrevalenceIncidence prevalence(
      static_cast<int>(std::ceil(max_tr_time)), initial);
  TotalInfected total_infected(trajectory_cnt);
  int64_t it_trajectory_cnt=std::min(trajectory_cnt, infection_times_cnt);
  InfectionTimes infection_times(total_individuals, it_trajectory_cnt);

  std::string max_trajectory_name;
  int64_t max_trajectory_len{0};

  auto traj_names=file.Trajectories();
  for (auto traj_name : traj_names) {
    BOOST_LOG_TRIVIAL(debug)<<"Loading dataset "<<traj_name;
    int64_t count_cnt;
    int64_t time_cnt;
    bool tc_success=file.LoadTrajectoryCounts(traj_name, seirc, count_cnt);
    if (!tc_success) {
      return -1;
    }
    bool tt_success=file.LoadTrajectoryTimes(traj_name, time, time_cnt);
    if (!tt_success) {
      return -2;
    }
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
      prevalence.observe(entry, time[i]);
      infection_times.observe(entry, time[i]);
    }
    one_percent_clinical.done_trajectory();
    five_percent_clinical.done_trajectory();
    density.done_trajectory();
    prevalence.done_trajectory();
    infection_times.done_trajectory();

    total_infected.observe(seirc, time, time_cnt);
  }
  
  file.Close();

  HDFFile out(outfilename);
  out.Open();
  out.WriteImageAttribute(max_trajectory_name, "largesttrajectory");
  out.Save1DArray(times, "endtime");
  out.Save1DArray(one_percent_clinical.data(), "onepercentclinical");
  out.Save1DArray(five_percent_clinical.data(), "fivepercentclinical");
  out.Save1DArray(total_infected.data(), "totalinfected");
  out.Save1DArray(infection_times.data(), "infectiontimes");
  out.Save2DArray(prevalence.prevalence(), 5, "prevalence");
  out.Save2DArray(prevalence.incidence(), 5, "incidence");
  auto xy=density.xy();
  out.Save2DPDF(density.data(), xy.first, xy.second, "trajectorydensity");
  out.Close();
  return 0;
}
