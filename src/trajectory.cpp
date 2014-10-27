#include <cmath>
#include <cstdlib>
#include "trajectory.hpp"


TrajectorySave::TrajectorySave(int64_t cnt) : trajectory_(cnt), cnt_{0} {}
TrajectorySave::~TrajectorySave() {}
void TrajectorySave::Step(TrajectoryEntry seirt) {
  if (cnt_==trajectory_.size()) {
    trajectory_.resize(2*cnt_);
  }
  trajectory_[cnt_]=seirt;
  ++cnt_;
}
const std::vector<TrajectoryEntry>& TrajectorySave::Trajectory() {
  trajectory_.resize(cnt_);
  return trajectory_;
}

PercentTrajectorySave::PercentTrajectorySave() {}
PercentTrajectorySave::~PercentTrajectorySave() {}

void PercentTrajectorySave::Step(TrajectoryEntry seirt) {
  if (0==step_) {
    last_=seirt;
    threshhold_=std::floor(percent_*(seirt.s+seirt.e+seirt.i+seirt.r));
    trajectory_.emplace_back(seirt);
  } else {
    bool ps=std::labs(seirt.s-last_.s)>threshhold_;
    bool pe=std::labs(seirt.e-last_.e)>threshhold_;
    bool pi=std::labs(seirt.i-last_.i)>threshhold_;
    bool pr=std::labs(seirt.r-last_.r)>threshhold_;
    if (ps||pe||pi||pr) {
      trajectory_.emplace_back(seirt);
      last_=seirt;
    }
  }
  ++step_;
}

const std::vector<TrajectoryEntry>&
PercentTrajectorySave::Trajectory() {
  return trajectory_;
}

PenTrajectorySave::PenTrajectorySave(size_t ind_cnt) : trajectory_{ind_cnt*4} {}
PenTrajectorySave::~PenTrajectorySave() {}
void PenTrajectorySave::SetInitial(const std::vector<TrajectoryEntry>& init) {
  initial_=init;
}
void PenTrajectorySave::Step(PenTrajectory entry)  {
  if (cnt_==trajectory_.size()) {
    trajectory_.resize(2*cnt_);
  }
  trajectory_[cnt_]=entry;
  ++cnt_;
}
const std::vector<PenTrajectory>& PenTrajectorySave::Trajectory() {
  trajectory_.resize(cnt_);
  return trajectory_;
}
const std::vector<TrajectoryEntry>& PenTrajectorySave::PenInitial() {
  return initial_;
}