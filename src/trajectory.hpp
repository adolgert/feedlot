#ifndef _TRAJECTORY_HPP_
#define _TRAJECTORY_HPP_ 1


struct TrajectoryEntry {
  int64_t s;
  int64_t e;
  int64_t i;
  int64_t r;
  double t;
  TrajectoryEntry(int64_t s, int64_t e, int64_t i, int64_t r, double t)
  : s(s), e(e), i(i), r(r), t(t) {}
  TrajectoryEntry()=default;
};

class TrajectoryObserver {
public:
  virtual void Step(TrajectoryEntry sirt)=0;
  virtual const std::vector<TrajectoryEntry>& Trajectory() const =0;
};

struct PenTrajectory {
  int64_t individual;
  int64_t pen;
  int64_t transition;
  double time;
};

class PenTrajectoryObserver {
public:
  virtual void SetInitial(const std::vector<TrajectoryEntry>& init)=0;
  virtual void Step(PenTrajectory pt)=0;
};

#endif
