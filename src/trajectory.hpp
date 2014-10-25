#ifndef _TRAJECTORY_HPP_
#define _TRAJECTORY_HPP_ 1
#include <vector>
#include <cstdint>
#include <cstddef>

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
  virtual void Step(TrajectoryEntry sirt)=0;
  virtual const std::vector<TrajectoryEntry>& Trajectory() =0;
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


/*! Save the whole trajectory.
 */
class TrajectorySave : public TrajectoryObserver
{
  std::vector<TrajectoryEntry> trajectory_;
  int64_t cnt_;
 public:
  TrajectorySave(int64_t cnt);
  virtual ~TrajectorySave();
  virtual void Step(TrajectoryEntry seirt) override;
  virtual const std::vector<TrajectoryEntry>& Trajectory();
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
  PercentTrajectorySave();
  virtual ~PercentTrajectorySave();

  virtual void Step(TrajectoryEntry seirt) override;
  virtual const std::vector<TrajectoryEntry>& Trajectory();
};


/*! Save the whole trajectory.
 */
class PenTrajectorySave : public PenTrajectoryObserver
{
  std::vector<PenTrajectory> trajectory_;
  std::vector<TrajectoryEntry> initial_;
  int64_t cnt_{0};
 public:
  PenTrajectorySave(size_t ind_cnt);
  virtual ~PenTrajectorySave();
  virtual void SetInitial(const std::vector<TrajectoryEntry>& init);
  virtual void Step(PenTrajectory entry) override;
  const std::vector<PenTrajectory>& Trajectory();
  const std::vector<TrajectoryEntry>& PenInitial();
};


#endif
