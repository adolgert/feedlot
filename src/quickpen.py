import logging
import numpy as np
import matplotlib.pyplot as plt
import h5py
from default_parser import DefaultArgumentParser

logger=logging.getLogger(__file__)

def trajectories(h5f):
    assert(issubclass(type(h5f), h5py._hl.files.File))
    return list(h5f["/trajectory"].keys())

def initial_values(h5f, traj_dset_name):
    assert(issubclass(type(h5f), h5py._hl.files.File))
    initial_ds=h5f["/trajectory/{0}/initialpencount".format(traj_dset_name)]
    logger.debug("Initial values dataset size {0}".format(initial_ds.shape))
    return initial_ds

def per_pen_trajectory(h5f, traj_dset_name):
    assert(issubclass(type(h5f), h5py._hl.files.File))
    initial=initial_values(h5f, traj_dset_name)
    total_individuals=np.sum(initial)
    pen_cnt=initial.shape[0]
    per_pen=total_individuals//pen_cnt;
    traj_dset=h5f["/trajectory/{0}/trajectory".format(traj_dset_name)]
    # The dataset is an array of tuples, not a 2d array.
    trajectory=np.zeros((pen_cnt, 4, 3*per_pen+1), np.int64)
    trajectory[:,:,0]=initial
    times=np.zeros((pen_cnt, 3*per_pen+1), np.float64)
    last_obs=np.zeros(pen_cnt, np.int64)
    for idx, val in enumerate(traj_dset):
        who, pen, transition, t=val
        pen_step_idx=last_obs[pen]
        s,e,i,r=trajectory[pen,:,pen_step_idx]
        logger.debug((pen, pen_step_idx, transition, (s,e,i,r)))
        # transition 0 takes S to E
        # transition 1 takes E to I
        # transition 2 takes I to R
        if transition==0:
            s-=1
            e+=1
        elif transition==1:
            e-=1
            i+=1
        elif transition==2:
            i-=1
            r+=1
        else:
            print("unknown transition")
        pen_step_idx+=1
        if (pen_step_idx<trajectory.shape[2]):
            trajectory[pen,:,pen_step_idx]=np.array((s,e,i,r),np.int64)
            times[pen,pen_step_idx]=t
            last_obs[pen]=pen_step_idx
        else:
            logger.error("Too many entries for pen {0}".format(pen))
            last_obs[pen]=pen_step_idx
    return trajectory, times, last_obs


def total_trajectory(h5f, traj_dset_name):
    initial=initial_values(h5f, traj_dset_name)
    traj_dset=h5f["/trajectory/{0}/trajectory".format(traj_dset_name)]
    seir_start=np.sum(initial, axis=0)
    seir=np.zeros((len(traj_dset)+1,4), np.int)
    times=np.zeros(len(traj_dset)+1, np.float64)
    seir[0,:]=seir_start
    times[0]=0.0
    seir_now=seir_start
    for idx, val in enumerate(traj_dset):
        who, pen, transition, t=val
        # transition 0 takes S to E
        # transition 1 takes E to I
        # transition 2 takes I to R
        if transition==0:
            seir_now[0]-=1
            seir_now[1]+=1
        elif transition==1:
            seir_now[1]-=1
            seir_now[2]+=1
        elif transition==2:
            seir_now[2]-=1
            seir_now[3]+=1
        else:
            print("unknown transition")
        times[idx+1]=t
        seir[idx+1,:]=seir_now
    return seir, times


def add_to_binned_trajectory(binned, seir, times):
    time_idx=0
    cur_val=seir[time_idx,:]
    for bin_idx in range(binned.shape[0]):
        while time_idx<len(times) and times[time_idx]<bin_idx:
            cur_val=seir[time_idx,:]
            time_idx+=1
        binned[bin_idx,:]+=cur_val
    return binned


class total_infected(object):
    def __init__(self, cnt):
        self.data=np.zeros(cnt, np.int)
        self.idx=0

    def observe(self, total_trajectory, times):
        self.data[self.idx]=total_trajectory[-1,3]
        self.idx+=1


class time_to_recovery(object):
    def __init__(self, cnt):
        self.data=np.zeros(cnt, np.float64)
        self.idx=0

    def observe(self, total_trajectory, times):
        self.data[self.idx]=times[-1]
        self.idx+=1

class all_states(object):
    def __init__(self, cnt):
        self.seir_=list()
        self.t=list()

    def observe(self, total_trajectory, times):
        self.seir_.append(total_trajectory)
        self.t.append(times)

    def seir(self):
        return np.vstack(self.seir_)

    def times(self):
        return np.hstack(self.t)

def summary_of_ensemble(f, max_cnt):
    cnt=len(trajectories(f))
    if max_cnt>cnt or max_cnt<1:
        max_cnt=cnt
    summaries=[total_infected(max_cnt), time_to_recovery(max_cnt),
        all_states(max_cnt)]
    for trajectory_name in trajectories(f)[0:max_cnt]:
        logger.debug(trajectory_name)
        total, times=total_trajectory(f, trajectory_name)
        for s in summaries:
            s.observe(total, times)
    return summaries

class FileTrajectories(object):
    def __init__(self, h5f):
        self.h5f=h5f
        self.trajectory_names=None

    def __len__(self):
        if self.trajectory_names is None:
            self.trajectory_names=trajectories(self.h5f)
        return len(self.trajectory_names)

    def __getitem__(self, intkey):
        if self.trajectory_names is None:
            self.trajectory_names=trajectories(self.h5f)
        total, times=total_trajectory(self.h5f, self.trajectory_names[intkey])
        return total, times


def foreach_trajectory(f, func):
    trajectories=f['/trajectory']
    for trajectory_name in trajectories:
        dset=trajectories[trajectory_name]
        func(dset)


def seir_initial(f):
    initial_values=f["/trajectory"].attrs["Initial Values"]
    if len(initial_values)<4:
        seir_values=np.zeros(4)
        for i in range(len(initial_values)):
            seir_values[i]=initial_values[i]
        initial_values=seir_values
    return initial_values

def showds(ds):
    print(ds)
    attributes=list()
    for x in ds.attrs:
        attributes.append("{0} {1}".format(x, ds.attrs[x][0]))
    print("  {0}".format(", ".join(attributes)))


def showproginfo(f):
    attrs=f['/trajectory'].attrs
    for x in attrs:
        if (attrs[x].dtype.type is np.bytes_):
            print("{0}: ".format(x))
            for v in attrs[x]:
                print(v.decode())
        else:
            print("{0}: {1}".format(x, attrs[x]))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    parser=DefaultArgumentParser(description="Quick look at an H5 file")
    parser.add_function("info", "Find what program made the file.")
    parser.add_function("trajectory", "Plot the trajectory")
    parser.add_function("dir", "List datasets")
    parser.add_argument("--file", dest="file", action="store",
        default="rider.h5", help="data file to read")

    args=parser.parse_args()

    filename=args.file
    f=h5py.File(filename, "r")

    if args.info:
        showproginfo(f)
    if args.trajectory:
        foreach_trajectory(f, plot_single)
    if args.dir:
        foreach_trajectory(f, showds)
    #foreach_trajectory("sirexp.h5", plot_single)
    if not parser.any_function():
        parser.print_help()

