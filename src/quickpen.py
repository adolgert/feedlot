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

