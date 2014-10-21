import logging
import numpy as np
import matplotlib.pyplot as plt
import h5py
from default_parser import DefaultArgumentParser

logger=logging.getLogger(__file__)

def per_pen_trajectory(traj_dset, initial_values, pen_cnt):
    cnt=len(traj_dset)
    logger.debug("{0} values in dataset".format(cnt))
    # The dataset is an array of tuples, not a 2d array.
    pen=np.zeros((pen_cnt, 4, 4*per_pen), np.int64)
    t=np.zeros((pen_cnt, 4*per_pen), np.float64)
    obs_cnt=np.ones(pen_cnt, np.int64)
    for idx, val in enumerate(traj_dset):
        who, pen, transition, t=val
        pen_idx=obs_cnt[pen]
        pen[pen,]

    # transition 0 takes S to E
    # transition 1 takes E to I
    # transition 2 takes I to R

def foreach_trajectory(f, func):
    trajectories=f['/trajectory']
    for trajectory_name in trajectories:
        dset=trajectories[trajectory_name]
        func(dset)

def seir_initial(f):
    initial_values=f["/trajectory"].attrs["Initial Values"]
    if len(initial_values)<4
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

