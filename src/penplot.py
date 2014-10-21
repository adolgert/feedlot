import logging
import numpy as np
import matplotlib.pyplot as plt
import h5py
import quickpen

logger=logging.getLogger(__file__)


def single_trajectory(seir, times):
    colors=["blue", "green", "brown", "black"]
    for compartment_idx in range(4):
        plt.plot(times, seir[:,compartment_idx],
                color=colors[compartment_idx])
    plt.title("Prevalence for a Single Instance")
    labels=("susceptible", "latent", "infectious", "removed")
    legend=plt.legend(labels)
    plt.xlabel("Time since initial [days]")
    plt.ylabel("Individual count")
    plt.savefig("out.pdf", format='pdf')


def small_multiples(seir_pen, times_pen, obs_pen, column_cnt):
    pen_cnt=seir_pen.shape[0]
    row_cnt=pen_cnt//column_cnt;
    max_time=np.max(times_pen)
    max_cnt=np.max(seir_pen)

    font={"weight" : "medium",
        "size" : 5}
    plt.rc("font", **font)

    fig=plt.figure(1)
    fig, axes=plt.subplots(nrows=row_cnt, ncols=column_cnt,
        figsize=(6,6))
    plt.subplots_adjust(hspace=0.5)
    colors=["blue", "green", "purple", "black"]
    for pen_idx in range(pen_cnt):
        axes=plt.subplot(row_cnt, column_cnt, pen_idx+1)
        axes.set_xlim([0, max_time])
        axes.set_ylim([0, max_cnt])
        plt.title("pen {0}".format(pen_idx))
        mn=0
        mx=obs_pen[pen_idx]+1
        for compartment_idx in range(1,3):
            plt.plot(times_pen[pen_idx,mn:mx],
                    seir_pen[pen_idx,compartment_idx,mn:mx],
                    color=colors[compartment_idx])
    fig.suptitle("Exposed and Infected Across All Pens", fontsize=14)
    plt.savefig("out.pdf", format='pdf')


def test_one():
    f=h5py.File("rider.h5","r")
    tr,times,obs=quickpen.per_pen_trajectory(f, "dset1-1")
    small_multiples(tr, times, obs, 4)

if __name__=='__main__':
    test_one()
