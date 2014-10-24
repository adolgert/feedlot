import logging
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import h5py
import quickpen

logger=logging.getLogger(__file__)


def plot_total_infected(data):
    plt.plot(times, seir[:,compartment_idx],
                color=colors[compartment_idx])
    plt.title("Prevalence for a Single Instance")
    labels=("susceptible", "latent", "infectious", "removed")
    legend=plt.legend(labels)
    plt.xlabel("Time since initial [days]")
    plt.ylabel("Individual count")
    plt.savefig("seir_single.pdf", format='pdf')


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
    plt.savefig("seir_single.pdf", format='pdf')


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
    plt.savefig("multiples.pdf", format='pdf')


def end_time_plot(end_times):
    end_times.sort()
    cnt=len(end_times)

    n, bins, patches=plt.hist(end_times, int(cnt/5), normed=1,
        facecolor='orange', alpha=0.75)
    plt.xlabel("End Time [days]")
    plt.ylabel("Probability")
    plt.title("Time of Last Recovery")
    plt.grid(True)
    plt.savefig("end_time_hist.pdf", format="pdf")

def total_infected_plot(total_infected):
    total_infected.sort()
    cnt=len(total_infected)
    n, bins, patches=plt.hist(total_infected, 40, normed=1,
        facecolor='orange', alpha=0.75)
    plt.xlabel("Total Infected")
    plt.ylabel("Probability")
    plt.title("Total Infected Over Course")
    plt.grid(True)
    plt.savefig("total_infected_hist.pdf", format="pdf")

def trajectory_density_plot(m1, m2, name):
    xmin=min(m1)
    xmax=max(m1)
    ymin=min(m2)
    ymax=max(m2)
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = scipy.stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(Z, cmap=plt.cm.gist_earth_r,
              extent=[xmin, xmax, ymin, ymax])
    # ax.plot(m1, m2, 'k.', markersize=1)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    plt.savefig("trajectory_density_{0}.pdf".format(name), format="pdf")




def test_one():
    f=h5py.File("rider.h5","r")
    tr,times,obs=quickpen.per_pen_trajectory(f, "dset1-1")
    small_multiples(tr, times, obs, 4)

if __name__=='__main__':
    test_one()
