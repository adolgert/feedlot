import logging
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import h5py
import quickpen
import lifelines

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
    colors=["blue", "green", "brown", "black", "red"]
    for compartment_idx in range(5):
        plt.plot(times, seir[:,compartment_idx],
                color=colors[compartment_idx])
    plt.title("Prevalence for a Single Instance")
    labels=("susceptible", "latent", "infectious", "removed", "clinical")
    legend=plt.legend(labels)
    plt.xlabel("Time since initial [days]")
    plt.ylabel("Individual count")
    plt.savefig("seir_single.pdf", format='pdf')


def small_multiples(seir_pen, times_pen, obs_pen, column_cnt):
    pen_cnt=seir_pen.shape[0]
    column_cnt=min(pen_cnt, column_cnt)
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
    plt.clf()
    end_times.sort()
    kmf=lifelines.KaplanMeierFitter()
    kmf.fit(end_times, label="Last Recovery")
    kmf.plot()
    plt.xlabel("End Time [days]")
    plt.ylabel("Probability of Survival")
    plt.title("Time of Last Recovery")
    plt.grid(True)
    plt.savefig("end_time_hist.pdf", format="pdf")

def total_infected_plot(total_infected):
    plt.clf()
    total_infected.sort()
    kmf=lifelines.KaplanMeierFitter()
    kmf.fit(total_infected, label="Total Infected")
    kmf.plot()
    plt.xlabel("Total Infected [count]")
    plt.ylabel("Probability of at Least This Many")
    plt.title("Total Infected By End of Outbreak")
    plt.grid(True)
    plt.savefig("total_infected_hist.pdf", format="pdf")

def survival_susceptible(when_infected):
    plt.clf()
    kmf=lifelines.KaplanMeierFitter()
    kmf.fit(when_infected, label="Total Infected")
    kmf.plot()
    plt.xlabel("Time of Infection [days]")
    plt.ylabel("Probability of Survival to This Time")
    plt.title("Survival of Infected During Outbreak")
    plt.grid(True)
    plt.savefig("survival_infected_hist.pdf", format="pdf")

def total_infected_count_plot(total_infected):
    plt.clf()
    maxval=np.max(total_infected)
    data=np.zeros(maxval+1, dtype=np.double)
    x=np.arange(0, maxval+1)
    for icnt in total_infected:
        data[icnt]+=1
    logger.debug("total_points data {0} sum {1}".format(np.sum(data),
        np.sum(total_infected)))
    data/=np.sum(total_infected)
    plt.plot(x, data, 'go')
    plt.title("Distribution of Total Infected")
    plt.xlabel("Individuals [count]")
    plt.ylabel("Probability")
    plt.savefig("total_infected_points.pdf", format="pdf")

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

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(Z, cmap=plt.cm.YlGn, origin="lower",
              extent=[xmin, xmax, ymin, ymax],
              aspect=0.6*(xmax-xmin)/(ymax-ymin))
    # ax.plot(m1, m2, 'k.', markersize=1)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Individuals [count]")
    ax.set_title("Kernel Density Estimate of {0}".format(name))
    plt.savefig("trajectory_density_{0}.pdf".format(name), format="pdf")

def trajectory_image(hist, x, y, name):
    if type(hist) is h5py.Dataset:
        hist=hist[:,:] # Take it out of hdf5
        x=x[:]
        y=y[:]
    # Looking for 95%, 50%, and 5%
    reverse_percentile=np.cumsum(hist[::-1,:], axis=0)
    logger.debug(reverse_percentile[-1,0:200])
    # max_accum=np.repeat(np.expand_dims(reverse_percentile[-1,:], axis=0),
    #     y.shape[0], axis=0)
    max_accum=np.max(reverse_percentile[-1,:])
    logger.debug("ti: x n {0} x {1}, y n {2} x {3}".format(min(x),
        max(x), min(y), max(y)))
    logger.debug("trajectory_image max_accum {0} {1} {2}".format(3,
        reverse_percentile.shape, y.shape))
    ptiles=list()
    for ptile in [0.05, 0.5, 0.95]:
        diff=np.abs(reverse_percentile-ptile*max_accum)
        rev95=np.argmin(diff, axis=0)
        logger.debug("trajectory_image shape {0} diff {1} rev95 {2}".format(
            rev95.shape, diff.shape, rev95.dtype))
        ptiles.append(y[-rev95])

    project_onto_y=np.sum(hist, axis=1)
    accum=np.cumsum(project_onto_y[::-1])
    cutoff=0.0001*accum[-1]
    excise=len(accum[accum<cutoff])
    hist=hist[:-excise,:]
    y=y[:-excise]
    logger.debug("ti: x n {0} x {1}, y n {2} x {3}".format(min(x),
        max(x), min(y), max(y)))
    xmin=min(x)
    xmax=max(x)
    ymin=min(y)
    ymax=max(y)
    logger.debug("hist from min {0} to max {1}.".format(np.min(hist),
            np.max(hist)))
    hist=np.log(hist)
    mn=np.max(hist)-10
    hist[hist<mn]=-np.Inf
    #hist=hist/np.max(hist)
    plt.clf()
    fig=plt.figure()
    plt.plot(x, ptiles[0], color="blue")
    plt.plot(x, ptiles[1], color="black")
    plt.plot(x, ptiles[2], color="blue")
    ax=fig.add_subplot(111)
    ax.imshow(hist, cmap=plt.cm.YlGn, origin="lower",
        extent=[xmin, xmax, ymin, ymax],
        aspect=0.6*(xmax-xmin)/(ymax-ymin))
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Infected Individuals [count]")
    ax.set_title("Infecteds")
    plt.savefig("trajectory_lines_{0}.pdf".format(name), format="pdf")

def plot_trajectory_lines(file_trajectories):
    plt.clf()
    cnt=min(len(file_trajectories), 1000)
    for i in range(cnt):
        total, times=file_trajectories[i]
        infected=total[:,1]+total[:,2] # e + i
        plt.plot(times, infected, alpha=20/cnt, color="green")

    plt.xlabel("Time [days]")
    plt.ylabel("Individuals [count]")
    plt.title("Exposed and Infectious for Multiple Trajectories")
    plt.savefig("trajectory_lines.pdf", format="pdf")
    plt.clf()

def prevalence_by_day(binned):
    plt.clf()
    logger.debug("nonzeros {0}".format(np.nonzero(binned[:,1:3])))
    daycnt=np.nonzero(binned)[0][-1]+1
    logger.debug("last nonzero day {0}".format(daycnt))
    xl=np.linspace(0, daycnt-1, daycnt)
    xr=np.linspace(1, daycnt, daycnt)
    colors=["blue", "green", "black", "red"]
    for compartment_idx in [1, 2]:
        plt.hlines(binned[:daycnt,compartment_idx], xl, xr,
            color=colors[compartment_idx-1])
    plt.xlabel("Time [days]")
    plt.ylabel("Prevalence [individuals]")
    plt.title("Average Prevalence Over All Realizations")
    plt.savefig("prevalencebyday.pdf", format="pdf")

def test_one():
    f=h5py.File("rider.h5","r")
    tr,times,obs=quickpen.per_pen_trajectory(f, "dset1-1")
    small_multiples(tr, times, obs, 4)

if __name__=='__main__':
    test_one()
