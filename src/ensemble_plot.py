'''
Trying to find a way to plot contours that represent the density of
realizations of a stochastic system.

Assume that we generate some 1000 sets of y values over a series
of x locations and store them in a matrix runs[y instance,x location].

We can accumulate how many of those values fall in a particular 
square between (y0, y0+dy) and (x0, x0+dx). That gives a raw 
probability, but the resulting plot graphically darkens the 
parts of the plot where nothing is happening. It doesn't seem
to show the cumulative percentiles well. I think the reason is
that the total probability is smeared out in regions of interest.


'''
import logging
import numpy as np
import custom_colormap
from pylab import *


logger=logging.getLogger('ensemble_plot')


def plot_it():
    res, xaxis=sample_sim(250, 1000)
    logger.debug('The data array has shape %s.' % str(res.shape))

    percentiles=pull_percentiles(res)

    cont, x, y=runs_to_contour(res, xaxis, (1,1))
    logger.debug('The array to contour has shape %s.' % str(cont.shape))
    logger.debug('And there are new x and y axis arrays to go with it'+\
            ' of size %s and %s.' % (str(x.shape), str(y.shape)))
    dx=0.1
    levels=np.arange(0,1+dx, dx)
    levels=levels*np.max(cont)
    levels=[0.5, 1, 1.5, 2, 2.5, 3]

    figure()
    helix_map=custom_colormap.light_gray_cmap(0.7, gamma=.8)
    CS=contourf(x, y, cont, levels, cmap=helix_map, origin='lower')
    colorbar(CS)
    plot(xaxis, percentiles[:,0], color='black')
    plot(xaxis, percentiles[:,1], color='black', linewidth=2)
    plot(xaxis, percentiles[:,2], color='black')
    show()


def runs_to_contour(runs, x_axis, subset=(1,1), dtype=np.int):
    # runs.shape is (1000,250), (runs, y values)
    # map multiple y(x) plots into a grid.
    # x values are the left-hand side of each box.
    x_lim=(np.min(x_axis), np.max(x_axis))
    x=np.arange(np.min(x_axis),np.max(x_axis)+1,subset[0])
    y=np.arange(np.min(runs),  np.max(runs)+1,  subset[1])

    x_target=((x_axis-x_lim[0])/subset[0]).astype(dtype)
    data_min=np.min(runs)

    target=np.zeros((len(y), len(x)), dtype=np.float)
    # Accumulate counts into a 2D matrix.
    for i in range(runs.shape[0]):
        y_target=((runs[i,:]-data_min)/subset[1]).astype(dtype)
        target[y_target, x_target]+=1

    target=np.log10(1+target)
    # Don't rescale target because log10 counts should be in legend.
    # Normalize across x values.
    #for i in range(len(x)):
    #    target[:,i]/=np.max(target[:,i])

    return target, x, y



def pull_percentiles(runs, levels=[0.05, 0.5, 0.95]):
    '''
    Given an array of 1000 runs of 250 steps in 
    an array runs[1000,250], at what value of y
    are levels percent of the runs less than that y?
    This returns an array that's percentiles[250,3],
    for the given example.
    '''
    levels=(np.array(levels)*runs.shape[0]).astype(np.int)
    logger.debug('percentile levels have index %s' % str(levels))
    percentiles=np.zeros((runs.shape[1],len(levels)))

    # Sort y values for each timestep and just count
    # until you reach 5%, 50%, or 95% of them.
    runs_sort=np.sort(runs, 0)
    for i in range(len(levels)):
        percentiles[:,i]=runs_sort[levels[i], :]
    return percentiles


def other_percentile_try():
    cdf=np.cumsum(target,0)
    logger.debug('cdf.shape %s' % str(cdf.shape))
    logger.debug(np.max(cdf,0).shape)
    cdfmax=np.max(cdf,0)
    logger.debug(cdfmax)
    logger.debug(target.shape)
    lower=np.zeros(target.shape[1])
    upper=np.zeros(target.shape[1])
    for i in range(target.shape[1]):
        col=cdf[:,i]/float(cdfmax[i])
        logger.debug(col)
        ll=np.where(col>percent)[0]
        print(ll)
        if len(ll)>0:
            lower[i]=y[ll[-1]]
        else:
            logger.debug('lower is the max')
            y[-1]
        ul=np.where(col<(1-percent))[0]
        print(ul)
        if len(ul)>0:
            upper[i]=y[ul[0]]
        else:
            upper[i]=y[0]
            logger.debug('upper is zero')
    logger.debug('upper and lower')
    logger.debug(upper)
    logger.debug(lower)
    return lower, upper


def sample_sim(day_cnt, sim_cnt):
    def logistic(t,t0,s=5):
        return 1/(1+np.exp(-(t-t0)/s))

    res=np.zeros((sim_cnt,day_cnt),np.int)
    base_prob=lambda d: 0.04*logistic(d,100,7)
    #base_prob=lambda d: 0.05
    
    sim_days=np.arange(day_cnt,dtype=np.float)
    prob=base_prob(sim_days)
    for sim_idx in range(sim_cnt):
        population=963
        for d in sim_days:
            infected=int(np.random.binomial(population,prob[d]))
            res[sim_idx,d]=infected
            population -= infected
            if population<=0: break;
    return(res, sim_days)



if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    plot_it()

