import logging
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import h5py

logger=logging.getLogger(__file__)

def smoothed2d():
    imf=h5py.File("image.h5")
    x=imf['images/ensemble2dx']
    y=imf['images/ensemble2dy']
    ds=imf['images/ensemble2d']

    extents=[np.min(x), np.max(x), np.min(y), np.max(y)]
    logger.debug(extents)
    logger.debug("dataset shape {0}".format(ds.shape))
    logger.debug("x shape {0}".format(x.shape))
    logger.debug("y shape {0}".format(y.shape))

    fig, ax=plt.subplots(1, 1, figsize=(6,6))
    print(type(fig))
    #ax=fig.add_axes(extents)
    axes=ax.imshow(ds, origin="lower",
        extent=extents, interpolation="none",
        aspect=0.6*(extents[1]-extents[0])/(extents[3]-extents[2]),
        cmap=plt.cm.YlGn)
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Individuals [count]")
    ax.set_title("Prevalence of Infection (E+I)")
    logger.debug(type(axes))
    plt.savefig("ensemble_smoothed.pdf", format="pdf")


logging.basicConfig(level=logging.DEBUG)
smoothed2d()
