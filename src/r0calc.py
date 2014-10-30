"""
This estimates the infection hazard necessary to find a certain R0.
It integrates over the survival of the infectious state.
"""

import numpy as np
import scipy
import scipy.integrate
import scipy.stats
import scipy.optimize


class Integrand(object):
    def __init__(self, beta):
        self.rv=scipy.stats.gamma(a=3.969, scale=1/1.107)
        self.beta=beta

    def __call__(self, t):
        return self.rv.sf(t)*t*self.beta
    
def determine_accuracy():
    integrand=Integrand(0.077)
    for i in [10, 20, 1000]:
        value, uncertainty=scipy.integrate.quad(integrand, 0, i)
        print("{0} value {1}, uncertainty {2}".format(i, value, uncertainty))

class BetaEstimate(object):
    def __init__(self, goal):
        self.r0=goal
    def __call__(self, beta):
        integrand=Integrand(beta)
        value, uncertainty=scipy.integrate.quad(integrand, 0, 20)
        return (value-self.r0)**2

for R0 in [2.0, 4.0, 8.0, 112.0]:
    R02=BetaEstimate(R0)
    res=scipy.optimize.minimize(R02, 0.12*R0)
    print("Minimum for R0={0} is beta={1}".format(R0, res.x[0]))
    print(res)

determine_accuracy()
