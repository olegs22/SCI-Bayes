from __future__ import division
import emcee
import corner
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
from model_fit import *

model = 'tanh'

#freqs, temp = np.loadtxt('scihi_fore.txt', unpack = True)
freqs, temp, err = np.loadtxt('dataset7.txt', unpack = True)
#err = radiometer(temp, 0.83, np.diff(freqs)[0])
#prior list for the gaussian model
#                           T         nu     sigma     c0       c1      c2     c3
priors_gauss = np.array([0.,1000., 40.,100., 0.,40., 0.,10., -10.,0., -1.,1.])
#stuff needed for the mcmc
nwalkers, b_steps, steps = 200, 100, 1000

#prior list for the tanh model

priors_tanh = np.array([0.,1., 0.,20., 0.,5, 1000.,10000., 0.,20., 0.,5., 0.,25., 0,30., 0.,5.,
                        0., 10, -10.,0., -1.,1.])

#finding the midpoint of each prior range
ndim = int(np.size(priors_tanh)/2)
mid = np.zeros(ndim)
for i in range(ndim):
    mid[i] = priors_tanh[2*i] + priors_tanh[2*i + 1]
mid = 0.5 * mid

#maximazing the likelihood
fun = lambda *args: -lnhood(*args)
result = op.minimize(fun, mid, args = (model, temp, freqs, err), method = 'L-BFGS-B')
mu = result.x[9:]

#disperse the walkers around the midpoint for each parameter
z = np.zeros((ndim, nwalkers))
h1 = 1e-1
pos1=[]
for k in range(ndim):
    z[k,:] = mid[k] + h1 * np.random.rand(nwalkers)

for i in range(nwalkers):
    pos1.append(np.array([z[0,i],z[1,i],z[2,i],z[3,i],z[4,i],z[5,i],z[6,i],z[7,i],z[8,i],z[9,i],z[10,i],z[11,i]]))


#for i in range(nwalkers):
#    pos1.append(np.array([z[0,i],z[1,i],z[2,i],z[3,i],z[4,i],z[5,i],z[6,i],
#                        z[7,i],z[8,i],z[9,i],z[10,i],z[11,i],z[12,i],
#                        z[13,i],z[14,i],z[15,i],z[16,i]]))



#burn_in phase
print('this is just a test')
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(model, temp, freqs, err, priors_tanh, mu))
pos, prob, state = sampler.run_mcmc(pos1, b_steps)
print('the test has ended')
"""
sampler.reset()

sampler.run_mcmc(pos, steps,rstate0=state)
print sampler.acceptance_fraction.mean()

#save the chain.
sample = sampler.chai.reshape((-1,ndim))
np.savetxt('test_normal.txt', sample)

fig = corner.corner(sample)
fig.savefig('fit_test.pdf')
"""
