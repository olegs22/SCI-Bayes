"""

synthetic_data.py

Author: Jordan Mirocha
Affiliation: UCLA
Created on: Sat Apr 16 14:00:42 PDT 2016

Description:

"""

import numpy as np
from models import Global21cm, Galaxy, radiometer_eq
import matplotlib.pyplot as plt

freqs = np.arange(40, 161)

signal = Global21cm(freqs)
galaxy = Galaxy(freqs)


seed = 36194

galaxy = Galaxy(freqs, 'logpoly')
fg = galaxy(7.5, -2.55, 0.1)
g21 = signal(-100., 80., 10.) / 1e3
error = radiometer_eq(fg+g21, tint=0.83, channel=np.diff(freqs)[0])

np.random.seed(seed)
noise = np.random.normal(loc=0., scale=error, size=freqs.size)

#np.random.seed(seed)
data = fg + g21 + noise
direct = '/Users/Oleg/Documents/SCI-Bayes/'
np.savetxt(direct+'mock3.txt', np.array([freqs, data, error]).T)

"""
Dataset 1.


seed = 24691

fg = galaxy(5e2, -2.5)
g21 = signal(-100, 80., 10.) / 1e3
error = radiometer_eq(fg+g21, tint=1.0, channel=np.diff(freqs)[0])

np.random.seed(seed)
noise = np.random.normal(loc=0., scale=error, size=freqs.size)

np.random.seed(seed)
data = fg + g21 + noise

np.savetxt('dataset1.txt', np.array([freqs, data, error]).T)



#Dataset 2.


seed = 36194

galaxy = Galaxy(freqs, 'logpoly')
fg = galaxy(7.5, -2.55, 0.1, 0.4)
g21 = signal(-100., 80., 10.) / 1e3
error = radiometer_eq(fg+g21, tint=10.0, channel=np.diff(freqs)[0])

np.random.seed(seed)
noise = np.random.normal(loc=0., scale=error, size=freqs.size)

#np.random.seed(seed)
data = fg + g21 + noise
direct = '/Users/Oleg/Documents/SCI-Bayes/'
np.savetxt(direct+'mock2.txt', np.array([freqs, data, error]).T)





seed = 190345

signal = Global21cm(freqs, 'tanh')
galaxy = Galaxy(freqs, 'logpoly')
fg = galaxy(7.5, -2.55, 0.1)
g21 = signal(0.8, 12., 4., 2000., 6., 2., 1.0, 4., 2.4) / 1e3
error = radiometer_eq(fg+g21, tint=10.0, channel=np.diff(freqs)[0])

np.random.seed(seed)
noise = np.random.normal(loc=0., scale=error, size=freqs.size)

#np.random.seed(seed)
data = fg + g21 + noise

np.savetxt('datas1.txt', np.array([freqs, data, error]).T)

#Dataset 4.


seed = 10187

signal = Global21cm(freqs, 'tanh')
galaxy = Galaxy(freqs, 'logpoly')

fg1 = galaxy(7.5, -2.55, 0.1, 0.05)
fg2 = galaxy(7.2, -2.45, 0.2, 0.01)

g21 = signal(18., 24.5, 3.9, 800, 9.7, 4.5, 0.25, 10., 1.) / 1e3

error1 = radiometer_eq(fg1+g21, tint=10.0, channel=np.diff(freqs)[0])
error2 = radiometer_eq(fg2+g21, tint=10.0, channel=np.diff(freqs)[0])

np.random.seed(seed)
noise1 = np.random.normal(loc=0., scale=error1, size=freqs.size)
noise2 = np.random.normal(loc=0., scale=error2, size=freqs.size)

data1 = fg1 + g21 + noise1
data2 = fg2 + g21 + noise2

np.savetxt('dataset4.txt', np.array([freqs, data1, error1, data2, error2]).T)


seed = 190345

signal = Global21cm(freqs, 'tanh')
galaxy = Galaxy(freqs, 'logpoly')
fg = galaxy(7.5, -2.55, 0.1, 0.05)
#fg = np.loadtxt('foreground.txt')
g21 = signal(15., 18.5, 3.7, 1000., 8., 3., 1., 10., 1.5) / 1e3
error = radiometer_eq(fg+g21, tint=10.0, channel=np.diff(freqs)[0])

np.random.seed(seed)
noise = np.random.normal(loc=0., scale=error, size=freqs.size)

#np.random.seed(seed)
data = fg + g21 + noise

np.savetxt('dataset6.txt', np.array([freqs, data, error]).T)


frqs, temp, sigma = np.loadtxt('dataset1.txt', unpack = True)
frqs2,temp2, err2 = np.loadtxt('mock1.txt', unpack = True)
fg = galaxy(5e2, -2.5)
T_b = (temp - fg) * 1000.

galaxy = Galaxy(freqs, 'logpoly')
fg2 = galaxy(7.50008337, -2.55, 0.0998)
T_b2 = (temp2-fg2) * 1000.
plt.figure(1)
plt.plot(frqs, T_b)
plt.figure(2)
plt.plot(frqs2, T_b2)
plt.show()
"""
