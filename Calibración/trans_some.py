import numpy as np
from glob import glob
from power_to_temperature import *

direct = '/Users/Oleg/Documents/SCI-Bayes/antenna_beam/'
d_out = '/Users/Oleg/Documents/gsm/antenna_beam/'
freq = np.arange(50,91)
for i in range(len(freq)):
    data_p=  glob(direct+'0**MHz.txt')

data_b = [np.loadtxt(p) for p in data_p]

print np.shape(data_b)
for i in range(41):
    data_b[i][:,2] = Radio_source_trans(data_b[i][:,2], freq[i], 1e6)
    np.savetxt(d_out+'freq_%2s.txt'% str(i), data_b[i])
