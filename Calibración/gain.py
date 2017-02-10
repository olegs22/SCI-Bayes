import numpy as np
from test import bandwidth, calibration
from glob import glob
import matplotlib.pyplot as plt

#The frequency band
freqs, vals, BW = bandwidth(50.0, 90.0)


#firts block is for the hour 23
data1 = glob('/Users/Oleg/desktop/Data/2013-06-13-23/2013-06-13-23-*-*_Ch_2_antenna.dat')
data2 = glob('/Users/Oleg/desktop/Data/2013-06-13-23/2013-06-13-23-*-*_Ch_2_noise.dat')
data3 = glob('/Users/Oleg/desktop/Data/2013-06-13-23/2013-06-13-23-*-*_Ch_2_short.dat')
data4 = glob('/Users/Oleg/desktop/Data/2013-06-13-23/2013-06-13-23-*-*_Ch_2_50ohm.dat')

num1 = np.size(data1)
num2 = np.size(data2)
x = num1 / num2
index = np.arange(1, num2 + 1)
T_amb = np.zeros(num2)
T_meas = np.zeros(num2)
y = num1 + num2 + np.size(data3) + np.size(data4)
print y  , x
"""
d = {}
for i in range(np.size(index)):
    d['asec_%s'% str(i)] = data1[i*x : x*index[i]]


val={}
for i in range(np.size(index)):
    val['antenna_%s'% str(i)] = [np.loadtxt(j) for j in d['asec_%s'% str(i)]]

noise = [np.loadtxt(j) for j in data2]
short = [np.loadtxt(k) for k in data3]
ohm = [np.loadtxt(l) for l in data4]
for i in range(np.size(index)):
    for k in range(x):
        val['antenna_%s'% str(i)][k] = val['antenna_%s'% str(i)][k][vals[0]:vals[1]] - 93.0

    noise[i] = noise[i][vals[0]:vals[1]]
    short[i] = short[i][vals[0]:vals[1]]
    T_amb[i] = ohm[i][0] + 273.15
    ohm[i] = ohm[i][vals[0]:vals[1]]



for i in range(np.size(index)):
    val['antenna_%s'% str(i)] = calibration(val['antenna_%s'% str(i)],
                                            noise[i][:], short[i][:], ohm[i][:], T_amb[i], freqs, 1.0)

cal_list = np.concatenate(val.values(), axis=0)
T_cal = np.mean(cal_list, axis = 0)
plt.plot(freqs, T_cal)
plt.show()

#np.savetxt('/Users/Oleg/documents/21cm_analysis/calibration/2013-06-13-23.dat', T_cal)

#dict.get() gets me the valeus in said key
#dict.values() prints the values of the whole dictionary

"""
