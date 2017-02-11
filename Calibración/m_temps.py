import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd
z=np.zeros(21)
#data = pd.read_csv('c_temps.csv')
#make the 2MHz bins
inde = np.linspace(50, 90, 5243)
for i in range(len(inde)):
    inde[i] = round(inde[i], 2)

for j in range(21):
    for i in range(len(inde)):
        if inde[i] == 50.00 + (2 * j):
            z[j] = i
