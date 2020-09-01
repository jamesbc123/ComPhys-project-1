# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 12:56:53 2020

@author: jamesbc
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

toi = pd.read_csv("./Results/toi.csv")
n_schedule = [10, 100, 1000, 10000, int(1e5), int(1e6), int(1e7)]

#plot maximum relative error against h.
plt.figure(figsize=(10,10))
plt.xlabel("log h")
plt.ylabel("maximum relative error")

for n in n_schedule:
    filter_n = toi['n'] == n
    #max_err_spec = (np.float64(toi[filter_n]['relative error spec'])).max()
    max_err_gen = (np.float64(toi[filter_n]['relative error gen'])).max()
    print("general case relative error is", max_err_gen)
    #plt.scatter(np.log(toi[filter_n]['h'][0]), max_err_spec,
    #         color = 'r', label ='special')
    h = 1/(n+1)
    print("h is:", h)
    print("log10 h :", np.log10(h))
    plt.scatter(np.log10(h), max_err_gen,
             color = 'b', label ='general')
plt.legend()
plt.savefig("./Results/h_vs_error.png")
plt.show()