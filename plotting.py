# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 12:56:53 2020

@author: jamesbc
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

plotNumVsSol = False

toi = pd.read_csv("./Results/toi.csv")
n_schedule = np.array([10, 100, 1e3, 1e4, 1e5], dtype=int)

#plot maximum relative error against h.
plt.figure(figsize=(10,10))
plt.xlabel("Log(h)")
plt.ylabel("Maximum Relative Error")

label_added = False
for n in n_schedule:
    filter_n = toi['n'] == n
    h = 1/(n+1)
    
    
    max_err_spec = (np.float64(toi[filter_n]['relative error spec'])).max()
    max_err_gen = (np.float64(toi[filter_n]['relative error gen'])).max()
    
    if not label_added:
        plt.scatter(np.log10(h), max_err_spec,
                 color = 'tab:orange', marker = "^", label ='Special solution')
        plt.scatter(np.log10(h), max_err_gen,
                 color = 'b', label ='General solution')
        label_added = True
    else:
        plt.scatter(np.log10(h), max_err_spec,
                 color = 'tab:orange', marker = "^",)
        plt.scatter(np.log10(h), max_err_gen,
                 color = 'b')
plt.legend()
plt.savefig("./Results/error_vs_h.png")
plt.show()

if plotNumVsSol:
    # plot x against exact and v
    for n in n_schedule:
        plt.figure(figsize=(10,10))
        plt.xlabel("x")
        plt.ylabel("numerical or exact solution")
    
        filter_n = toi['n'] == n
        plt.scatter(toi[filter_n]['x'], toi[filter_n]['exact'],
                 color = 'black', label ='exact')
        plt.scatter(toi[filter_n]['x'], toi[filter_n]['v general'],
                 color = 'b', label ='General solution')
        plt.scatter(toi[filter_n]['x'], toi[filter_n]['v special'],
                 color = 'r', label ='Special solution')
        plt.legend()
        plt.savefig("./Results/solution_vs_exact_n="+ str(n) +".png")
        plt.close()