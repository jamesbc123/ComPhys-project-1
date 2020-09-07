# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 12:56:53 2020

@author: jamesbc
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Change the font size of all figures.
plt.rcParams.update({'font.size': 22})

# Decide what to plot
plot_error = False
plot_num_and_exact = False
print_timing = True

if plot_error or plot_num_and_exact:
    toi = pd.read_csv("./Results/1.6GBtoi.csv")
n_schedule = np.array([10, 100, 1e3, 1e4, 1e5, 1e6, 1e7], dtype=int)
if plot_error:
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
                     color = 'tab:orange', marker = "^", label ='Special algorithm')
            plt.scatter(np.log10(h), max_err_gen,
                     color = 'b', label ='General algorithm')
            label_added = True
        else:
            plt.scatter(np.log10(h), max_err_spec,
                     color = 'tab:orange', marker = "^",)
            plt.scatter(np.log10(h), max_err_gen,
                     color = 'b')
    plt.legend()
    plt.savefig("./Results/error_vs_h.png")
    plt.show()

if plot_num_and_exact:
    # plot x against exact and v
    for n in n_schedule:
        plt.figure(figsize=(10,10))
        plt.xlabel("x")
        plt.ylabel("numerical or exact solution")
    
        filter_n = toi['n'] == n
        plt.scatter(toi[filter_n]['x'], toi[filter_n]['exact'],
                 color = 'k', marker = "1", label ='exact')
        plt.scatter(toi[filter_n]['x'], toi[filter_n]['v general'],
                 color = 'b', label ='General solution')
        plt.scatter(toi[filter_n]['x'], toi[filter_n]['v special'],
                 color = 'tab:orange', marker = "^", label ='Special solution')
        plt.legend()
        plt.savefig("./Results/solution_vs_exact_n="+ str(n) +".png")
        plt.close()
        
if print_timing:
    toi_timing = pd.read_csv("./Results/1.6GBtoiTiming.csv")
    print(toi_timing.to_latex(index = False))