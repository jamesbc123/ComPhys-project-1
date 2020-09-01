# -*- coding: utf-8 -*-
"""

@author: jamesbc
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def generate_x_and_fx(n, h):
    '''
    A function to generate x points and f(x)
    Generates arrays of length n+2, since x_0 and 
    x_n+1 exist, but equal 0.
    Inputs: 
        h: step length
        n: nbr of points
    Outputs
        x: array of length n+2
        fx: array of length n+2
    '''
    x = np.zeros((n+2))
    fx = np.zeros((n+2))
    for i in range(0, n+1):
        x[i] = i*h
        fx[i] = f_x(x[i])
        
    return x, fx

def f_x(x):
    '''
    returns f(x) from input x
    '''
    return 100*np.exp(-10*x)

def calc_exact(x):
    len_x = len(x)
    u = np.zeros((len_x))
    for i in range(1, len_x-1):
        u[i] = exact(x[i])
    return u 

def exact(x):
    '''
    This function returns the exact value
    '''
    return 1.0-(1-np.exp(-10))*x-np.exp(-10*x)

def general_algo(n, hh, fx, a, b, c):
    '''
    This function calculates v for the general case, where
    a, b and c are different values. This is done by forward 
    and backward substitution. 
    
    Inputs: 
        a, b, c are arrays of length n+1
        n: the nbr of points
        hh: the step length squared
    Outputs:
        v: solution, array of length n+2 
    '''
    b_sub = np.zeros((n+1))
    g_sub = np.zeros((n+1))
    g = hh*fx
    b_sub[1] = b[1]
    g_sub[1] = g[1]
    v = np.zeros((n+2))
    print("general")
    for i in range(2, n+1):
        ab = a[i-1]/b_sub[i-1]
        b_sub[i] = b[i] - (c[i-1]*ab)
        g_sub[i] = g[i] - (g_sub[i-1]*ab)
        print("g is:", g_sub[i])
    v[n] = g_sub[n]/b_sub[n]
    
    for i in range(n-1, 0, -1):
        v[i] = (g_sub[i] - c[i]*v[i+1])/b_sub[i] 
        print("v is:", v[i])
    return v

def special_algo(n, hh, fx, a, b):
    '''
    This function calculates v for a 
    tri-diagonal matrix. g is equal to h**2*f(x_i), in the report 
    it is written as b_tilda, but this makes it confusing since 
    here we have b. Therefore it will be g here.
    
    Inputs: 
        a, b are arrays of length n+1
        n: the nbr of points
        hh: the step length squared
    Outputs:
        v: solution, array of length n+2 
    '''
    
    b_sub = np.zeros((n+1))
    g_sub = np.zeros((n+1))
    g = hh*fx
    b_sub[1] = b
    g_sub[1] = g[1]
    v = np.zeros((n+2))
    
    for i in range(2, n+1):
        b_sub[i] = (i+1)/i
        #print("b", b_sub[i])
        g_sub[i] = g[i] - (a*g_sub[i-1]/b_sub[i-1])
        print("g", g_sub[i])
    v[n] = g_sub[n] / b_sub[n]
    
    for i in range(n-1, 0, -1):
        v[i] = (g_sub[i] - a*v[i+1])*b_sub[i]
        print("v is", v[i])
    return v

def rel_error(u, v):
    '''
    This function compares the exact value to the 
    numerical value.
    
    Inputs:
        u: an array of the exact values.
        v: an array of numerical values.
    Outputs:
        rel_err: an array of relative errors.
    '''
    rel_err = np.zeros((n+2))
    for i in range(1, n+1):
        rel_err[i] = np.log10(np.abs((u[i] - v[i])/u[i]))
        
    # there was a problem with the boundary points if the error
    # is equal to 0 since it is logged. Therefore manually
    # it has been set to a large negative value for 
    # practical reasons. 
    rel_err[0] = rel_err[-1] = -1e6
    return rel_err

# first we should initialise the component of the diagonal
# and off-diagonal elements and the number of points n.
a = -1
b = 2
c = -1
n_schedule = [10]

toi = pd.DataFrame(columns=["n", "h", "x", "exact", "v special", "v general",
                            "relative error spec",
                            "relative error gen"])

for n in n_schedule:
    h = np.float64(1/(n+1))
    hh = h*h
    
    x, fx = generate_x_and_fx(n, h)
    u = calc_exact(x)
    
    v_spec = special_algo(n, hh, fx, a, b)
    
    v_gen = general_algo(n, hh, fx, a*np.ones((n+1)), 
                     b*np.ones((n+1)), c*np.ones((n+1))
                     )
    rel_err_spec = rel_error(u, v_spec)
    rel_err_gen = rel_error(u, v_gen)
    
    # add all this info to a csv file. 
    dat = np.array([[n for i in range(n+2)],
                    [h for i in range(n+2)],
                    x, u, v_spec, v_gen,
                    rel_err_spec, rel_err_gen
                    ])
    
    temp = pd.DataFrame(dat.T, columns = ["n", "h", "x", "exact", "v special",
                                          "v general",
                            "relative error spec",
                            "relative error gen"])
    
    toi = toi.append(temp)


# save to csv
toi.to_csv('./Results/toi.csv')

#plot maximum relative error against h.
plt.figure(figsize=(10,10))
plt.xlabel("log h")
plt.ylabel("maximum relative error")

for n in n_schedule:
    filter_n = toi['n'] == n
    max_err_spec = (np.float64(toi[filter_n]['relative error spec'])).max()
    max_err_gen = (np.float64(toi[filter_n]['relative error gen'])).max()
    print("general case relative error is", max_err_gen)
    plt.scatter(np.log10(toi[filter_n]['h'][0]), max_err_spec,
             color = 'r', label ='special')
    plt.scatter(np.log10(toi[filter_n]['h'][0]), max_err_gen,
             color = 'b', label ='general')
plt.legend()
plt.savefig("./Results/h_vs_error.png")
plt.show()


# plot x against exact and v
for n in n_schedule:
    plt.figure(figsize=(10,10))
    plt.xlabel("x")
    plt.ylabel("numerical or exact solution")

    filter_n = toi['n'] == n
    plt.scatter(toi[filter_n]['x'], toi[filter_n]['exact'],
             color = 'black', label ='exact')
    plt.scatter(toi[filter_n]['x'], toi[filter_n]['v general'],
             color = 'b', label ='general')
    plt.scatter(toi[filter_n]['x'], toi[filter_n]['v special'],
             color = 'r', label ='special')
    plt.legend()
    plt.savefig("./Results/solution_vs_exact_n="+ str(n) +".png")
    plt.close()
 
