# -*- coding: utf-8 -*-
"""

@author: jamesbc
"""
import numpy as np

def generate_x_and_fx(h, n):
    '''
    A function to generate x points and f(x)
    
    Inputs: 
        h: step length
        n: nbr of points
    Outputs
        x: array of length n+1
        fx: array of length n+1
    '''
    x = np.zeros((n+1))
    fx = np.zeros((n+1))
    for i in range(0, n):
        x[i] = i*h
        fx[i] = f_x(x[i])
        
    return x, fx

def f_x(x):
    '''
    returns f(x) from input x, float
    '''
    
    return 100*np.exp(-10*x)
    
def forward_sub_special_case(d_1):
    '''
    This function is incomplete!
    
    This function returns d_tilda in the case of 
    d_1 equals 2 and e_1 equals -1.
    
    Inputs:
        d_1 : an integer
    Outputs:
        d_tilda : an array of length n-1
    '''
    # initialise d tilda and first element
    d_tilda = np.zeros((n+1))
    d_tilda[1] = d_1
    for i in range(2, n-1):
        d_tilda[i] = (i+1)/i
    
    return d_tilda

def forward_sub(d_1, e_1, n, h, hh):
    '''
    This function iteratively calculates d tilda and g tilda
    g_tilda is equal to h**2*f(x_i). In the report it is written as 
    b_tilda. 
    
    Inputs: 
        d_1 : an integer
        e_1 : an integer
        n : an integer, the nbr of points.
    Outputs:
        d_tilda : an array of length n+1
        g_tilda : an array of length n+1
    '''
    x, fx = generate_x_and_fx(h, n)
    g = hh * fx
    
    d = d_1 * np.ones((n+1)) 
    e = e_1 * np.ones((n+1)) 
    d_tilda = np.zeros((n+1))
    g_tilda = np.zeros((n+1))
    
    d_tilda[1] = d[1]
    g_tilda[1] = g[1]
    
    for i in range(2, n-1):
        d_tilda[i] = d[i] - e[i-1]**2 * 1/(d_tilda[i-1])
        g_tilda[i] = g[i] - (e[i-1]*g_tilda[i-1]/d_tilda[i-1])
        
    return d_tilda, g_tilda, e
    
def backward_sub(d_tilda, g_tilda, e):
    '''
    This function backwards substitutes to calculate the 
    discretized v(x).
    
    '''
    v = np.zeros((n+1))
    v[n] = g_tilda[n] / d_tilda[n]
    
    for i in range(n-2, 1):
        v[i] = (g_tilda[i] - e[i]*v[i+1])*d_tilda[i]
        
    return v
# first we should initialise the component of the diagonal
# and off-diagonal elements and the number of points n.
d_1 = 2
e_1 = -1
n = 10

# calculate the step length, h from n.
h = 1/(n+1)
hh = h*h

d_tilda, g_tilda, e = forward_sub(d_1, e_1, n, h, hh)
print('d tilda is',d_tilda)
print('g tilda is', g_tilda)

v = backward_sub(d_tilda, g_tilda, e)
print(v)







