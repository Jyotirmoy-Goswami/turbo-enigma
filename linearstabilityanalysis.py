# -*- coding: utf-8 -*-
"""
Created on Thu May  5 18:00:45 2022

@author: Jyotirmoy Goswami
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import pylab as pl

a, b, c, d = -1, -1, 1, -1
def dx_dt(x, t):
    return [a*x[0] + b*x[1], c*x[0] + d*x[1]]

ts = np.linspace(0, 4, 100)
ic = np.linspace(-3, 3, 5)
for r in ic:
    for s in ic:
        x0 = [r, s]
        xs = odeint(dx_dt, x0, ts)
        plt.plot(xs[:,0], xs[:,1], "r-")
ts = np.linspace(0, -4, 100)
ic = np.linspace(-3, 3, 5)
for r in ic:
    for s in ic:
        x0 = [r, s]
        xs = odeint(dx_dt, x0, ts)
        plt.plot(xs[:,0], xs[:,1], "b-")
        
plt.xlabel('x', fontsize=15)
plt.ylabel('y', fontsize=15)
plt.tick_params(labelsize=15)
plt.xlim(-3, 3)
plt.ylim(-3, 3)

X,Y = np.mgrid[-3:3:10j, -3:3:10j]
u= a*X + b*Y
v= c*X + d*Y
pl.quiver(X, Y, u, v, color = 'g')
plt.show()