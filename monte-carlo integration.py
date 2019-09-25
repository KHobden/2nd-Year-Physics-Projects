# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 14:34:19 2019

@author: kjh57
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import math

def calculate(N): 

    total_sines = 0
    total_sines_squared = 0

    """Volume element"""
    V = (np.pi/8)**8

    x_coords = np.random.random((8, N))
    for i in range(N):
        """Pick random coords for a point in 8D space"""
        x_coords = np.random.random(8) * (np.pi / 8)
        
        """Finds <f> * N and then <f^2> * N
        No need to divide by N each time"""
        total_sines += np.sin(np.sum(x_coords))
        total_sines_squared += (np.sin(np.sum(x_coords))**2)
    
    """Finds <f> and then <f^2> (written as total_sines
    and total_sines_squared respectively"""
    total_sines = total_sines / N
    total_sines_squared = total_sines_squared / N

    """Produces mean value and theoretical error"""
    mean = (10**6) * total_sines * V
    sigma = (10**6) * V * np.sqrt((total_sines_squared - (total_sines ** 2)) / N)
    return mean, sigma

"""Logarithmically even spaced points""" 
num_samples = 10
Ns = np.logspace(1, 4, num_samples, dtype=int)

"""Calculate the mean of the Monte-Carlo simulations and
the error in its scatter and the theoretical error"""
means = np.zeros(num_samples)
error_in_mc_scatter = np.zeros(num_samples)
theor_error = np.zeros(num_samples)
for i in range(num_samples):
    calculation_results = np.zeros(25)
    monte_means = np.zeros(25)
    monte_theor_error = np.zeros(25)
    for j in range(25):
        monte_means[j] = calculate(Ns[i])[0]
        monte_theor_error[j] = calculate(Ns[i])[1] 
    means[i] = np.mean(monte_means)
    error_in_mc_scatter[i] = np.std(monte_means)
    theor_error[i] = np.mean(monte_theor_error)

"""Plot graph with error bars"""
plt.errorbar(Ns, means, yerr = error_in_mc_scatter, color = '#000000', ecolor = '#000000', elinewidth = 0.5, capsize = 2)
plt.xscale('log')
plt.title('Monte-Carlo integration')
plt.xlabel('Number of iterations')
plt.ylabel('Integral value')
plt.savefig(fname = "plot.pdf")
plt.show()

"""Plot of error of Monte-Carlo scatter against theoretical error"""
plt.plot(Ns, error_in_mc_scatter, "b", label="Error in Monte_Carlo scatter")
plt.plot(Ns, theor_error, "g", label="Theoretical error")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="best")
plt.title("Theoretical error and error in Monte-Carlo scatter")
plt.xlabel("Number of iterations")
plt.ylabel("Error in calculation")
plt.savefig(fname = "errorplot.pdf")
plt.show()

print("Final value =", means[num_samples - 1], "+/-", error_in_mc_scatter[num_samples - 1])
print(np.polyfit(np.log(Ns), np.log(error_in_mc_scatter), 1)[0])


"""Cornu Spiral"""
"""No. of steps in the plot"""
n = 1000

"""Upper limit of integral"""
u_max = 5

"""U values"""
u_vals = np.linspace(-u_max, u_max, num=n)

"""Set up coordinate value arrays"""
x_val = np.zeros(n)
y_val = np.zeros(n)

"""Compute coordinate values"""
tally = 0
for i in u_vals:
    fresnel_c = integrate.quad(lambda x: np.cos(np.pi * x**2/2), 0, i)
    fresnel_s = integrate.quad(lambda x: np.sin(np.pi * x**2/2), 0, i)
    x_val[tally] = fresnel_c[0]
    y_val[tally] = fresnel_s[0]
    tally += 1

"""Set up marker coordinate values"""
x_marker_val = np.zeros(2*u_max+1)
y_marker_val = np.zeros(2*u_max+1)

"""Compute markers"""
tally2 = 0
for i in range(-u_max, u_max+1):
    fresnel_c = integrate.quad(lambda x: np.cos(np.pi * x**2/2), 0, i)
    fresnel_s = integrate.quad(lambda x: np.sin(np.pi * x**2/2), 0, i)
    x_marker_val[tally2] = fresnel_c[0]
    y_marker_val[tally2] = fresnel_s[0]
    tally2 += 1


"""Plot Cornu Spiral"""
plt.plot(x_val, y_val, linewidth=1, color="#000000")
plt.plot(x_marker_val, y_marker_val, linewidth=0, marker="x", color="#000000")
plt.title("Cornu Spiral")
plt.xlabel("C(u)")
plt.ylabel("S(u)")
plt.savefig(fname = "fresnel_plot.pdf")
plt.show()

