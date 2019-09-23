#Fresnel Integrals
#Kieran Hobden
#19-Feb-'19

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

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
