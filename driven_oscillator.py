#The Driven Pendulum
#Kieran Hobden
#24-Feb-'19

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

"""Pendulum function"""
def pend(y, t, q, F):
    theta, omega = y
    dydt = [omega, -np.sin(theta) - q*omega + F*np.sin(2*t/3)]
    return dydt

"""Define constants"""
q = 0.5
F = 1.465

"""Set initial conditions"""
#General starting position in phase space - useful in theoretical calculation
y_0 = [F * np.sin(np.arctan((2/3) * q / ((2/3)**2 - 1))) / np.sqrt(((2/3) * q)**2 + (1 - (2/3)*(2/3))**2), F*(2/3)*np.cos(np.arctan((2/3)*q/((2/3)**2-1)))/np.sqrt(((2/3)*q)**2+(1-(2/3)**2)**2)]
#Specific starting position in phase space
#y_0 = [0.01, 0.0]

"""Coordinates to plot"""
time_period = 2*np.pi
t = np.linspace(0, 20*time_period, 10000)

"""Integrate ODE"""
sol = odeint(pend, y_0, t, args=(q, F))

"""Theoretical solution"""
theor_theta = F * np.sin((2/3)*t+np.arctan((2/3)*q/((2/3)**2-1)))/np.sqrt(((2/3)*q)**2+(1-(2/3)**2)**2)
theor_omega = F*(2/3)*np.cos((2/3)*t+np.arctan((2/3)*q/((2/3)**2-1)))/np.sqrt(((2/3)*q)**2+(1-(2/3)**2)**2)
#theor_energy = theor_theta**2 + theor_omega**2

"""Plot solution"""

plt.plot(t, sol[:, 0], "b", label = "Angle")
plt.plot(t, sol[:, 1], "g", label = "Angular Velocity")
plt.plot(t, 0.5*sol[:, 1]**2 + (1 - np.cos(sol[:, 0])), "r", label = "Energy")
plt.legend(loc = "best")
plt.title("Oscillatory motion with time (F=1.465)")
plt.xlabel("Time (seconds)")
plt.grid()
plt.show()

"""Theoretical solution"""
"""
plt.plot(t, theor_theta, label = "Theoretical angle")
plt.plot(t, theor_omega, label = "Theoretical angular velocity")
#plt.plot(t, theor_energy, label = "Theoretical energy")
plt.xlabel("Time (seconds)")
plt.legend(loc="best")
plt.show()
"""

"""Solution and theoretical solution"""
"""
plt.plot(t, sol[:, 0], "b", label = "Angle")
plt.plot(t, theor_theta, "r", label = "Theoretical angle")
plt.title("Angle of Oscillation with time (10 oscillations)")
plt.xlabel("Time (seconds)")
plt.legend(loc="best")
plt.show()
"""

"""Plot of the energy of the system with time"""
"""
plt.plot(t, 0.5*sol[:, 1]**2 + (1 - np.cos(sol[:, 0])), "r", label = "Energy")
plt.title("Energy of the system with time")
plt.xlabel("Time (seconds)")
plt.ylabel("Energy")
plt.show()
"""

"""Determining how amplitude affects the period"""
#As described in the readme.pdf, read the time periods and insert into the plot below
"""for N in range(0,10000):
	print(solution[N,1],N)
	N += 1"""
"""Plot of the change in period with amplitude"""
"""
x = np.linspace(1,8,8)*np.pi/16
y = [6.2983, 6.3442, 6.4222, 6.5344, 6.6841, 6.8760, 7.1165, 7.4154]
plt.plot(x , y)
plt.title("Period of oscillations against amplitude")
plt.xlabel("Initial displacement (radians)")
plt.ylabel("Period (seconds)")
plt.show()
"""
