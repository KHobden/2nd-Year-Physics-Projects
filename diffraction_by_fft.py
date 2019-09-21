#Diffraction by FFT
#Kieran Hobden
#01-Mar-'19

import numpy as np
import matplotlib.pyplot as plt


"""Task 1"""

def fourier_transform(func, domain, dx):
    """Returns array with DFT and the frequency sampling of the DFT"""
    dft = np.fft.fft(func) #Discrete Fourier Transform
    dft = np.fft.fftshift(dft)
    dft_freq = np.fft.fftfreq(len(domain), dx) #Discrete Fourier Transform sample frequencies
    dft_freq = np.fft.fftshift(dft_freq)
    return dft, dft_freq

def grating(n, d, dx, k=0.0, D=1.0):
    """Generates the modified aperture of a grating"""
    aperture = np.zeros(2 ** n, dtype=np.complex64)
    centre = len(aperture) // 2
    width_of_slit = int(d / dx)
    for i in range(width_of_slit):
            correction = np.exp(1j * k * (dx*(i-(width_of_slit/2)))**2/(2*D)) #Correction for the near field
            aperture[centre - width_of_slit // 2 + i] = correction
    return aperture


"""Define constants"""
wavelength = 5*10**-7
k = 2*np.pi/wavelength
d = 10**-4
L = 0.005
D = 1
#We can add multiple slits here for some general aperture
n = 10 #Number of elements is 2**n for divide and conquer algorithm
dx = L / (2**n)

x = np.linspace(-L/2, L/2, 2**n)
j = np.linspace(-len(x)/2, len(x)/2, 1)
aperture = grating(n, d, dx)

"""Plot the results in with the one above the other. We plot the top graph first"""
top_graph_1 = plt.subplot(211) #Warnings about reusing axes labels appear without this line
top_graph_1.plot(x, abs(aperture), lw=2)
top_graph_1.set_title("Aperture function")
top_graph_1.set_xlabel("x (metres)")
top_graph_1.set_ylabel("Aperture function")

[psi, u] = fourier_transform(aperture, x, dx) #Psi is the DFT of the aperture and 2πu = ky/D
y = u * 2 * np.pi * D / k

"Bottom graph"""
bottom_graph_1 = plt.subplot(212)
bottom_graph_1.plot(y,(abs(psi)/np.amax(abs(psi)))**2, "b", lw=2, label='FFT intensity')
bottom_graph_1.plot(y,(abs(np.sinc(d*k*y/(2*np.pi*D))))**2, "g", lw=2, label = 'Theoretical intensity')
bottom_graph_1.set_title("FFT and theoretical intensities")
bottom_graph_1.set_xlabel("y (metres)")
bottom_graph_1.set_ylabel("I / I max")
bottom_graph_1.legend(loc="best")

plt.tight_layout()
plt.show()
plt.savefig(fname="general_fourier_transform.pdf")



"""Task 2"""

def sinusoidal_phase_grating(x, m, s, k=0, D=1):
    """Generates the modified aperture of a grating with sinusoidally varying phase"""
    phi = (m/2) * np.sin(2 * np.pi * x / s) #Phase changes sinusoidally with x
    aperture = np.exp(1j * phi)* np.exp(1j * k * x**2 / D) #Near field correction can be ignored
    return aperture

"""Define constants"""
L = 0.002
D = 10
n = 9 
m = 8
s = 10**-4

x = np.linspace(-L/2, L/2, 2 ** n)
dx = L / (2**n)
aperture = sinusoidal_phase_grating(x, m, s)

[psi, u] = fourier_transform(aperture,x,dx)
y = 2 * np.pi * (D / k) * u

"""Top graph"""
top_graph_2 = plt.subplot(211)
top_graph_2.plot(x, np.real(aperture), "b", lw=2, label='Real component')
top_graph_2.plot(x, np.imag(aperture), "g", lw=2, label='Imaginary component')
top_graph_2.set_title("Aperture function")
top_graph_2.set_xlabel("x (metres)")
top_graph_2.set_ylabel("Aperture function phase")
top_graph_2.legend(loc="best")

"""Bottom graph"""
bottom_graph_2 = plt.subplot(212)
bottom_graph_2.plot(y,(abs(psi)/np.amax(abs(psi)))**2, lw=2)
bottom_graph_2.set_title("Fraunhofer relative intensity")
bottom_graph_2.set_xlabel("y (metres)")
bottom_graph_2.set_ylabel("I / I max")

plt.tight_layout()
plt.show()
plt.savefig(fname="sinusoidal_phase_grating_plot.pdf")



"""Supplementary Task"""

"""Near field intensity patterns with D=5mm"""
d = 10**-4
L = 0.005
D = 0.005
n = 9
x = np.linspace(-L/2, L/2, 2 ** n)
dx = L/(2**n)
aperture = grating(n, d, dx, k, D)
[psi, u] = fourier_transform(aperture, x, dx)
y = 2 * np.pi * (D / k) * u

"""Top graph"""
top_graph_3 = plt.subplot(211)
top_graph_3.plot(y,(abs(psi)/np.amax(abs(psi)))**2, lw=2)
top_graph_3.set_title("Near-field finite slit relative intensity pattern")
top_graph_3.set_xlabel("y (metres)")
top_graph_3.set_ylabel("I / I max")

"""Near field intensity patterns with D=0.5m"""
D = 0.5
L = 0.002
m = 8
s = 10**-4
n = 9
x = np.linspace(-L/2, L/2, 2**n)
dx = L/(2**n)
aperture = sinusoidal_phase_grating(x, m, s, k, D)
[psi, u] = fourier_transform(aperture, x, dx)
y = 2 * np.pi * (D / k) * u

"""Bottom graph"""
bottom_graph_3 = plt.subplot(212)
bottom_graph_3.plot(y, (abs(psi) / np.amax(abs(psi))) ** 2, lw=2)
bottom_graph_3.set_title('Near-field phase grating relative intensity pattern')
bottom_graph_3.set_xlabel("y (metres)")
bottom_graph_3.set_ylabel("I / I max")

plt.tight_layout()
plt.show()
plt.savefig(fname="near_field_grating.pdf")
