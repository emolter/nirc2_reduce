import numpy as np

sun = np.loadtxt('sun_spectrum.txt')
w = sun.T[0]
f = sun.T[1]

f = f/1e4

