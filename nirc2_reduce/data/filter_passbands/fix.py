#!/usr/bin/env python
import numpy as np

stuff = np.loadtxt('he1a_u.csv', delimiter = ',').T
#stuff[1] *= 100.0
stuff[0] *= 0.001
np.savetxt('he1_a.csv', stuff.T, delimiter = ',')