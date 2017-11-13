#!/usr/bin/env python
import numpy as np

stuff = np.loadtxt('chrs_new.csv', delimiter = ',').T
stuff[1] *= 100.0
np.savetxt('ch4s.csv', stuff.T, delimiter = ',')