 #! /usr/bin/env python

import math
import scipy.constants
import numpy

#------Physical Constants-------------------------------------------------------

c0 = scipy.constants.speed_of_light
mu0  = scipy.constants.mu_0
eps0 = scipy.constants.epsilon_0

#------Grid Conditions----------------------------------------------------------

L = 10
dx = 0.5
nP = [int(L // dx),int(L // dx)-1]
sP = nP[0]/2
CFLN = 0.8
dt = CFLN*dx/(math.sqrt(3)*c0)
totalTime = L/c0*2
timeSteps = int(totalTime / dt)
spread = 1.0/math.sqrt(2.0)

#------Update Coefficients------------------------------------------------------

cE = dt / (eps0 * dx)
cH = dt / (mu0  * dx)

#------Initial and PEC Conditions-------------------------------------------------------

e = numpy.zeros(nP[0])
h = numpy.zeros(nP[1])

#------Temporal Loop------------------------------------------------------------
t0 = 0.0
t  = 0.0

for i in range(0, timeSteps):
   
    for j in range (0,nP[1]):

        h[j] = h[j] + cH * (e[j+1] - e[j])     

    for j in range (0,nP[1]):

        e[j] = e[j] + cE * (h[j] - h[j-1])    

    e[sP] = e[sP] + math.exp(-0.5*math.pow((t-t0)/spread, 2))   
    t += dt
    

