 #! /usr/bin/env python

import math
import scipy.constants

#------Physical Constants-------------------------------------------------------

c0 = scipy.constants.speed_of_light
mu0  = scipy.constants.mu_0
eps0 = scipy.constants.epsilon_0

#------Grid Conditions----------------------------------------------------------

L = 10
dx = 0.5
nP = [L/dx, L/dx-1]
sP = numPoints/2
CLFN = 0.8
dt = CFLN*dx/(math.sqrt(3)*c0)
totalTime = L/c0*2
timeSteps = totalTime // dt
spread = 1.0/math.sqrt(2.0)

#------Update Coefficients------------------------------------------------------

cE = dt / (eps0 * dx)
cH = dt / (mu0  * dx)

#------Initial Conditions-------------------------------------------------------

for i in range(0,numPoints[0]):

    e[i] = 0.0

for i in range(0, numPoints[1]):
    h[i] = 0.0

#------Temporal Loop------------------------------------------------------------

t = 0.0

for i in range(0, timeSteps):
   
    for j in range (0,nP[1]):

        h[j] = h[j] + cH * (e[j+1] - e[j])     

    for j in range (0,nP[1]):

        e[j] = e[j] + cE * (h[j] - h[j-1])    

    e[sP] = e[sP] + math.exp(-0.5*math.pow((t-t0)/spread, 2))   


    t += dt
    

