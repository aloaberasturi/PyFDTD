 #! /usr/bin/env python

import math
import scipy.constants
import numpy 
import matplotlib.pyplot as plt
import matplotlib.animation as animation


#------Physical Constants-------------------------------------------------------

c0   = scipy.constants.speed_of_light
mu0  = scipy.constants.mu_0
eps0 = scipy.constants.epsilon_0
imp0 = math.sqrt(mu0 / eps0)

#-------Source Conditions-------------------------------------------------------

delay  = 8e-9
spread = 2e-9

#------Grid Conditions----------------------------------------------------------

L         = 10.0
dx        = 0.05
nP        = int(L / dx)+1
sP        = int((L/8)/dx)
CFLN      = 0.99
dt        = CFLN * dx / c0
totalTime = L / c0 * 2 
timeSteps = int(totalTime / dt)
nSamples  = int( math.floor(totalTime/dt) )


#------Update Coefficients------------------------------------------------------

cE = dt / (eps0 * dx)
cH = dt / (mu0  * dx)

#------Initial Conditions-----------------------------------------------

gridE     = numpy.linspace(0,      L,        num=L/dx+1, endpoint=True)
gridH     = numpy.linspace(dx/2.0, L-dx/2.0, num=L/dx,   endpoint=True)
e         = numpy.zeros(nP)
h         = numpy.zeros(nP-1)

#------Probes-------------------------------------------------------------------
probeE    = numpy.zeros((nP, nSamples))
probeH    = numpy.zeros((nP-1, nSamples))
probeTime = numpy.zeros(nSamples)

#------Temporal Loop------------------------------------------------------------

t  = 0.0

for i in range(timeSteps):
   
    for j in range (nP-1):

        h[j] = h[j] + cH * (e[j+1] - e[j])   

    # -- Magnetic Source --

    h[sP-1] = h[sP-1] - cH * math.exp(-0.5*math.pow((t-delay)/spread, 2)) 

    for j in range (1, nP-1):

        e[j] = e[j] + cE * (h[j] - h[j-1])        
    
    # -- Electric Source --  

    e[sP] = e[sP] + cE * math.exp(-0.5*math.pow((t- delay + dx/(2*c0) + dt/2)/spread, 2)) / imp0

    # PEC
    e[ 0] = 0.0
    e[-1] = 0.0

    probeE[:,i]  = e[:]
    probeH[:,i]  = h[:]
    probeTime[i] = t 
    t           += dt
    


# ==== Post-processing ========================================================

# --- Creates animation ---
fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(1, 2, 1)
ax1 = plt.axes(xlim=(gridE[0], gridE[-1]), ylim=(-1.1, 1.1))
ax1.grid(color='gray', linestyle='--', linewidth=.2)
ax1.set_xlabel('X coordinate [m]')
ax1.set_ylabel('Field')
line1,    = ax1.plot([], [], 'o', markersize=1)
timeText1 = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)

ax2 = fig.add_subplot(2, 2, 2)
ax2 = plt.axes(xlim=(gridE[0], gridE[-1]), ylim=(-1.1, 1.1))
ax2.grid(color='gray', linestyle='--', linewidth=.2)
# ax2.set_xlabel('X coordinate [m]')
# ax2.set_ylabel('Magnetic field [T]')
line2,    = ax2.plot([], [], 'o', markersize=1)
timeText2 = ax2.text(0.02, 0.95, '', transform=ax2.transAxes)

def init():
    line1.set_data([], [])
    timeText1.set_text('')
    line2.set_data([], [])
    timeText2.set_text('')
    return line1, timeText1, line2, timeText2

def animate(i):
    line1.set_data(gridE, probeE[:,i])
    timeText1.set_text('Time = %2.1f [ns]' % (probeTime[i]*1e9))
    line2.set_data(gridH, probeH[:,i]*100)
    timeText2.set_text('Time = %2.1f [ns]' % (probeTime[i]*1e9))
    return line1, timeText1, line2, timeText2

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=nSamples, interval=50, blit=True)

plt.show()

print('=== Program finished ===')