#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint

# Pendulum animation - solve nonlinear pendulum equation in polar coordinates, then animate solution, including velocity and acceleration vectors. 

# Pendulum parameters
L = 7 
g = 9.8

# simulation parameters
fps = 30
T = 20
n_frames = fps*T
dt = 1/fps
t = np.arange(0, T, dt)

theta0 = 45*np.pi/180
omega0 = 0
u0 = [theta0, omega0]

# angle theta is measured anticlockwise from vertical (straight down is zero)
def pendulum(u, t):
    # unpack variables
    theta, omega = u
    return omega, -g / L * np.sin(theta)

# get position, velocity, acceleration vectors
def polar_to_cart(theta, omega):
    x = L*np.array([np.sin(theta), -np.cos(theta)])
    v = L*omega*np.array([np.cos(theta), np.sin(theta)])
    a = L*omega**2*np.array([-np.sin(theta), np.cos(theta)]) - g*np.sin(theta)*np.array([np.cos(theta), np.sin(theta)])
    return x, v, a

# Perform simulation
U = odeint(pendulum, u0, t)

# Create figure
fig = plt.figure(figsize=(12, 9))
ax = plt.axes(xlim=(-10, 10), ylim = [-12, 0])

shaft, = plt.plot([0, 0], [0, 0], linewidth=3)
mass, = plt.plot([0], [0], 'ko', markersize=20)
velocity = plt.arrow(0, 0, 0, 0, color='r', lw=3, head_width=0.15)
acceleration = plt.arrow(0, 0, 0, 0, color='g', lw=3, head_width=0.15)

# Define animation function
def animate(i):
    x, v, a = polar_to_cart(U[i, 0], U[i, 1])
    shaft.set_xdata([0, x[0]])
    shaft.set_ydata([0, x[1]])
    mass.set_xdata([x[0]])
    mass.set_ydata([x[1]])
    velocity.set_data(x=x[0], y=x[1], dx=v[0], dy=v[1])
    acceleration.set_data(x=x[0], y=x[1], dx=a[0], dy=a[1])
    return shaft, mass, velocity

ani = FuncAnimation(fig, animate, n_frames, interval=1000/fps, blit=False)
plt.show()
    




