#!/usr/bin/env python3
# Pendulum animation - solve nonlinear pendulum equation in polar coordinates, then animate solution, including velocity and acceleration vectors. 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint

plt.style.use("dark_background")

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
ax = plt.axes(xlim=(-10, 10), ylim = [-12, 0], aspect='equal')

shaft, = plt.plot([0, 0], [0, 0], color='C0', linewidth=3)
mass, = plt.plot([0], [0], 'o', color='C0', markersize=20)
velocity = plt.arrow(0, 0, 0, 0, lw=3, color='C1', head_width=0.15, label='velocity')
acceleration = plt.arrow(0, 0, 0, 0, lw=3, color='C3', head_width=0.15,label='acceleration')
v_text = plt.text(0, 0, 'v', color='C1', fontweight='bold', fontsize=20)
a_text = plt.text(0, 0, 'a', color='C3', fontweight='bold', fontsize=20)

plt.axis('off')


# Define animation function
def animate(i):
    x, v, a = polar_to_cart(U[i, 0], U[i, 1])
    shaft.set_xdata([0, x[0]])
    shaft.set_ydata([0, x[1]])
    mass.set_xdata([x[0]])
    mass.set_ydata([x[1]])
    velocity.set_data(x=x[0], y=x[1], dx=v[0], dy=v[1])
    acceleration.set_data(x=x[0], y=x[1], dx=a[0], dy=a[1])
    v_text.set_position(x + (np.linalg.norm(v) + 1)/np.linalg.norm(v)*v)
    a_text.set_position(x + (np.linalg.norm(a) + 1)/np.linalg.norm(a)*a)
    return shaft, mass, velocity, v_text, a_text

ani = FuncAnimation(fig, animate, n_frames, interval=1000/fps, blit=False)
plt.show()
    




