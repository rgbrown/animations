#!/usr/bin/env python3
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint

plt.style.use("dark_background")
# This file creates an animation movie of the Lorenz attractor. Black background with developing trails

# Simulation parameters
fps = 30
pts_per_frame = 5
sim_time = 40 # simulation time (seconds)
speed_factor = 0.5
trail_time = 1 # seconds

# Define Lorenz equations
rho = 28.0
sigma = 10.0
beta = 8.0 / 3.0


def f(u, t):
    x, y, z = u
    return speed_factor * sigma * (y - x), speed_factor * (x * (rho - z) - y), speed_factor * (x * y - beta * z)


# The plan. Start with initial conditions scattered throughout and solve Lorenz equations. Solutions for each trace will be in X[k], Y[k], Z[k]
# [-20, 20] x [-20, 20] x [0, 40]
n_traces = 30
X0 = -20 + 40 * np.random.rand(n_traces)
Y0 = -20 + 40 * np.random.rand(n_traces)
Z0 = 40 * np.random.rand(n_traces)

solns = []
t = np.arange(0, sim_time, 1/(pts_per_frame * fps))
for k in range(n_traces):
    u0 = [X0[k], Y0[k], Z0[k]]
    solns.append(odeint(f, u0, t))

# Add nans to the start for initial traces
n_trail = pts_per_frame*fps*trail_time
nans = np.empty(n_trail)
nans.fill(np.nan)
X = np.array([np.concatenate((nans, solns[k][:, 0])) for k in range(n_traces)])
Y = np.array([np.concatenate((nans, solns[k][:, 1])) for k in range(n_traces)])
Z = np.array([np.concatenate((nans, solns[k][:, 2])) for k in range(n_traces)])

# We're going to have trails one second long
n_frames = fps*sim_time

# Draw initial traces
fig = plt.figure(figsize=(12, 9))
ax = plt.axes((0, 0, 1, 1), projection="3d", xlim=(-20, 20), ylim=(-20, 20), zlim=(0, 40))
traces = []
for k in range(n_traces):
    traces.append(ax.plot(nans, nans, nans)[0])
plt.axis('off')

def animate(i):
    i_start = pts_per_frame*i
    for k in range(n_traces):
        traces[k].set_xdata(X[k][i_start:(n_trail+i_start)])
        traces[k].set_ydata(Y[k][i_start:(n_trail+i_start)])
        traces[k].set_3d_properties(Z[k][i_start:(n_trail+i_start)])
    ax.view_init(45*np.sin(0.2*t[i_start]), 180*np.cos(0.1*t[i_start]))
    ax.dist = 9

    return traces

ani = FuncAnimation(fig, animate, n_frames, interval=1000/fps, blit=False)

ani.save("lorenz.mp4", fps=fps)



plt.draw()
plt.show()
