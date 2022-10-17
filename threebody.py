#!/usr/bin/env/python3
from mpl_toolkits import mplot3d
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint

plt.style.use("dark_background")
# Simulation parameters
fps = 30
pts_per_frame = 10 
sim_time = 180 # simulation time (seconds)
trail_time = 2 # seconds

# Three body problem parameters
m1, m2, m3 = 1., 1., 1.
G = 0.1 

def f(u, t):
    du = np.zeros(12)
    x1, x2, x3 = u[0:2], u[2:4], u[4:6]
    v1, v2, v3 = u[6:8], u[8:10], u[10:12]

    du[0:2], du[2:4], du[4:6] = v1, v2, v3
    du[6:8] = G*(m2 / norm(x2 - x1)**3 * (x2 - x1) + m3/norm(x3 - x1)**3 * (x3-x1))
    du[8:10] = G*(m1 / norm(x1 - x2)**3 * (x1 - x2) + m3/norm(x3 - x2)**3 * (x3-x2))
    du[10:12] = G*(m1 / norm(x1 - x3)**3 * (x1 - x3) + m2/norm(x2 - x3)**3 * (x2-x3))
    return du
    
# Initial conditions
x1_0 = np.array([0.5, 0.3])
x2_0 = np.array([-0.7, 0.9])
x3_0 = np.array([-0.1, -0.1])
v1_0 = 0*np.array([-0.1, -0.1])
v2_0 = 0*np.array([0.2, -0.1])
v3_0 = 0*np.array([-0.2, 0.1])

t = np.arange(0, sim_time , 1/(pts_per_frame * fps))
u0 = np.concatenate((x1_0, x2_0, x3_0, v1_0, v2_0, v3_0))
soln = odeint(f, u0, t)

# Add nans to the start for initial traces
n_trail = pts_per_frame * fps * trail_time
nans = np.empty((n_trail, 6))
nans.fill(np.nan)
U = np.vstack((nans, soln[:, 0:6]))

# Draw initial plots
fig = plt.figure(figsize=(12, 9))
ax = plt.axes(xlim = (-1.2, 1.2), ylim = (-0.5, 1.3), aspect='equal')
colors = ['C0', 'C1', 'C3']
traces = []
pts = []
for k in range(3):
    traces.append(ax.plot(nans[:, 2*k], nans[:, 2*k+1], color=colors[k], alpha=0.7)[0])
    pts.append(ax.plot(np.nan, np.nan, 'o', color=colors[k], markersize=10)[0])
plt.axis('off')

n_frames = fps*sim_time
def animate(i):
    i_start = pts_per_frame*i
    for k in range(3):
        traces[k].set_xdata(U[i_start:(n_trail+i_start), 2*k])
        traces[k].set_ydata(U[i_start:(n_trail+i_start), 2*k + 1])
        pts[k].set_xdata(U[n_trail+i_start, 2*k])
        pts[k].set_ydata(U[n_trail+i_start, 2*k + 1])

    return traces, pts

ani = FuncAnimation(fig, animate, n_frames, interval=1000/fps, blit=False)
ani.save("threebody.mp4", fps=fps)




plt.draw()
plt.show()



