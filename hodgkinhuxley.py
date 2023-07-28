#!/usr/bin/env python3
# Animation of Hodgkin-Huxley equations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection
from scipy.integrate import solve_ivp

plt.style.use("dark_background")

# Simulation parameters
dt = 0.05 #ms
t_max = 300
t_max_ax = 50
nt_ax = round(50/dt) + 1
t_eval = np.arange(0, t_max, dt)
nt = len(t_eval)

fps = 30

# Define Hodgkin Huxley equations
# Conductances
gNa = 120
gK = 36
gL = 0.3
Cm = 1
# Nernst potentials
vNa = 115
vK = -12
vL = 10.6 

def fun(t, u, Iapp):
    v, m, n, h = u

    alpha_m = 0.1*(25 - v) / (np.exp((25 - v)/10) - 1)
    beta_m  = 4*np.exp(-v/18)
    alpha_h = 0.07*np.exp(-v/20)
    beta_h = 1/(np.exp((30-v)/10)+1)
    alpha_n = 0.01*(10-v)/(np.exp((10-v)/10) - 1)
    beta_n = 0.125*np.exp(-v/80)
    
    dv = 1/Cm * (-gK*n**4*(v - vK) - gNa*m**3*h*(v - vNa) - gL*(v - vL)) + Iapp
    dm = alpha_m * (1 - m) - beta_m * m
    dn = alpha_n * (1 - n) - beta_n * n
    dh = alpha_h * (1 - h) - beta_h * h
    du = [dv, dm, dn, dh]
    return du

# Now simulate
Iapp = 20
t_span = [0, t_max]
u_0 = [0, 0, 0, 0]
foo = solve_ivp(fun, t_span, u_0, args=(Iapp,), t_eval=t_eval)
t = foo.t
v = foo.y[0]

# Add nans to the start for initial traces
nans = np.empty(nt_ax)
nans.fill(np.nan)
t = np.concatenate((nans, t))
v = np.concatenate((nans, v))

def make_cmap():
# OK, now for the tricky bit. We want most of the trace to be opaque (maybe 80%), and then for opacity to linearly decay away for the last 20% of the entries

    # main trace colour
    col = (1.0, 0.0, 0.0)

    colours = [col + (0,), col + (1,), col + (1,)]
    nodes = [0.0, 0.2, 1.0]

    cm = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colours)))
    return cm

# Get colourmap that deals with transparency
cm = make_cmap()

# Create plot with segments

pts = np.array([t[0:nt_ax], v[0:nt_ax]]).T.reshape(-1, 1, 2)
segments = np.concatenate([pts[:-1], pts[1:]], axis=1)

fig, ax = plt.subplots()

norm = plt.Normalize(0, t_max_ax)
lc = LineCollection(segments, cmap=cm, norm=norm)
lc.set_array(t[0:nt_ax])
line = ax.add_collection(lc)
ax.set_xlim(0, t_max_ax)
ax.set_ylim(-20, 100)

def animate(i):
    t_data = np.mod(t[i:(nt_ax+1)], t_max_ax)
    pts = np.array([t_data, v[i:(nt_ax + i)]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc.set_segments(segments)
    return line

ani = FuncAnimation(fig, animate, nt, 1000/fps, blit=False)



plt.draw()
plt.show()


