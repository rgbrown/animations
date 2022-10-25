import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

plt.style.use("dark_background")

nx = 1001
fps = 10
x = np.linspace(-4, 4, nx)
L = 1
n_terms = 60
f_true = np.mod(x, 2)/(2*L)

F = np.zeros((n_terms, nx))
y = np.zeros(nx)
y.fill(0.5)
F[0] = y.copy()
for n in range(1, n_terms):
    y += -1/(n*np.pi)*np.sin(n*np.pi*x/L)
    F[n] = y.copy()

# Create figure
fig = plt.figure(figsize=(16, 9))
ax = plt.axes((0, 0, 1, 1), xlim=(-4, 4), ylim=(-1, 2))

nans = np.empty(nx)
nans.fill(np.nan)
yplot, = plt.plot(x, nans)
plt.plot(x, f_true)

plt.axis('off')

def animate(i):
    yplot.set_ydata(F[i])
    return yplot

ani = FuncAnimation(fig, animate, n_terms, interval=1000/fps, blit=False)

ani.save("fourier.mp4", fps=fps)


plt.show()