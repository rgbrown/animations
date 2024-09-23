import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

fps = 30

# Function to rotate and x limits
fun = lambda x: x**3
fun = lambda x: np.sqrt(4 - x**2)
xlim = [-2, 2]

# Construct the plots
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

nx, ntheta = 51, 51 
x = np.linspace(xlim[0], xlim[1], nx)
theta = np.linspace(0, 2*np.pi, ntheta)
X, T = np.meshgrid(x, theta)

Y = fun(X) * np.cos(T)
Z = fun(X) * np.sin(T)

surf = ax.plot_surface(X, Y, Z, color=('blue', 0.3))
lp, = ax.plot(x, fun(x), 0,'r')
ax.plot(x, fun(x), 0, 'b')
ax.plot(xlim, [0, 0], 0, 'k')
plt.xlabel('x')
plt.ylabel('y')

def animate(k):
    th = theta[np.mod(k, ntheta)]
    lp.set_data_3d(x, fun(x)*np.cos(th), fun(x)*np.sin(th))
    return lp

ani = FuncAnimation(fig, animate, 1000, interval = 1000/fps, blit=False)
plt.show()




