# animations

## Overview

This repository is an educational collection of small animations for teaching ideas from applied mathematics, physics, dynamical systems, numerical methods, and linear algebra. Most of the examples are written as self-contained Python scripts, with two older MATLAB scripts for matrix factorization demos.

The code is designed to make the underlying concept visible rather than to provide a reusable package. Each script solves or constructs a mathematical model, renders the result with simple plotting tools, and in several cases saves the animation to an `.mp4` file for classroom use.

## Repository contents

### Python animations

- `fourierseries.py` — Fourier series approximation visualization
- `hodgkinhuxley.py` — Hodgkin-Huxley neuron membrane dynamics
- `lorenz.py` — Lorenz attractor with many chaotic trajectories
- `pendulum.py` — nonlinear pendulum with velocity and acceleration vectors
- `threebody.py` — planar gravitational three-body simulation
- `volrot.py` — volume of revolution visualization

### MATLAB animations

- `plu.m` — PLU matrix decomposition walkthrough
- `qranim.m` — QR decomposition via Givens rotations

## Requirements and installation

### Python requirements

The Python scripts use standard scientific Python tools:

- Python 3
- `numpy`
- `matplotlib`
- `scipy`

Install the dependencies with:

```bash
python -m pip install numpy matplotlib scipy
```

### MATLAB requirements

The MATLAB scripts are plain `.m` files and rely on standard MATLAB matrix and display functionality. No additional toolboxes are documented in the repository.

### Optional video writer support

Several Python scripts save `.mp4` output. Depending on your Matplotlib setup, saving may require an available video writer such as `ffmpeg`.

## Running instructions

### Python

Run any script directly from the repository root:

```bash
python fourierseries.py
python hodgkinhuxley.py
python lorenz.py
python pendulum.py
python threebody.py
python volrot.py
```

What to expect:

- Most scripts open an interactive Matplotlib window.
- `fourierseries.py`, `lorenz.py`, `pendulum.py`, and `threebody.py` save an `.mp4` file in the current working directory.
- `hodgkinhuxley.py` displays the animation but does not currently save a movie file.
- `volrot.py` displays the rotating 3D visualization but does not currently save a movie file.

### MATLAB

From MATLAB, change into the repository directory and run:

```matlab
plu
qranim
```

The MATLAB scripts are intended as live visual or step-by-step demonstrations in the command window and figure environment.

## Python animation details

### `fourierseries.py`

**Concept being taught**

- Fourier series approximation of a periodic function
- Convergence of partial sums
- Relationship between a target waveform and its harmonic components

**What the animation visualizes**

The script builds successive partial sums of a Fourier series and animates how the approximation develops term by term against the target periodic function. The true function is plotted once, while the animated curve shows the approximation improving as more sine modes are added.

**Key parameters**

- `nx = 1001` — spatial resolution of the plotted function
- `fps = 10` — animation playback and saved video frame rate
- `L = 1` — half-period scaling used inside the Fourier modes
- `n_terms = 60` — number of partial sums shown
- `x = np.linspace(-4, 4, nx)` — plotted domain

**Teaching value**

This is useful for showing students how global approximations are assembled from simple oscillatory pieces. It also gives a concrete visual explanation of why early partial sums capture broad structure first and finer corrections later.

**Notable features**

- Dark background for projection or lecture use
- Precomputes all partial sums before animating
- Overlays the animated approximation on the target function
- Minimal axes for a clean visual presentation

**Output details**

- Opens a Matplotlib figure
- Saves `fourier.mp4`
- Frame count is tied to `n_terms`

### `hodgkinhuxley.py`

**Concept being taught**

- Hodgkin-Huxley membrane dynamics
- Excitable systems and action potential generation
- Coupling between voltage and ion-channel gating variables

**What the animation visualizes**

The script solves the Hodgkin-Huxley ODE system under constant applied current and animates the membrane voltage trace in time. The displayed trace uses a fading tail effect so recent motion remains prominent while older segments become transparent.

**Key parameters**

- `dt = 0.05` ms — simulation step used for evaluation points
- `t_max = 300` — total simulation time
- `t_max_ax = 50` — width of the visible time window
- `fps = 30` — animation speed
- `Iapp = 20` — constant applied current
- Conductances: `gNa = 120`, `gK = 36`, `gL = 0.3`
- Reversal potentials: `vNa = 115`, `vK = -12`, `vL = 10.6`

**Teaching value**

This makes electrophysiology easier to teach because the animation shows the membrane voltage as a living process rather than a static plot. The moving time window is particularly helpful for discussing repetitive spiking and the effect of sustained input current.

**Notable features**

- Uses `scipy.integrate.solve_ivp` for time integration
- Includes a custom transparent colormap for the fading trace
- Displays only a moving recent-history window instead of the full signal at once
- Encodes the standard four Hodgkin-Huxley state variables internally (`v`, `m`, `n`, `h`)

**Output details**

- Opens an animated Matplotlib plot of membrane voltage versus time
- Does not currently call `save`, so no video file is written by default
- Uses y-limits from `-20` to `100` for the voltage display

### `lorenz.py`

**Concept being taught**

- Chaotic dynamics in the Lorenz system
- Sensitivity to initial conditions
- Strange attractors and trajectory geometry in phase space

**What the animation visualizes**

The script samples many initial conditions, integrates the Lorenz equations for each one, and animates short moving trails in 3D. The camera angle also changes over time, which makes the butterfly-shaped attractor structure easier to appreciate from different perspectives.

**Key parameters**

- `fps = 30` — animation and video frame rate
- `pts_per_frame = 5` — temporal resolution per rendered frame
- `sim_time = 40` seconds — total simulated duration
- `speed_factor = 0.5` — scales the Lorenz dynamics
- `trail_time = 1` second — visible trail length
- Lorenz parameters: `rho = 28`, `sigma = 10`, `beta = 8/3`
- `n_traces = 30` — number of simultaneous trajectories

**Teaching value**

This is effective for introducing deterministic chaos because students can see multiple nearby trajectories evolve into complex motion on the same attractor. The rotating viewpoint reinforces that chaotic structure is geometric, not just a complicated time series.

**Notable features**

- Uses many simultaneous trajectories rather than a single orbit
- Adds finite-length trails instead of plotting full paths
- Continuously rotates the 3D camera during playback
- Uses randomized initial conditions in a large box of phase space

**Output details**

- Opens a 3D Matplotlib animation
- Saves `lorenz.mp4`
- Uses axis limits `x, y in [-20, 20]` and `z in [0, 40]`

### `pendulum.py`

**Concept being taught**

- Nonlinear pendulum dynamics
- Relationship between position, velocity, and acceleration
- Polar-to-Cartesian interpretation of mechanical motion

**What the animation visualizes**

The script solves the nonlinear pendulum equation and animates the bob position together with velocity and acceleration vectors attached to the mass. Labels move with the vector tips so the diagram behaves like a dynamic mechanics sketch.

**Key parameters**

- `L = 4` — pendulum length
- `g = 9.8` — gravitational acceleration
- `arrow_scale = 0.3` — scales vector lengths for visibility
- `fps = 30` — animation and video frame rate
- `T = 20` — total animation duration
- `theta0 = 45*np.pi/180` — initial angular displacement
- `omega0 = 0` — initial angular velocity

**Teaching value**

This is especially useful in mechanics teaching because it connects the differential equation, the geometric motion, and the kinematic vectors in one visual. It can support discussion of tangential versus centripetal contributions and of how acceleration differs from velocity in direction and magnitude.

**Notable features**

- Solves the full nonlinear equation, not the small-angle approximation
- Computes Cartesian position, velocity, and acceleration from angular variables
- Draws labeled velocity and acceleration arrows in real time
- Uses an uncluttered full-figure layout with hidden axes

**Output details**

- Opens a Matplotlib animation window
- Saves `pendulum.mp4`
- Uses a fixed 2D frame with equal aspect ratio

### `threebody.py`

**Concept being taught**

- Gravitational N-body interactions in a simplified three-body setting
- Coupled nonlinear ODE systems
- Orbit complexity and sensitivity in classical mechanics

**What the animation visualizes**

The script integrates a planar three-body gravitational system with three equal masses and animates the particle positions together with recent trajectory trails. Each body is shown in a different color so the evolving orbital interaction is easy to follow.

**Key parameters**

- `fps = 30` — animation and video frame rate
- `pts_per_frame = 10` — temporal resolution per rendered frame
- `sim_time = 180` seconds — total simulated duration
- `trail_time = 2` seconds — visible history length
- `m1 = m2 = m3 = 1` — equal masses
- `G = 0.1` — gravitational constant used in the model
- Initial positions: `x1_0`, `x2_0`, `x3_0`
- Initial velocities are currently all set to zero-valued vectors

**Teaching value**

This provides a concrete demonstration of why even small gravitational systems can be difficult to analyze exactly. It works well for showing students the difference between tractable two-body motion and richer multi-body behavior.

**Notable features**

- Uses colored trailing paths for all three bodies
- Integrates the full coupled state of positions and velocities
- Pads the start of the animation with `NaN` values so trails grow in naturally
- Keeps the view fixed and clean for presentation use

**Output details**

- Opens a 2D Matplotlib animation
- Saves `threebody.mp4`
- Uses plot limits `x in [-1.2, 1.2]` and `y in [-0.5, 1.3]`

### `volrot.py`

**Concept being taught**

- Volume of revolution
- Building a 3D surface by rotating a planar curve
- Relationship between a generating function and the resulting solid

**What the animation visualizes**

The script constructs a surface of revolution by rotating a profile curve around the x-axis, then animates a highlighted generating curve as it sweeps through angle. In the current file, the active profile is the semicircular curve `sqrt(4 - x^2)` over `[-2, 2]`.

**Key parameters**

- `fps = 30` — animation playback rate
- `fun` — generating function being revolved
- `xlim = [-2, 2]` — x-domain used for the profile
- `nx = 51` — x-resolution of the mesh
- `ntheta = 51` — angular resolution of the surface

**Teaching value**

This is useful for calculus or geometry classes because it gives a direct visual connection between a 2D generating graph and a 3D solid. It can support lessons on washers, shells, parameterizations, and geometric intuition for integral formulas.

**Notable features**

- Uses a 3D Matplotlib surface plot
- Highlights a moving copy of the generating curve during rotation
- Leaves an alternative cubic profile in the file as an easy starting point for experimentation
- Simple structure makes it easy to adapt to other functions

**Output details**

- Opens an interactive 3D Matplotlib figure
- Does not currently save a movie file
- Displays both the revolved surface and the profile curve

## MATLAB animation details

### `plu.m`

**Concept being taught**

- PLU factorization with pivoting
- Elimination and row permutations in linear algebra
- Relationship between an original matrix and its triangular factors

**What the animation visualizes**

This script performs a PLU-style elimination process on a random square matrix while tracking the permutation vector. It is less of a graphical animation and more of a compact computational walkthrough that reveals how the lower-triangular multipliers and upper-triangular matrix are produced.

**Key parameters**

- `n = 4` — matrix size
- `A = rand(n)` — random dense test matrix
- `p` — permutation index vector
- `s` — row scaling values used for pivot selection

**Teaching value**

This is a useful teaching script for explaining partial or scaled pivoting and for showing that factorization is an algorithmic process, not just a symbolic theorem. Because the matrix is small, it can be stepped through carefully in class.

**Notable features**

- Stores multipliers directly in the lower part of `A`
- Tracks row permutations separately with `p`
- Reconstructs `L` and `U` at the end
- Prints `L*U` and the original matrix for comparison

**Output details**

- Runs in MATLAB as a script
- Produces command-window output for the factorization result
- Does not create a saved movie file

### `qranim.m`

**Concept being taught**

- QR factorization
- Givens rotations for zeroing subdiagonal entries
- Orthogonal transformations in numerical linear algebra

**What the animation visualizes**

The script repeatedly applies Givens rotations to transform a random matrix toward upper-triangular form. After each rotation, it clears the command window and prints the evolving `Q` and `R` factors, creating a step-by-step animated effect.

**Key parameters**

- `n = 15` — matrix size
- `A = randn(n)` — random initial matrix
- `Q = eye(n)` — accumulated orthogonal factor
- `pause(0.15)` — delay between displayed updates

**Teaching value**

This is well suited to demonstrating how QR factorization is assembled from simple plane rotations. It helps students see that orthogonal transformations act locally on row pairs while building a global decomposition.

**Notable features**

- Uses Givens rotations explicitly rather than hiding them inside a library call
- Prints intermediate `Q` and `R` after every elimination step
- Clears the command window between frames for a live-animation effect
- Forces tiny roundoff-level values to zero for cleaner display

**Output details**

- Runs as a MATLAB script
- Produces an animated command-window progression of `Q` and `R`
- Does not save a movie file

## Customization guide

Most scripts are intentionally compact, so the easiest way to customize them is to edit top-level parameters near the start of each file.

### Common Python customizations

- Change `fps` to control smoothness and output file size
- Change simulation duration variables such as `T`, `sim_time`, or `t_max`
- Adjust resolution variables like `nx`, `ntheta`, or `pts_per_frame`
- Modify initial conditions to generate different trajectories or responses
- Tune axis limits if you want tighter framing for lecture slides or exported videos

### Concept-specific ideas

- **Fourier series**: increase `n_terms`, change the target function, or narrow the x-range to focus on local convergence behavior
- **Hodgkin-Huxley**: vary `Iapp` to compare subthreshold response, spiking onset, and sustained firing
- **Lorenz attractor**: change `rho`, `sigma`, `beta`, or `n_traces` to explore different chaotic or non-chaotic regimes
- **Pendulum**: change `theta0`, `omega0`, `L`, or `g` to compare amplitude effects and time scales
- **Three-body problem**: change masses, gravitational constant, or initial conditions to create qualitatively different orbits
- **Volume of revolution**: replace `fun` with another profile and update `xlim` accordingly
- **PLU / QR**: change the matrix size `n` to balance clarity against computational complexity in class

## Use cases

This repository can be useful for:

- classroom demonstrations in mathematics, physics, neuroscience, and numerical analysis
- lecture slides or recorded teaching videos
- exploratory coding exercises for students learning scientific Python or MATLAB
- quick prototypes for visual explanations of dynamical systems
- inspiration for building more polished educational simulations

## General notes

- These scripts are self-contained examples rather than a packaged library.
- Output files are written to the current working directory when `save(...)` is used.
- Some animations may be computationally heavier at higher frame rates or longer simulation times.
- If `.mp4` saving fails, check that your Matplotlib animation backend can access a video writer.
- The code favors readability and directness, which makes it well suited for adaptation in teaching settings.
