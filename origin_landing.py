import numpy as np
import matplotlib.pyplot as plt

# Load and format data from files
def load_data(vw, dt):
    file = f"t_vw_{vw:.4f}_t_dt_{dt:.3f}.dat"
    t, x, y, vx, vy = np.loadtxt(file, unpack=True)
    return t, x, y, vx, vy

# Plotting function
def plot_decor(indep, dep, title, ylim, xlim, ):
    plt.plot(indep, dep, label = fr'Trial {i+1}: wind $\vec v = {vw}$ m/s' , linestyle = '--')
    plt.title(title)
    plt.ylabel('Y (m)')
    plt.xlabel('X (m)')
    plt.ylim(*ylim)
    plt.xlim(*xlim)
    plt.axhline(0, color='k')
    plt.axvline(0, color='k')
    plt.grid(True)
    plt.legend()
    
# Parameters
v0 = 6.9  
th0 = 34
h0 = 7.7
dt = 0.005
wind_velocities = [-1.0, -0.5, -0.7, -0.68, -0.6775]

# Set figure dimension and title
plt.figure(figsize=(8, 12))
plt.suptitle('Problem 1 - Part B', bbox={"facecolor":"white", "alpha":0.5, "pad":5})

# Loop through and plot
for i, vw in enumerate(wind_velocities):
    t, x, y, vx, vy = load_data(vw, dt)

    plt.subplot(2, 1, 1) # First subplot for trajectories
    plot_decor(x, y,
             fr'Falling object trajectory with $\theta_0 = {th0}^\circ, h_0 = {h0}$ m, $v_0 = {v0}$ m/s',
             (-0.5, 9),
             (-2, 3))
    plt.scatter(x[0], y[0], s = 20, marker = 'o', color = 'black', zorder = 5)

    plt.subplot(2, 1, 2) # Second subplot for (0,0) zoom-in
    plot_decor(x, y,
             fr'Zooming in on (0,0) ...',
             (-0.25, 0.25),
             (-0.25, 0.25))

plt.figtext(0.5, 0.01, fr"wind $\vec v = {vw}$ m/s seems to be optimal", ha="center", fontsize=14, bbox={"facecolor":"orange", "alpha":0.5, "pad":5})

# Box the subplots
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig("TimothyIvanov_FE_P1-partB.pdf")

plt.show()