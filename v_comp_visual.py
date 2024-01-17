import numpy as np
import matplotlib.pyplot as plt

# Load and format data from files
def load_data(vw, dt):
    file = f"vw_{vw:.1f}_dt_{dt:.3f}.dat"
    t, x, y, vx, vy = np.loadtxt(file, unpack=True)
    return t, x, y, vx, vy

# Plotting function
def plot_decor(indep, dep, vw, dt, title, ylab, xlab):
    plt.plot(indep, dep, label = f'vw={vw}, dt={dt}')
    plt.title(title)
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.grid(True)
    plt.legend()

# Parameters
v0 = 6.9  
th0 = 34
h01 = 7.7
h02 = 2.1
dt_values = [0.01, 0.005]
wind_velocities = [2.1, 0, -2.1]

# Set figure dimension and title
plt.figure(figsize=(8, 12))
plt.suptitle('Problem 1 - Part A', bbox={"facecolor":"white", "alpha":0.5, "pad":5})

# Loop through and plot
for vw in wind_velocities:
    for dt in dt_values:
        t, x, y, vx, vy = load_data(vw, dt)

        plt.subplot(3, 1, 1)  # First subplot for vx
        plot_decor(t, vx, vw, dt,
                   fr'Falling object $v_x$ component with $\theta_0 = {th0}^\circ, h_0 = {h01}$ m, $v_0 = {v0}$ m/s',
                   fr'$v_x(t)$',
                   'Time (s)')

        plt.subplot(3, 1, 2)  # Second subplot for vy
        plot_decor(t, vy, vw, dt,
                   fr'Falling object $v_y$ component with $\theta_0 = {th0}^\circ, h_0 = {h01}$ m, $v_0 = {v0}$ m/s',
                   fr'$v_y(t)$',
                   'Time (s)')

        plt.subplot(3, 1, 3)  # Third subplot for trajectories
        plot_decor(x, y, vw, dt,
                   fr'Falling object trajectory with $\theta_0 = {th0}^\circ, h_0 = {h02}$ m, $v_0 = {v0}$ m/s',
                   'Y (m)',
                   'X (m)')

# Extra details for third plot
plt.xlim(-1, 6)
plt.ylim(0, 4)
plt.xticks(np.arange(-1, 7, 1))
plt.axhline(0, color='k')
plt.scatter(x[0], y[0], s=30, marker='o', color='black', zorder=2)

# Box the subplots
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig("TimothyIvanov_FE_P1-partA.pdf")

plt.show()