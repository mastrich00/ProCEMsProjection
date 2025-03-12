import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Configure plot style
plt.style.use('seaborn-whitegrid')
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'serif',
    'axes.labelpad': 12
})

# Number of facets
n = 6

# Base radius and height
base_radius = 3
base_height = 3

# Seed for reproducibility
np.random.seed(0)

# Generate the rays
generators = []
for i in range(n):
    base_angle = 2 * np.pi * i / n
    # Introduce a small random variation in angle
    theta = base_angle + np.random.uniform(-0.1, 0.1)
    # Introduce slight variations in radius and height
    r = base_radius * (1 + np.random.uniform(-0.1, 0.1))
    h = base_height * (1 + np.random.uniform(-0.1, 0.1))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    generators.append([x, y, h])
generators = np.array(generators)

# Rotate the cone to tilt it 45 degrees
alpha = np.deg2rad(45)
R_y = np.array([
    [np.cos(alpha), 0, np.sin(alpha)],
    [0, 1, 0],
    [-np.sin(alpha), 0, np.cos(alpha)]
])
generators = np.dot(generators, R_y.T)

# Translate the cone to fit within [0,10] x [0,10] x [0,10]
translation = np.array([5, 5, 5])  # Move the vertex to the center of the [0,10] cube
generators += translation

# Create the plot
fig = plt.figure(figsize=(7, 5))
ax = fig.add_subplot(111, projection='3d', proj_type='ortho')
ax.set_title(r'$\mathbf{3D\ Polyhedral\ Cone}$', fontsize=14, pad=20)

# Plot the vertex
ax.scatter(0, 0,0, color='#e74c3c', s=80, label='Vertex', zorder=5)
# ax.text(translation[0], translation[1], translation[2] - 0.5, 'Vertex', color='#e74c3c', ha='center')

# Plot the extreme rays
for vec in generators:
    ax.quiver(0, 0, 0, vec[0]*.7, vec[1]*.7, vec[2]*.7,
              color='#2c3e50', lw=2, arrow_length_ratio=0.1)

# Plot the facets
t = np.linspace(0, 1, 40)
s = np.linspace(0, 1, 40)
T, S = np.meshgrid(t, s)

for i in range(n):
    v1 = generators[i]
    v2 = generators[(i + 1) % n]
    X = T * ((1 - S) * v1[0] + S * v2[0])
    Y = T * ((1 - S) * v1[1] + S * v2[1])
    Z = T * ((1 - S) * v1[2] + S * v2[2])
    ax.plot_surface(X, Y, Z, color='#3498db', alpha=0.25, edgecolor='none', antialiased=True)

# Configure axes and view
ax.set(xlabel='X Axis', ylabel='Y Axis', zlabel='Z Axis',
       xlim=[0, 10], ylim=[0, 10], zlim=[0, 10])
ax.view_init(elev=24, azim=-135)
ax.set_box_aspect([1, 1, 1])

plt.tight_layout()
plt.savefig('polyhedral_cone.pdf', bbox_inches='tight', dpi=300)
plt.show()