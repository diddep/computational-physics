import numpy as np
import matplotlib.pyplot as plt

# Define the simple roots of the sl(2) Lie algebra
alpha1 = np.array([1, -1])
alpha2 = np.array([-2, 1])

# Define the simple coroots of the sl(2) Lie algebra
alpha1_star = np.array([1, 0])
alpha2_star = np.array([0, 1])

# Define the Weyl group of sl(2), which consists of all permutations and sign flips of the simple roots
weyl_group = [np.array([1, 1]), np.array([1, -1]), np.array([-1, 1]), np.array([-1, -1])]

# Compute the set of all roots by applying the Weyl group to the simple roots
roots = [alpha1, alpha2]
for alpha in roots:
    for s in weyl_group:
        new_root = alpha * s
        if not np.any([np.allclose(new_root, alpha_prime) for alpha_prime in roots]):
            roots.append(new_root)

# Plot the root lattice
fig = plt.figure()
ax = fig.add_subplot(111)

# Plot the root lattice points
for alpha in roots:
    ax.scatter(alpha[0], alpha[1])

# Plot the lines connecting the root lattice points
for alpha1 in roots:
    for alpha2 in roots:
        ax.plot([alpha1[0], alpha2[0]], [alpha1[1], alpha2[1]], 'k-')

# Set the axes labels
ax.set_xlabel('$\\alpha_1$')
ax.set_ylabel('$\\alpha_2$')

# Show the plot
plt.savefig('sl(2).png')
plt.show()

