import numpy as np
import matplotlib.pyplot as plt

print("Creating graphene property plots...")

# 1. Honeycomb Lattice
print("1. Creating honeycomb lattice plot...")
a = 1.0  # scaled for visualization
a1 = a * np.array([np.sqrt(3)/2, 1/2])
a2 = a * np.array([np.sqrt(3)/2, -1/2])

n1_vals = np.arange(-3, 4)
n2_vals = np.arange(-3, 4)

atoms_A = []
atoms_B = []

for n1 in n1_vals:
    for n2 in n2_vals:
        R = n1 * a1 + n2 * a2
        atoms_A.append(R)
        atoms_B.append(R + (a/np.sqrt(3)) * np.array([0, 1]))

atoms_A = np.array(atoms_A)
atoms_B = np.array(atoms_B)

fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(atoms_A[:, 0], atoms_A[:, 1], c='red', s=100, label='Sublattice A', zorder=2)
ax.scatter(atoms_B[:, 0], atoms_B[:, 1], c='blue', s=100, label='Sublattice B', zorder=2)

corners = np.array([[0, 0], a1, a1+a2, a2, [0, 0]])
ax.plot(corners[:, 0], corners[:, 1], 'k--', linewidth=2, alpha=0.7, label='Unit Cell')

ax.set_xlabel('x (scaled units)')
ax.set_ylabel('y (scaled units)')
ax.set_title('Honeycomb Lattice of Single-Layer Graphene')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig('honeycomb_lattice.png', dpi=300, bbox_inches='tight')
plt.close()  # Close figure to free memory
print("   Honeycomb lattice plot saved!")

# 2. Band Structure
print("2. Creating band structure plot...")
k = np.linspace(-0.5, 0.5, 1000)
E_plus = np.abs(k)  # Conduction band
E_minus = -np.abs(k)  # Valence band

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(k, E_plus, 'b-', linewidth=2, label='Conduction Band')
ax.plot(k, E_minus, 'r-', linewidth=2, label='Valence Band')

ax.set_xlabel('Wave Vector k (scaled units)')
ax.set_ylabel('Energy E (scaled units)')
ax.set_title('Band Structure of Graphene - Linear Dispersion at Dirac Point')
ax.grid(True, alpha=0.3)
ax.legend()
ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
ax.axvline(x=0, color='k', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('band_structure.png', dpi=300, bbox_inches='tight')
plt.close()
print("   Band structure plot saved!")

# 3. Density of States
print("3. Creating density of states plot...")
energies = np.linspace(-1.0, 1.0, 1000)
dos = np.abs(energies)

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(energies, dos, 'g-', linewidth=2)
ax.set_xlabel('Energy (scaled units)')
ax.set_ylabel('Density of States (a.u.)')
ax.set_title('Density of States of Single-Layer Graphene')
ax.grid(True, alpha=0.3)
ax.axvline(x=0, color='k', linestyle='--', alpha=0.5, label='Dirac Point')
ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('density_of_states.png', dpi=300, bbox_inches='tight')
plt.close()
print("   Density of states plot saved!")

# 4. Quantum Hall Effect
print("4. Creating quantum Hall effect plot...")
n_range = np.arange(-10, 11, 1)
conductance = n_range + 0.5

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(n_range, conductance, 'bo-', markersize=6, linewidth=1.5)

ax.set_xlabel('Landau Level Index (n)')
ax.set_ylabel('Hall Conductance (in units of e²/h)')
ax.set_title('Quantum Hall Effect in Graphene (Half-Integer Quantization)')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('quantum_hall_effect.png', dpi=300, bbox_inches='tight')
plt.close()
print("   Quantum Hall effect plot saved!")

# 5. Optical Absorption
print("5. Creating optical absorption plot...")
wavelengths = np.linspace(300, 2000, 1000)
absorption_percent = np.full_like(wavelengths, np.pi * 1/137 * 100)

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(wavelengths, absorption_percent, 'purple', linewidth=2)
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Absorption (%)')
ax.set_title('Universal Optical Absorption of Single-Layer Graphene (~2.3%)')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('optical_absorption.png', dpi=300, bbox_inches='tight')
plt.close()
print("   Optical absorption plot saved!")

print("\nAll plots generated successfully!")
print("Generated files:")
print("- honeycomb_lattice.png: Visualization of the honeycomb lattice structure")
print("- band_structure.png: Linear dispersion relation at Dirac points")
print("- density_of_states.png: Vanishing DOS at the Dirac point")
print("- quantum_hall_effect.png: Half-integer quantum Hall effect")
print("- optical_absorption.png: Universal ~2.3% optical absorption")

print("\nKey findings demonstrated:")
print("1. Honeycomb lattice structure with two sublattices (A and B)")
print("2. Linear dispersion relation at Dirac points (E = ±v_F |k|)")
print("3. Vanishing density of states at the Dirac point")
print("4. Half-integer quantum Hall effect")
print("5. Universal optical absorption of ~2.3% per layer")