import numpy as np
import matplotlib.pyplot as plt

def plot_honeycomb_lattice():
    """Simple visualization of graphene honeycomb lattice"""
    # Lattice constant
    a = 1.0  # scaled for visualization
    
    # Basis vectors
    a1 = a * np.array([np.sqrt(3)/2, 1/2])
    a2 = a * np.array([np.sqrt(3)/2, -1/2])
    
    # Generate lattice points
    n1_vals = np.arange(-5, 6)
    n2_vals = np.arange(-5, 6)
    
    atoms_A = []
    atoms_B = []
    
    for n1 in n1_vals:
        for n2 in n2_vals:
            R = n1 * a1 + n2 * a2
            atoms_A.append(R)
            atoms_B.append(R + (a/np.sqrt(3)) * np.array([0, 1]))
    
    atoms_A = np.array(atoms_A)
    atoms_B = np.array(atoms_B)
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot A and B sublattices
    ax.scatter(atoms_A[:, 0], atoms_A[:, 1], c='red', s=100, label='Sublattice A', zorder=2)
    ax.scatter(atoms_B[:, 0], atoms_B[:, 1], c='blue', s=100, label='Sublattice B', zorder=2)
    
    # Draw unit cell
    corners = np.array([[0, 0], a1, a1+a2, a2, [0, 0]])
    ax.plot(corners[:, 0], corners[:, 1], 'k--', linewidth=2, alpha=0.7, label='Unit Cell')
    
    ax.set_xlabel('x (scaled units)')
    ax.set_ylabel('y (scaled units)')
    ax.set_title('Honeycomb Lattice of Single-Layer Graphene')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig('honeycomb_lattice_simple.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Honeycomb lattice plot saved as honeycomb_lattice_simple.png")

def plot_graphene_band_structure():
    """Simple plot of graphene band structure with linear dispersion at Dirac points"""
    # Generate k-points along high symmetry path
    # Simplified path: just show the linear crossing at the Dirac point
    
    k = np.linspace(-0.5, 0.5, 1000)
    
    # Linear dispersion: E = ±v_F * |k| (simplified)
    E_plus = np.abs(k)  # Conduction band
    E_minus = -np.abs(k)  # Valence band
    
    fig, ax = plt.subplots(figsize=(10, 6))
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
    plt.savefig('band_structure_simple.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Band structure plot saved as band_structure_simple.png")

def plot_density_of_states():
    """Plot density of states for graphene (linear behavior at Dirac point)"""
    # Energy values near the Dirac point
    energies = np.linspace(-1.0, 1.0, 1000)
    
    # Density of states in graphene is proportional to |E| near the Dirac point
    dos = np.abs(energies)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(energies, dos, 'g-', linewidth=2)
    ax.set_xlabel('Energy (scaled units)')
    ax.set_ylabel('Density of States (a.u.)')
    ax.set_title('Density of States of Single-Layer Graphene')
    ax.grid(True, alpha=0.3)
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.5, label='Dirac Point')
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('density_of_states_simple.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Density of states plot saved as density_of_states_simple.png")

def plot_quantum_hall_effect():
    """Plot the half-integer quantum Hall effect in graphene"""
    # Landau level index
    n_range = np.arange(-10, 11, 1)  # Include negative indices for hole states
    
    # Quantum Hall conductance in graphene: σ_xy = (4e²/h)(n + 1/2)
    # Simplified for visualization
    conductance = n_range + 0.5
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(n_range, conductance, 'bo-', markersize=8, linewidth=2)
    
    ax.set_xlabel('Landau Level Index (n)')
    ax.set_ylabel('Hall Conductance (in units of e²/h)')
    ax.set_title('Quantum Hall Effect in Graphene (Half-Integer Quantization)')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('quantum_hall_simple.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Quantum Hall effect plot saved as quantum_hall_simple.png")

def plot_optical_absorption():
    """Plot universal optical absorption of graphene"""
    # Wavelength range
    wavelengths = np.linspace(300, 2000, 1000)  # nm
    
    # Graphene absorbs πα ≈ 2.3% of incident light per layer (universal value)
    absorption_percent = np.full_like(wavelengths, np.pi * 1/137 * 100)  # π times fine structure constant
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wavelengths, absorption_percent, 'purple', linewidth=2)
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Absorption (%)')
    ax.set_title('Universal Optical Absorption of Single-Layer Graphene (~2.3%)')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('optical_absorption_simple.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Optical absorption plot saved as optical_absorption_simple.png")

if __name__ == "__main__":
    print("Generating simple graphene property plots...")
    
    plot_honeycomb_lattice()
    plot_graphene_band_structure()
    plot_density_of_states()
    plot_quantum_hall_effect()
    plot_optical_absorption()
    
    print("All plots generated successfully!")
    print("\nKey findings demonstrated:")
    print("1. Honeycomb lattice structure with two sublattices (A and B)")
    print("2. Linear dispersion relation at Dirac points (E = ±v_F |k|)")
    print("3. Vanishing density of states at the Dirac point")
    print("4. Half-integer quantum Hall effect")
    print("5. Universal optical absorption of ~2.3% per layer")