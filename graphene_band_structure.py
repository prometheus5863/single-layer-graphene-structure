import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set up plotting parameters
rcParams['font.size'] = 12
rcParams['figure.figsize'] = (10, 8)

def graphene_lattice_vectors():
    """
    Define the lattice vectors for graphene
    """
    # Lattice constant for graphene (in Angstrom)
    a = 2.46  # graphene lattice constant
    
    # Real space lattice vectors
    a1 = a * np.array([np.sqrt(3)/2, 1/2, 0])
    a2 = a * np.array([np.sqrt(3)/2, -1/2, 0])
    
    # Reciprocal lattice vectors
    b1 = (2*np.pi/a) * np.array([1/np.sqrt(3), 1, 0])
    b2 = (2*np.pi/a) * np.array([1/np.sqrt(3), -1, 0])
    
    return a1, a2, b1, b2

def k_path_graphene(n_points=100):
    """
    Generate k-points along high symmetry path in graphene Brillouin zone
    Path: Γ -> K -> M -> Γ
    """
    # High symmetry points in reciprocal space (in units of 2π/a)
    gamma = np.array([0.0, 0.0])  # Γ point
    k_point = np.array([4.0/(3*np.sqrt(3)), 0.0])  # K point  
    m_point = np.array([np.pi/(3*np.sqrt(3)), np.pi/3])  # M point
    
    # Interpolate between points
    k_path = []
    labels = []
    
    # Γ to K
    for i in range(n_points//3):
        t = i / (n_points//3 - 1) if n_points//3 > 1 else 0
        k = gamma + t * (k_point - gamma)
        k_path.append(k)
    labels.extend(['Γ'] + ['']*(n_points//3-1))
    
    # K to M
    for i in range(n_points//3):
        t = i / (n_points//3 - 1) if n_points//3 > 1 else 0
        k = k_point + t * (m_point - k_point)
        k_path.append(k)
    labels.extend(['K'] + ['']*(n_points//3-1))
    
    # M to Γ
    for i in range(n_points//3):
        t = i / (n_points//3 - 1) if n_points//3 > 1 else 0
        k = m_point + t * (gamma - m_point)
        k_path.append(k)
    labels.extend(['M'] + ['']*(n_points//3-1))
    
    return np.array(k_path), labels

def graphene_hamiltonian(k, hopping=2.8):  # hopping parameter in eV
    """
    Calculate the 2x2 Hamiltonian for graphene at wavevector k
    Based on tight-binding model with nearest neighbor hopping
    """
    # Nearest neighbor vectors in real space (in units of lattice constant)
    delta1 = np.array([0, 1])
    delta2 = np.array([np.sqrt(3)/2, -1/2])
    delta3 = np.array([-np.sqrt(3)/2, -1/2])
    
    # Phase factors
    phi_k = (
        np.exp(1j * np.dot(k, delta1)) + 
        np.exp(1j * np.dot(k, delta2)) + 
        np.exp(1j * np.dot(k, delta3))
    )
    
    # Hamiltonian matrix
    H = np.array([
        [0, phi_k],
        [np.conj(phi_k), 0]
    ])
    
    # Scale by hopping parameter
    H *= hopping
    
    return H

def calculate_band_structure():
    """
    Calculate the band structure of graphene along high symmetry path
    """
    k_points, labels = k_path_graphene()
    energies = []
    
    for k in k_points:
        H = graphene_hamiltonian(k)
        eigenvals = np.sort(np.linalg.eigvalsh(H))
        energies.append(eigenvals)
    
    energies = np.array(energies)
    
    return k_points, energies, labels

def plot_band_structure():
    """
    Plot the band structure of graphene
    """
    k_points, energies, labels = calculate_band_structure()
    
    # Calculate distances along the path
    distances = np.zeros(len(k_points))
    for i in range(1, len(k_points)):
        distances[i] = distances[i-1] + np.linalg.norm(k_points[i] - k_points[i-1])
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot bands
    ax.plot(distances, energies[:, 0], 'r-', linewidth=2, label='Valence Band')
    ax.plot(distances, energies[:, 1], 'b-', linewidth=2, label='Conduction Band')
    
    # Mark high symmetry points
    # Find positions of high symmetry points
    gamma_pos = distances[0]
    k_pos = distances[len(k_points)//3]
    m_pos = distances[2*len(k_points)//3]
    
    ax.axvline(x=gamma_pos, color='k', linestyle='--', alpha=0.5)
    ax.axvline(x=k_pos, color='k', linestyle='--', alpha=0.5)
    ax.axvline(x=m_pos, color='k', linestyle='--', alpha=0.5)
    
    # Labels for high symmetry points
    ax.text(gamma_pos, energies.min(), 'Γ', ha='center', va='top', fontsize=14)
    ax.text(k_pos, energies.min(), 'K', ha='center', va='top', fontsize=14)
    ax.text(m_pos, energies.min(), 'M', ha='center', va='top', fontsize=14)
    
    ax.set_xlabel('Wave Vector (k)')
    ax.set_ylabel('Energy (eV)')
    ax.set_title('Band Structure of Single-Layer Graphene')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('band_structure.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_fermi_surface():
    """
    Plot the Fermi surface of graphene around the Dirac points
    """
    # Create a grid of k-points around the K point
    kx_range = np.linspace(-0.5, 0.5, 200)
    ky_range = np.linspace(-0.5, 0.5, 200)
    K_x, K_y = np.meshgrid(kx_range, ky_range)
    
    # Shift to center around K point (approximately)
    k_K = 4.0/(3*np.sqrt(3))  # x-coordinate of K point
    K_x_shifted = K_x + k_K
    
    # Calculate energy at each k-point
    energies = np.zeros_like(K_x)
    for i in range(K_x.shape[0]):
        for j in range(K_x.shape[1]):
            k = np.array([K_x_shifted[i,j], K_y[i,j]])
            H = graphene_hamiltonian(k)
            eigenvals = np.sort(np.linalg.eigvalsh(H))
            # Take the difference between conduction and valence bands
            energies[i,j] = eigenvals[1] - eigenvals[0]
    
    fig, ax = plt.subplots(figsize=(8, 8))
    contour = ax.contour(K_x_shifted, K_y, energies, levels=20, cmap='viridis')
    ax.contour(K_x_shifted, K_y, energies, levels=[0.01], colors='red', linewidths=2)  # Fermi level
    ax.set_xlabel('k_x (2π/a)')
    ax.set_ylabel('k_y (2π/a)')
    ax.set_title('Fermi Surface of Single-Layer Graphene (Around K Point)')
    plt.colorbar(contour)
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('fermi_surface.png', dpi=300, bbox_inches='tight')
    plt.show()

def calculate_density_of_states():
    """
    Calculate the density of states for graphene (analytical approximation near Dirac point)
    """
    # Energy values near the Dirac point
    energies = np.linspace(-3.0, 3.0, 1000)
    
    # Linear dispersion relation near Dirac point: E = ±v_F * |k|
    # Density of states is proportional to |E| near the Dirac point
    dos = np.abs(energies) / (np.pi * (2.8)**2)  # Simplified form
    
    # Near-zero values need special treatment due to linear dispersion
    # In graphene, DOS goes as |E|, so it vanishes at the Dirac point
    return energies, dos

def plot_density_of_states():
    """
    Plot the density of states for graphene
    """
    energies, dos = calculate_density_of_states()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(energies, dos, 'g-', linewidth=2)
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Density of States (a.u.)')
    ax.set_title('Density of States of Single-Layer Graphene')
    ax.grid(True, alpha=0.3)
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.5, label='Dirac Point')
    ax.legend()
    plt.tight_layout()
    plt.savefig('density_of_states.png', dpi=300, bbox_inches='tight')
    plt.show()

def visualize_unit_cell():
    """
    Visualize the graphene unit cell and honeycomb lattice
    """
    # Lattice constant
    a = 2.46  # Angstrom
    
    # Basis vectors
    a1 = a * np.array([np.sqrt(3)/2, 1/2])
    a2 = a * np.array([np.sqrt(3)/2, -1/2])
    
    # Basis atoms in the unit cell
    basis_A = np.array([0, 0])
    basis_B = (a/np.sqrt(3)) * np.array([0, 1])  # Distance between layers in hexagonal structure
    
    # Generate lattice points
    n1_vals = np.arange(-3, 4)
    n2_vals = np.arange(-3, 4)
    
    atoms_A = []
    atoms_B = []
    
    for n1 in n1_vals:
        for n2 in n2_vals:
            R = n1 * a1 + n2 * a2
            atoms_A.append(R + basis_A)
            atoms_B.append(R + basis_B)
    
    atoms_A = np.array(atoms_A)
    atoms_B = np.array(atoms_B)
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot A and B sublattices
    ax.scatter(atoms_A[:, 0], atoms_A[:, 1], c='red', s=100, label='Sublattice A', zorder=2)
    ax.scatter(atoms_B[:, 0], atoms_B[:, 1], c='blue', s=100, label='Sublattice B', zorder=2)
    
    # Draw unit cell
    corners = np.array([[0, 0], a1, a1+a2, a2, [0, 0]])
    ax.plot(corners[:, 0], corners[:, 1], 'k--', linewidth=2, alpha=0.7, label='Unit Cell')
    
    # Draw bonds between nearest neighbors
    bond_distance = a / np.sqrt(3)  # distance between A and B atoms
    for atom_A in atoms_A[:25]:  # Just plot for a subset to avoid overcrowding
        for atom_B in atoms_B:
            dist = np.linalg.norm(atom_A - atom_B)
            if abs(dist - bond_distance) < 0.1:  # tolerance for bonding distance
                ax.plot([atom_A[0], atom_B[0]], [atom_A[1], atom_B[1]], 'k-', alpha=0.3, linewidth=1)
    
    ax.set_xlabel('x (Å)')
    ax.set_ylabel('y (Å)')
    ax.set_title('Honeycomb Lattice of Single-Layer Graphene')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig('honeycomb_lattice.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    print("Calculating and plotting graphene properties...")
    
    # Calculate and plot band structure
    print("Plotting band structure...")
    plot_band_structure()
    
    # Plot Fermi surface
    print("Plotting Fermi surface...")
    plot_fermi_surface()
    
    # Plot density of states
    print("Plotting density of states...")
    plot_density_of_states()
    
    # Visualize unit cell
    print("Plotting honeycomb lattice...")
    visualize_unit_cell()
    
    print("All plots generated successfully!")