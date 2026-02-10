import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, e, Boltzmann as kB
import matplotlib.patches as patches

def calculate_fermi_velocity():
    """
    Calculate the Fermi velocity of graphene
    The linear dispersion relation near Dirac point gives rise to constant Fermi velocity
    """
    # Tight binding parameter (hopping integral) in eV
    t = 2.8  # eV
    
    # Lattice constant in meters
    a = 2.46e-10  # meters
    
    # Fermi velocity calculation
    # v_F = (3*a*t)/(2*hbar) where hbar is in J*s
    hbar_ev = 6.582119569e-16  # eV*s
    v_F = (3 * a * t) / (2 * hbar_ev)
    
    return v_F  # in m/s

def plot_linear_dispersion():
    """
    Plot the linear dispersion relation near the Dirac point
    """
    # Wave vector range (near Dirac point)
    k_range = np.linspace(-0.5, 0.5, 1000)  # in units of inverse lattice constant
    
    # Fermi velocity (in eV*angstrom)
    v_F = 1.0  # scaled for visualization
    
    # Linear dispersion: E = ±v_F * |k| 
    E_plus = v_F * np.abs(k_range)
    E_minus = -v_F * np.abs(k_range)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(k_range, E_plus, 'b-', linewidth=2, label='Conduction Band')
    ax.plot(k_range, E_minus, 'r-', linewidth=2, label='Valence Band')
    
    ax.set_xlabel('Wave Vector k (1/Å)')
    ax.set_ylabel('Energy E (eV)')
    ax.set_title('Linear Dispersion Relation of Graphene (Near Dirac Point)')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('linear_dispersion.png', dpi=300, bbox_inches='tight')
    plt.show()

def quantum_hall_conductance(B_field_range, n_carriers_range):
    """
    Calculate the quantum Hall conductance for graphene
    In graphene, the quantum Hall effect has unusual half-integer quantization
    σ_xy = (4e²/h)(n + 1/2) where n = 0, 1, 2, ...
    """
    # Conductance quantum
    G0 = 4 * e**2 / hbar  # Note: factor of 4 comes from spin and valley degeneracy
    
    # Calculate conductance for different Landau level indices
    n_values = n_carriers_range  # Landau level index
    conductance_values = G0 * (n_values + 0.5)
    
    return n_values, conductance_values

def plot_quantum_hall_effect():
    """
    Plot the quantum Hall effect in graphene
    """
    n_range = np.arange(-10, 11, 1)  # Include negative indices for hole states
    
    n_values, conductance_values = quantum_hall_conductance(None, n_range)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(n_values, conductance_values, 'bo-', markersize=8, linewidth=2)
    
    ax.set_xlabel('Landau Level Index (n)')
    ax.set_ylabel('Hall Conductance (e²/h)')
    ax.set_title('Quantum Hall Effect in Graphene (Half-Integer Quantization)')
    ax.grid(True, alpha=0.3)
    
    # Normalize by conductance quantum for better visualization
    ax_twin = ax.twinx()
    ax_twin.set_ylabel('Hall Conductance (units of 4e²/h)')
    ax_twin.set_ylim(ax.get_ylim()[0] / (4*e**2/hbar), ax.get_ylim()[1] / (4*e**2/hbar))
    
    plt.tight_layout()
    plt.savefig('quantum_hall_effect.png', dpi=300, bbox_inches='tight')
    plt.show()

def calculate_mobility_vs_temperature(T_range, scattering_mechanisms=True):
    """
    Calculate carrier mobility in graphene as a function of temperature
    Different scattering mechanisms dominate at different temperatures
    """
    T = T_range
    
    # Scattering mechanisms in graphene:
    # 1. Acoustic phonon scattering: μ ~ T^(-1.5) at high T
    # 2. Optical phonon scattering: becomes significant at higher temperatures
    # 3. Impurity scattering: relatively temperature independent
    # 4. Remote interfacial scattering: temperature dependent
    
    # At room temperature and high mobility substrates, 
    # acoustic phonon scattering often dominates
    mu_acoustic = 10000 * (300/T)**1.5  # in cm²/Vs, normalized to 10000 at 300K
    
    # At low temperatures, impurity scattering dominates
    mu_impurity = 15000 * np.ones_like(T)  # Temperature independent
    
    # Combined mobility using Matthiessen's rule: 1/μ_total = 1/μ_acoustic + 1/μ_impurity
    mu_total = 1.0 / (1.0/mu_acoustic + 1.0/mu_impurity)
    
    return T, mu_total, mu_acoustic, mu_impurity

def plot_mobility_vs_temperature():
    """
    Plot carrier mobility vs temperature for graphene
    """
    T_range = np.linspace(10, 500, 100)  # Temperature in Kelvin
    
    T, mu_total, mu_acoustic, mu_impurity = calculate_mobility_vs_temperature(T_range)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(T, mu_total, 'b-', linewidth=2, label='Total Mobility')
    ax.semilogy(T, mu_acoustic, 'r--', linewidth=2, label='Acoustic Phonon Scattering')
    ax.semilogy(T, mu_impurity, 'g:', linewidth=2, label='Impurity Scattering')
    
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Carrier Mobility (cm²/Vs)')
    ax.set_title('Carrier Mobility vs Temperature in Single-Layer Graphene')
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.savefig('mobility_vs_temperature.png', dpi=300, bbox_inches='tight')
    plt.show()

def calculate_quantum_conductance():
    """
    Calculate the quantum of conductance in graphene
    Due to Klein tunneling and unique band structure, graphene shows unique transport properties
    """
    # Conductance quantum
    G0 = 2 * e**2 / hbar  # Factor of 2 for spin degeneracy in graphene
    
    # At the charge neutrality point (Dirac point), minimum conductivity is 4e²/h
    G_min = 4 * e**2 / hbar  # This is a hallmark of graphene
    
    return G0, G_min

def plot_conductance_map(Vg_range, Vd_range):
    """
    Plot a conductance map showing the characteristic "pinch-off" behavior of graphene FETs
    """
    Vg, Vd = np.meshgrid(Vg_range, Vd_range)
    
    # Simulate graphene FET characteristics
    # Charge neutrality point occurs when gate voltage balances doping
    V_neutral = 0  # Assume undoped graphene for simplicity
    
    # Carrier density as function of gate voltage
    n_carrier = (Vg - V_neutral) * 1e12  # Simplified linear relationship
    
    # Conductance as function of carrier density
    # Minimum conductance at neutrality point
    conductance = 1.0 + np.abs(n_carrier) * 0.1  # Add minimum conductance of 1 in arbitrary units
    
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.contourf(Vg, Vd, conductance, levels=50, cmap='viridis')
    contour_lines = ax.contour(Vg, Vd, conductance, levels=10, colors='white', alpha=0.3, linewidths=0.5)
    ax.clabel(contour_lines, inline=True, fontsize=8)
    
    ax.set_xlabel('Gate Voltage (V)')
    ax.set_ylabel('Drain Voltage (V)')
    ax.set_title('Conductance Map of Graphene Field-Effect Transistor')
    plt.colorbar(im, label='Conductance (a.u.)')
    plt.tight_layout()
    plt.savefig('conductance_map.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_klein_tunneling():
    """
    Illustrate Klein tunneling phenomenon in graphene
    This is a unique quantum effect due to the massless Dirac fermions
    """
    # Transmission probability through a potential barrier
    # For normal incidence: T = 1 (perfect transmission - Klein tunneling)
    # For oblique incidence: T = cos²θ where θ is the incidence angle
    
    angles = np.linspace(0, np.pi/2, 1000)  # Incidence angles from 0 to π/2
    transmission_normal = np.ones_like(angles)  # Perfect transmission for normal incidence
    
    # For massive particles, transmission would be exponentially suppressed
    # But for massless Dirac fermions in graphene, transmission is angle-dependent
    # T = cos²θ for perfectly transparent barrier
    transmission_graphene = np.cos(angles)**2
    
    # For comparison, plot exponential suppression for massive particles
    barrier_height = 0.5  # eV
    barrier_width = 1e-9  # m
    hbar_ev = 6.582119569e-16  # eV*s
    m_eff = 0.01 * 9.10938356e-31  # Effective mass (fraction of electron mass)
    
    # Classical transmission through barrier (for massive particles)
    # This is just an approximation to show the contrast
    transmission_massive = np.exp(-2 * angles * 5)  # Exponential suppression
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(np.degrees(angles), transmission_normal, 'g-', linewidth=2, label='Perfect Transmission (Normal Incidence)')
    ax.plot(np.degrees(angles), transmission_graphene, 'b-', linewidth=2, label='Graphene (Massless Fermions)')
    ax.plot(np.degrees(angles), transmission_massive, 'r--', linewidth=2, label='Massive Particles (Classical)')
    
    ax.set_xlabel('Incidence Angle (degrees)')
    ax.set_ylabel('Transmission Probability')
    ax.set_title('Klein Tunneling in Graphene vs Classical Tunneling')
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.savefig('klein_tunneling.png', dpi=300, bbox_inches='tight')
    plt.show()

def calculate_optical_conductivity(freq_range):
    """
    Calculate optical conductivity of graphene
    Universal value at low frequencies: σ₀ = πe²/(2h) ≈ (4π/α)⁻¹
    where α is the fine structure constant
    """
    # Universal optical conductivity of graphene
    sigma0 = np.pi * e**2 / (2 * hbar)  # ≈ 4π/α⁻¹ where α is fine structure constant
    
    # Frequency-dependent correction (simplified)
    omega = 2 * np.pi * freq_range
    hbar_omega = hbar * omega  # Energy of photons
    
    # At finite chemical potential μ, there's interband contribution
    mu = 0.1  # Chemical potential in eV
    
    # Interband transitions become significant when photon energy > 2μ
    interband_factor = np.where(hbar_omega > 2 * mu, 1.0, 0.5)
    
    sigma_real = sigma0 * interband_factor
    
    return freq_range, sigma_real

def plot_optical_absorption():
    """
    Plot the universal optical absorption of graphene
    Graphene absorbs πα ≈ 2.3% of incident light per layer
    """
    # Wavelength range from UV to IR
    wavelengths = np.logspace(-1, 3, 1000)  # nm
    freqs = 3e8 / (wavelengths * 1e-9)  # Hz
    
    freq_range, sigma_real = calculate_optical_conductivity(freqs)
    
    # Absorption coefficient related to conductivity
    # For normal incidence: A = πα ≈ 2.3% per layer
    absorption_percent = np.full_like(wavelengths, np.pi * 1/137 * 100)  # πα in percent
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    # Plot absorption vs wavelength
    ax1.semilogx(wavelengths, absorption_percent, 'purple', linewidth=2)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Absorption (%)')
    ax1.set_title('Universal Optical Absorption of Single-Layer Graphene (~2.3%)')
    ax1.grid(True, alpha=0.3)
    
    # Plot conductivity vs frequency
    ax2.loglog(freqs, sigma_real, 'orange', linewidth=2)
    ax2.set_xlabel('Frequency (Hz)')
    ax2.set_ylabel('Optical Conductivity (S)')
    ax2.set_title('Optical Conductivity of Graphene')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('optical_properties.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    print("Calculating and plotting transport properties of graphene...")
    
    # Calculate Fermi velocity
    v_F = calculate_fermi_velocity()
    print(f"Fermi velocity in graphene: {v_F:.2e} m/s")
    
    # Plot linear dispersion
    print("Plotting linear dispersion relation...")
    plot_linear_dispersion()
    
    # Plot quantum Hall effect
    print("Plotting quantum Hall effect...")
    plot_quantum_hall_effect()
    
    # Plot mobility vs temperature
    print("Plotting mobility vs temperature...")
    plot_mobility_vs_temperature()
    
    # Plot conductance map
    print("Plotting conductance map...")
    Vg_range = np.linspace(-10, 10, 100)
    Vd_range = np.linspace(-1, 1, 50)
    plot_conductance_map(Vg_range, Vd_range)
    
    # Plot Klein tunneling
    print("Plotting Klein tunneling effect...")
    plot_klein_tunneling()
    
    # Plot optical properties
    print("Plotting optical properties...")
    plot_optical_absorption()
    
    print("All transport property plots generated successfully!")