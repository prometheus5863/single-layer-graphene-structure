"""
Main script to run all graphene calculations and generate plots
"""

import subprocess
import sys
import os

def install_requirements():
    """Install required packages"""
    print("Installing required packages...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"])

def run_band_structure():
    """Run the band structure calculations"""
    print("Running band structure calculations...")
    subprocess.check_call([sys.executable, "graphene_band_structure.py"])

def run_transport_properties():
    """Run the transport properties calculations"""
    print("Running transport properties calculations...")
    subprocess.check_call([sys.executable, "graphene_transport_properties.py"])

def main():
    """Main execution function"""
    print("Starting single-layer graphene analysis...")
    print("This will calculate and plot various fundamental properties of graphene.")
    
    try:
        install_requirements()
        run_band_structure()
        run_transport_properties()
        
        print("\nAll calculations completed successfully!")
        print("Generated plots:")
        print("- band_structure.png: Electronic band structure along high symmetry path")
        print("- fermi_surface.png: Fermi surface around the K point") 
        print("- density_of_states.png: Density of states with linear behavior at Dirac point")
        print("- honeycomb_lattice.png: Visualization of the honeycomb lattice structure")
        print("- linear_dispersion.png: Linear dispersion relation near Dirac point")
        print("- quantum_hall_effect.png: Half-integer quantum Hall effect")
        print("- mobility_vs_temperature.png: Carrier mobility vs temperature")
        print("- conductance_map.png: Conductance map of graphene FET")
        print("- klein_tunneling.png: Klein tunneling vs classical tunneling")
        print("- optical_properties.png: Optical absorption and conductivity")
        
    except Exception as e:
        print(f"An error occurred: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())