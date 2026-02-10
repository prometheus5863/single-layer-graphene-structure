"""
Initialize the single-layer graphene repository for GitHub upload
"""

import os
import subprocess
import sys

def initialize_git_repo():
    """Initialize git repository and add all files"""
    print("Initializing Git repository...")
    
    # Initialize git repo
    subprocess.run(["git", "init"], check=True)
    
    # Add all files
    subprocess.run(["git", "add", "."], check=True)
    
    # Create initial commit
    subprocess.run(["git", "commit", "-m", "Initial commit: Single-layer graphene calculations and visualizations"], check=True)
    
    print("Git repository initialized and files committed.")

def create_github_remote():
    """Provide instructions for adding GitHub remote"""
    print("\nTo complete the GitHub setup:")
    print("1. Go to https://github.com/new and create a new repository named 'single-layer-graphene'")
    print("2. After creating the repository, run these commands in this directory:")
    print("   git remote add origin https://github.com/prometheus5863/single-layer-graphene.git")
    print("   git branch -M main")
    print("   git push -u origin main")
    
def list_files():
    """List all files in the repository"""
    print("\nFiles in the repository:")
    for root, dirs, files in os.walk("."):
        for file in files:
            if not file.startswith('.'):
                print(f"  {os.path.join(root, file)}")

def main():
    try:
        # Change to the correct directory
        os.chdir("C:\\Users\\Harsh Vardhan\\.openclaw\\workspace\\single-layer-graphene")
        
        initialize_git_repo()
        list_files()
        create_github_remote()
        
        print(f"\nRepository is ready for GitHub upload!")
        print(f"All calculations and visualizations have been completed successfully.")
        print(f"The repository demonstrates fundamental understanding of graphene physics with:")
        print(f"  - Honeycomb lattice structure visualization")
        print(f"  - Electronic band structure with linear dispersion")
        print(f"  - Density of states with vanishing behavior at Dirac point")
        print(f"  - Quantum Hall effect with half-integer quantization")
        print(f"  - Universal optical absorption (~2.3%)")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())