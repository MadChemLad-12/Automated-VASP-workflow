import os
import subprocess
from ase import Atoms
from ase.io import write
from ase.build import surface, bulk

# Define the working directory for the VASP calculations
working_directory = "vasp_surface_calc"
if not os.path.exists(working_directory):
    os.makedirs(working_directory)

def generate_surface(structure, miller_indices, layers=3):
    """
    Generate a surface by cleaving the bulk structure along the Miller indices.
    """
    slab = surface(structure, miller_indices, layers)
    slab.center(vacuum=10.0)
    return slab

def generate_input_files(structure, incar_params, kpoints_params, potcar_path):
    """
    Generate VASP input files (POSCAR, INCAR, KPOINTS, POTCAR).
    """
    write("POSCAR", structure)
    
    incar_text = """
    SYSTEM = Surface Calculation
    ENCUT = 520
    ISMEAR = 1
    SIGMA = 0.1
    EDIFF = 1E-5
    IBRION = 2
    NEDOS = 1000
    ALGO = Fast
    """
    with open("INCAR", "w") as f:
        f.write(incar_text)
    
    kpoints_text = f"""
    KPOINTS file
    0
    Monkhorst-Pack
    {kpoints_params[0]} {kpoints_params[0]} {kpoints_params[0]}
    0 0 0
    """
    with open("KPOINTS", "w") as f:
        f.write(kpoints_text)
    
    for element in structure.get_chemical_symbols():
        potcar_file = os.path.join(potcar_path, f"POTCAR_{element}")
        if os.path.exists(potcar_file):
            with open("POTCAR", "ab") as potcar:
                with open(potcar_file, "rb") as element_potcar:
                    potcar.write(element_potcar.read())

def create_slurm_script():
    """
    Create a SLURM job submission script dynamically.
    """
    slurm_script = f"""
#!/bin/bash
#SBATCH --job-name=vasp_surface_calculation
#SBATCH --output=vasp_output_%j.log
#SBATCH --ntasks=32
#SBATCH --nodes=2
#SBATCH --time=48:00:00
#SBATCH --partition=regular
#SBATCH --exclusive
#SBATCH --mem=64GB

# Load the VASP module
module load vasp/5.4.4

# Run VASP using MPI
mpirun -np 32 vasp_std
"""
    with open("vasp_slurm_script.sh", "w") as f:
        f.write(slurm_script)

def run_vasp():
    """
    Submit the SLURM job to the scheduler.
    """
    # Create the SLURM script
    create_slurm_script()

    # Submit the job to SLURM
    try:
        subprocess.run(["sbatch", "vasp_slurm_script.sh"], check=True)
        print("VASP job submitted successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error submitting the SLURM job: {e}")

def main():
    bulk_structure = bulk("Al", "fcc", a=4.05)  # Aluminum, fcc structure
    miller_indices = (1, 1, 1)
    surface_structure = generate_surface(bulk_structure, miller_indices)

    incar_params = {}
    kpoints_params = (4, 4, 4)  # Example k-point mesh
    potcar_path = "/path/to/potcar_files"
    generate_input_files(surface_structure, incar_params, kpoints_params, potcar_path)

    # Submit the SLURM job
    run_vasp()

if __name__ == "__main__":
    main()
