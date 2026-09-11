# OpenAWSEM constant pH

## AWSEM coarse-grained protein simulations in OpenMM with constant-pH sampling

This repository extends [OpenAWSEM](https://github.com/npschafer/openawsem), an implementation of the AWSEM (Associative memory, Water-mediated Structure, and Energy Model) coarse-grained protein force field for OpenMM, with constant-pH sampling.

The standard OpenAWSEM molecular-dynamics workflow is preserved. The constant-pH extension periodically interrupts the dynamics and attempts a Monte Carlo change in the protonation state of a titratable residue. When a move is accepted, its charge is updated in the Debye-Huckel force without rebuilding the OpenMM system.

## How this fork differs from OpenAWSEM

The constant-pH implementation adds:

- `mm_run_pH.py`, which alternates blocks of molecular dynamics with protonation-state Monte Carlo attempts;
- `Montecarlo_2.py`, which selects residues and evaluates the pH, electrostatic, and local polar-environment contributions used in the acceptance criterion;
- `pH_debyeHuckelTerms.py`, which defines a Debye-Huckel force whose per-particle charges can be changed during a simulation;
- `forces_setup_pH.py`, a force configuration that enables the constant-pH Debye-Huckel term in force group 30; and
- `charge.txt`, generated from the input sequence, which defines the initial residue charges and the residues eligible for constant-pH moves.

The remaining preparation, force-field, simulation, and analysis machinery comes from OpenAWSEM.

## Installation

### Conda

The Conda package installs standard OpenAWSEM, not the constant-pH extension:

```bash
conda install -c conda-forge openawsem
```

To use constant-pH sampling, install this repository from source as described below.

### Git

This installation mode is recommended for users that want to contribute to the code and Wolynes lab members.

```bash
# Clone the constant-pH repository
git clone https://github.com/JoseCuellarBio/openawsem_constant_pH.git
cd openawsem_constant_pH

# Create a new conda environment
conda create -n openawsem -c conda-forge --file requirements.txt
conda activate openawsem

# Install the package in editable mode
pip install -e .
```

## Requirements

### STRIDE
STRIDE is used for secondary structure prediction. 
Download and install STRIDE and add it to your PATH:
https://webclu.bio.wzw.tum.de/stride/
```bash
mkdir stride
cd stride
wget https://webclu.bio.wzw.tum.de/stride/stride.tar.gz
tar -xvzf stride.tar.gz
make
echo 'export PATH=$PATH:'`pwd` >> ~/.bashrc
```

Note: If the webpage above becomes unavailable, please use an alternative repository like https://github.com/MDAnalysis/stride/tree/master/src .

### PSIBLAST    
Install psiblast using the distribution from bioconda:

```bash
conda install -c conda-forge -c bioconda blast
```

Alternatively Download and install psiblast and add it to your PATH: 
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/$(curl -s "https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/" | grep -o 'ncbi-blast-[0-9.]*+-x64-linux.tar.gz'| head -n 1)
tar -xvzf ncbi-*.tar.gz
cd ncbi*/bin
echo 'export PATH=$PATH:'`pwd` >> ~/.bashrc
```

### PDB_SEQRES
* Download pdb_seqres.txt and put it in the cloned openawsem repository location

```bash
wget https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt
OPENAWSEM_LOCATION=$(python -c "import openawsem; print(openawsem.__location__)")
cp pdb_seqres.txt $OPENAWSEM_LOCATION/data
```

### Predict_Property

For secondary structure prediction from the fasta file OpenAWSEM can use "Predict_Property.sh -i {name}.fasta".
Install it from https://github.com/realbigws/Predict_Property.
After installation, add Predict_property.sh to $PATH so it can be executed
For example add 'export PATH = $PATH:/Users/weilu/Research/Build/Predict_Property/' inside the ~/.bash_profile file.

## Configuration
OpenAWSEM allows users to configure data storage paths. To do this:

Create a .awsem directory in your home folder.
Inside .awsem, create a configuration file named config.ini to specify data paths. 
The default paths point to the local data directory inside the OpenAWSEM module.
Example config.ini:

```ini
[Data Paths]
blast = /home/USER/data/database/cullpdb_pc80_res3.0_R1.0_d160504_chains29712
gro = /home/USER/data/Gros
pdb = /home/USER/data/PDBs
index = /home/USER/data/Indices
pdbfail = /home/USER/data/notExistPDBsList
pdbseqres = /home/USER/data/pdb_seqres.txt
topology = /home/USER/topology
```

## Example
Simulation of the amino terminal domain of Phage 434 repressor (1r69)

1. **Activate the OpenMM Environment:**
   Activate the required environment for running simulations.
   ```bash
   source activate openmm
   ```

2. **Set Up the Simulation Folder:**
   Create a simulation folder using the `awsem_create` command. The awsem_create command will automatically download the corresponding pdb.
   ```bash
   awsem_create 1r69 --frag
   ```
   Alternatively, if you have the `1r69.pdb` file:
   ```bash
   awsem_create 1r69.pdb --frag
   ```

3. **Modify the forces_setup.py**

   The `forces_setup.py` script determines which force (energy) terms are included in the simulation. 
   To activate the fragment memory term uncomment the fragment memory term and comment the single memory term.
   ```python
      # templateTerms.fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),
        templateTerms.fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=False),
   ```
   It should look like this:
   ```python
        templateTerms.fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=False),
      #  templateTerms.fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=False),
   ```
4. **Run a standard OpenAWSEM simulation:**
   Execute the simulation using the `awsem_run` command, specifying the platform, number of steps, and start and end temperatures for the annealing simulation.
   As an example we are running 1e5 steps, but it is common to run from 5 to 30 million steps in a single run.
   
   ```bash
   awsem_run 1r69 --platform CPU --steps 1e5 --tempStart 800 --tempEnd 200 -f forces_setup.py
   ```
5. **Run a constant-pH simulation:**

   `awsem_create` also copies the constant-pH scripts and generates `charge.txt` in the simulation directory. Run the local constant-pH driver with the constant-pH force setup:

   ```bash
   ./mm_run_pH.py 1r69 \
       --platform CPU \
       --steps 1e5 \
       --tempStart 300 \
       --simulation_mode 0 \
       --pH 7.0 \
       --interruptFrequency 100 \
       -f forces_setup_pH.py
   ```

   The relevant options are:

   - `--simulation_mode 0`: runs constant-temperature dynamics and activates protonation-state sampling. This option is required; mode `1` performs temperature annealing without constant-pH moves.
   - `--pH`: sets the fixed solution pH used by the Monte Carlo acceptance criterion (default: `7.0`).
   - `--interruptFrequency`: sets the number of MD steps between protonation-state attempts (default: `1000`). Smaller values attempt moves more often and add more overhead.
   - `--tempStart`: sets the MD temperature in constant-temperature mode.
   - `-f forces_setup_pH.py`: enables the updateable Debye-Huckel force required by the constant-pH driver.

   At each interruption, the driver chooses one eligible residue, proposes a change between its charged and neutral states, and accepts or rejects it using a Metropolis criterion. The energy change combines the imposed-pH term, screened electrostatics, and a local polar/nonpolar environment term. Accepted charges are immediately propagated to the OpenMM context.

   Two additional files are written:

   - `protonation_attempts.dat`: one row per attempted move, including the residue index, acceptance flag, pH, and Debye-Huckel energy;
   - `Hawsem.state`: charge-state snapshots recorded when an interruption coincides with the trajectory reporting interval.

### Configuring titratable residues

`charge.txt` contains zero-based residue indices and initial charges:

```text
0 0.0
1 1.0
2 -1.0
```

The file is generated automatically from the FASTA sequence. In the current workflow, only entries with a nonzero initial charge are placed in the Monte Carlo candidate list. The generated defaults therefore sample Arg and Lys from `+1` to `0`, and Asp and Glu from `-1` to `0`. Although the Monte Carlo module contains parameters for additional residue types, neutral entries are not selected by the current driver. Review `charge.txt` before starting a production run; its residue numbering must match the OpenAWSEM system.

The constant-pH force must remain in force group 30 because `mm_run_pH.py` uses that group to locate the force, update its particle charges, and report its energy.

> **Current model assumptions:** protonation moves use residue-specific intrinsic pKa values and an internal Monte Carlo temperature of 300 K. That acceptance temperature is currently independent of `--tempStart`. Protonation is represented by changing coarse-grained charges; explicit protons are not added to the structure.

6. **Compute Energy and Q:**
   Analyze the simulation results and redirect the output to `info.dat`.
   ```bash
   awsem_analyze 1r69 > info.dat
   ```

7. **Run Local Scripts (Optional):**
   The scripts are copied to the project folder and can be modified as needed. To run the local scripts, use the following commands:
   ```bash
   ./mm_run.py 1r69 --platform CPU --steps 1e5 --tempStart 800 --tempEnd 200 -f forces_setup.py
   ./mm_analyze.py 1r69 > energy.dat
   ```

## CLI Tools

### Fix Amino Acid Names (`awsem fix_aminoacids`)

Converts AWSEM placeholder residue names (NGP, IGL, IPR) to standard amino acids using a sequence file.

```bash
awsem fix_aminoacids movie.pdb -f crystal_structure.fasta -o movie_standard.pdb
```

### All-Atom Reconstruction (`awsem reconstruct`)

Reconstructs all-atom structures from AWSEM+3SPN2 coarse-grained models using SCWRL4 (protein) and DNAbackmap (DNA).

```bash
# Protein + DNA
awsem reconstruct model.pdb -f protein.seq --scwrl /path/to/Scwrl4 --dnabackmap /path/to/DNAbackmap

# Protein only (no DNAbackmap needed)
awsem reconstruct protein.pdb -f protein.seq --scwrl /path/to/Scwrl4
```

**Requirements** [SCWRL4](https://dunbrack.fccc.edu/lab/scwrl), [DNAbackmap](https://www.cafemol.org/download/) (DNA only)

## Notes
AWSEM is capable of modeling protein-DNA interactions when used together with open3SPN2, which can be found in a separate package at https://github.com/cabb99/open3spn2.

For small proteins, the LAMMPS version may be faster than OpenAWSEM, especially if a GPU is unavailable. Consider using http://awsem-md.org for such cases.

A quick check of the stability of a protein in AWSEM can be done using the frustratometer server http://frustratometer.qb.fcen.uba.ar/

## Acknowledgements
This project is also supported by the Center for Theoretical Biological Physics (NSF Grants PHY-2019745 and PHY-1522550), with additional support from the D.R. Bullard Welch Chair at Rice University (Grant No. C-0016 to PGW).  We thank AMD (Advanced Micro Devices, Inc.) for the donation of high-performance computing hardware and HPC resources.  Carlos Bueno was supported by the MolSSI Software Fellowship. The skeleton of this project is based on the [Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.11.


## Data availability
Data related to the paper "OpenAWSEM with Open3SPN2: A fast, flexible, and accessible framework for large-scale coarse-grained biomolecular simulations" is available at https://app.globus.org/file-manager?origin_id=b4cef8ce-7773-4016-8513-829f388f7986&origin_path=%2FopenAWSEM_data%2F

## Citation
Please cite the following paper when using OpenAWSEM:
Lu, W., Bueno, C., Schafer, N. P., Moller, J., Jin, S., Chen, X., ... & Wolynes, P. G. (2021). OpenAWSEM with Open3SPN2: A fast, flexible, and accessible framework for large-scale coarse-grained biomolecular simulations. PLoS computational biology, 17(2), e1008308. https://doi.org/10.1371/journal.pcbi.1008308
