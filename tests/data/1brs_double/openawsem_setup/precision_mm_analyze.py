#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import fileinput
import platform
import importlib.util
import MDAnalysis as mda

# __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
# __author__ = 'Wei Lu'

from openawsem import *
from openawsem.helperFunctions.myFunctions import *

# load sequence and get atom indices in lammps topologies corresponding 
# to starts and ends of different chains
seq = []
with open('../1brs.seq','r') as f:
    for line in f:
        seq.append(line.strip())
concat_seq = ''.join(seq)
n_atoms_per_chain = np.array([3*len(chain_seq) for chain_seq in seq]) # all lammps awsemmd residues have exactly 3 particles
n_atoms = np.sum(n_atoms_per_chain)
chain_starts = [1]
for increment in n_atoms_per_chain[:-1]:
    chain_starts.append(chain_starts[-1]+increment)
chain_starts = np.array(chain_starts)
chain_ends = [n_atoms_per_chain[0]]
for increment in n_atoms_per_chain[1:]:
    chain_ends.append(chain_ends[-1]+increment)
chain_ends = np.array(chain_ends)


# load high-precision lammps coordinates into a correctly shaped array
#     lammps only logs CA, CB, and O, but we have to add positions for the virtual sites
coordses = []
coords = None
with open('../dump.lammpstrj','r') as f:
    for counter, line in enumerate(f):
        if counter % (n_atoms+9) < 9: # header lines
            if coords:
                coordses.append(coords)
            coords = []
        else:
            info = line.strip().split(" ")
            try:
                atom_id = int(info[0])
            except:
                breakpoint()
            #try:
            #    assert atom_id == ((counter-8) % (n_atoms)) + 0
            #except:
            #    breakpoint()
            atom_type = int(info[1])
            xyz = [float(number) for number in info[2:]]
            if atom_type == 5: # the "HB" atom type from awsemmd is never included in openawsem topologies
                continue
            elif np.any((0<=atom_id-chain_starts) & (atom_id-chain_starts<=2)): # first residue, which has a different virtual site arrangement than interior ones
                if atom_type == 1: 
                    coords.append(xyz) # these are our CA coords
                elif atom_type == 3:
                    coords.append([0,0,0]) # these are our C' virtual site coords, which will be fixed later
                    coords.append(xyz) # these are our O coords
                elif atom_type == 4:
                    coords.append(xyz) # these are our CB coords
                else:
                    raise AssertionError()
            elif np.any((chain_ends-2<=atom_id) & (atom_id<=chain_ends)): # last residue, which has a different virtual site arrangement than interior ones
                if atom_type == 1:
                    coords.append([0,0,0]) # these are our N virtual site coords, which will be fixed later
                    if concat_seq[(atom_id-1)//3] != "P": 
                        coords.append([0,0,0]) # H virtual site coords, which only exist if not proline
                    coords.append(xyz) # CA coords of the last residue
                elif atom_type == 3:
                    coords.append(xyz) # O coords of the last residue
                elif atom_type == 4:
                    coords.append(xyz) # CB coords of the last residue
                else:
                    raise AssertionError()
            else: # we on an interior residue
                #breakpoint()
                if atom_type == 1:
                    coords.append([0,0,0]) # these are our N virtual site coords, which will be fixed later
                    if concat_seq[(atom_id-1)//3] != "P":
                        coords.append([0,0,0]) # H virtual site coords, which only exist if not proline
                    coords.append(xyz) # CA coords 
                elif atom_type == 3:
                    coords.append([0,0,0]) # these are our C' virtual site coords, which will be fixed later
                    coords.append(xyz) # O coords
                elif atom_type == 4:
                    coords.append(xyz) # CB coords
                else:
                    raise AssertionError()
coordses.append(coords)
coordses = np.array(coordses)/10 # don't forget the unit conversion!



# normal analysis setup, except we're loading a lammps format trajectory
forceSetupFile = 'forces_setup.py'
platform = Platform.getPlatformByName('Reference')
chain = getAllChains("crystal_structure.pdb")
seq = read_fasta('crystal_structure.fasta')
u = mda.Universe('../dump.lammpstrj', topology_format='LAMMPSDUMP') 
oa = OpenMMAWSEMSystem('1brs-openmmawsem.pdb', chains=chain, k_awsem=1.0, xml_filename=openawsem.xml, seqFromPdb=seq)
print(f"using force setup file from {forceSetupFile}")
spec = importlib.util.spec_from_file_location("forces", forceSetupFile)
# print(spec)
forces = importlib.util.module_from_spec(spec)
spec.loader.exec_module(forces)
forces = forces.set_up_forces(oa, computeQ=False)
oa.addForcesWithDefaultForceGroup(forces)
# print(forces)
# start simulation
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
# apply forces
forceGroupTable = {"Con":0, "Chain":1, "Chi":2, "Excl":3, "Rama":4, "Rama_P":5, "Rama_ssweight":6, "Contact":22, "Fragment":23, 
                   "Debye_huckel":30, "Rama_Total":[4,5,6],
                    "Total":list(range(32))}
print("Please ensure the forceGroupTable in mm_analysis is set up correctly if you are adding new energy terms.")
showValue = []
# term in showEnergy will assume to take on the energy unit of kilojoule_per_mole, it will be shown in unit of kilocalories_per_mole(divided by 4.184) 
# term in showValue will not be converted.
showEnergy = ["Con", "Chain", "Chi", "Excl", "Rama", "Rama_P", "Rama_ssweight", "Contact", "Fragment", "Debye_huckel","Total", "Rama_Total"]
showAll = showValue + showEnergy



# print energies, using positions saved in the coordses array
print("Printing energies")
with open('info2.dat', "w") as out:
    line = " ".join(["{0:<8s}".format(i) for i in ["Steps"] + showAll])
    print(line)
    out.write(line+"\n")
    for frame_index in range(coordses.shape[0]):
        simulation.context.setPositions(coordses[frame_index,:,:])
        simulation.context.computeVirtualSites() 
        e = []
        for term in showAll:
            if type(forceGroupTable[term]) == list:
                g = set(forceGroupTable[term])
            elif forceGroupTable[term] == -1:
                g = -1
            else:
                g = {forceGroupTable[term]}
                #if g == {4}:
                #    print('4 found!')
            state = simulation.context.getState(getEnergy=True, groups=g)
            #if g== {4}:
            #    print(state.getPotentialEnergy())
            # if term == "Q" or term == "Rg" or term == "Qc" or term == "Q_wat" or term == "Q_mem":
            if term in showValue:
                termEnergy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
            else:
                termEnergy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
            e.append(termEnergy)
    #     print(*e)
        line = " ".join([f"{frame_index:<8}"] + ["{0:<8.20f}".format(i) for i in e])
        print(line)
        out.write(line+"\n")
    #         print(forceGroupTable[term], state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
