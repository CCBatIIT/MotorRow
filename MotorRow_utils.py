import os, shutil
import mdtraj as md
import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *


def get_positions_from_pdb(fname_pdb, lig_resname: str=None, lig_chain: str=None):
    nameMembrane = ['DPP', 'POP']
    with open(fname_pdb, 'r') as f_pdb:
        l_pdb = f_pdb.read().split('\n')

    coords = []
    prt_heavy_atoms = []
    mem_heavy_atoms = []
    lig_heavy_atoms = []
    iatom = 0
    
    for line in l_pdb[:-1]:
        if line[:6] in ['ATOM  ', 'HETATM']:
            resname = line[17:20]
            chain = line[21:23]

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = str(line[76:78].replace('' , ''))

            coords.append(Vec3(x, y, z))

            if line[17:20] in nameMembrane and element != 'H':
                mem_heavy_atoms.append(iatom)
            elif line[:6] in ['ATOM  '] and element != 'H':
                if lig_resname != None and resname == lig_resname and element != 'H':
                    lig_atom_name = line[12:16].strip().strip('x')
                    lig_heavy_atoms.append([iatom, lig_atom_name])
                else:
                    prt_heavy_atoms.append(iatom)
            elif lig_resname != None and resname == lig_resname and element != 'H':
                lig_atom_name = line[12:16].strip().strip('x')
                lig_heavy_atoms.append([iatom, lig_atom_name])

            iatom += 1

    return np.array(coords), prt_heavy_atoms, mem_heavy_atoms, lig_heavy_atoms

def restrain_atoms(system, crds, atom_inds, rst_name: str='fc_pos', rst_strength: float=20.0):

    rest = CustomExternalForce(f'{rst_name}*periodicdistance(x,y,z,x0,y0,z0)^2')
    rest.addGlobalParameter(rst_name, rst_strength)
    rest.addPerParticleParameter('x0')
    rest.addPerParticleParameter('y0')
    rest.addPerParticleParameter('z0')
    for atom_i in atom_inds:
        x, y, z = crds[int(atom_i)] / 10
        rest.addParticle(int(atom_i), [x, y, z])
    system.addForce(rest)

    return system

def unpack_infiles(xml, pdb):
    """
    Parse XML and PDB into Openmm System Topology adn Positions
    """
    print(f'Unpacking {xml}, {pdb}')
    pdb = PDBFile(pdb)
    with open(xml) as f:
        system = XmlSerializer.deserialize(f.read())
    return system, pdb.topology, pdb.positions

def parse_atom_inds(atom_inds, parse_atom_names, find_atom_names):
    
    parsed_atom_inds = np.empty(len(find_atom_names), dtype=int)
    for (atom_i, parse_atom_name) in zip(atom_inds, parse_atom_names):
        if parse_atom_name in find_atom_names:
            find_atom_name_ind = list(find_atom_names).index(parse_atom_name)
            parsed_atom_inds[find_atom_name_ind] = atom_i

    return parsed_atom_inds

def minimize_from_sys(sys, top, pos, temp=300.0, dt=2.0):
    
        integrator = LangevinMiddleIntegrator(temp*kelvin, 1/picosecond, dt*femtosecond)
        simulation = Simulation(top, sys, integrator)
        simulation.context.setPositions(pos)
        _ = _describe_state(simulation, "Original state")
        simulation.minimizeEnergy()
        min_PE = _describe_state(simulation, "Minimized state")

        return simulation, min_PE
