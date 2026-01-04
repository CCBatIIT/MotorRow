#Simulate the systems
#USAGE python simulate_from_filepair.py directory
#directory must contain "input_top.pdb" and "input_sys.xml"
#Simulate for 1 nanosecond at a 1 fs timestep recording frames every 1 ps at 300K
import os, sys, glob
import numpy as np
from openmm import *
from openmm.app import *
from openmm.unit import *

def simulate_from_filepair(pdb_top_fn, sys_xml_fn,
                           temp=300, ts=0.001, n_steps=1e6, #1ns
                           dcd_fn='output.dcd', dcd_freq=1000,
                           stdout_fn='output.stdout', stdout_freq=1000,
                           working_dir='./'):
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir, exist_ok=True)
    #Construct
    pdb = PDBFile(pdb_top_fn)
    with open(sys_xml_fn, 'r') as f:
        system = XmlSerializer.deserialize(f.read())
    integrator = LangevinMiddleIntegrator(temp*kelvin, 1/picosecond, ts*picoseconds)
    sim = Simulation(pdb.topology, system, integrator)
    _ = sim.context.setPositions(pdb.positions)
    print('Constructed')
    #Minimize
    _ = sim.minimizeEnergy()
    with open(os.path.join(working_dir,'minimized.pdb'), 'w') as f:
        _ = PDBFile.writeFile(sim.topology, sim.context.getState(getPositions=True).getPositions(), file=f, keepIds=True)
    print('Minimized')
    #Simulate
    _ = sim.reporters.append(DCDReporter(os.path.join(direc, dcd_fn), dcd_freq))
    _ = sim.reporters.append(StateDataReporter(os.path.join(direc, stdout_fn), stdout_freq, step=True, potentialEnergy=True, temperature=True, speed=True))
    _ = sim.step(n_steps)

    #Final
    state = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    contents = XmlSerializer.serialize(state)
    with open(os.path.join(working_dir, 'final_state.xml'), 'w') as f:
        _ = f.write(contents)
    with open(os.path.join(working_dir,'final_top.pdb'), 'w') as f:
        _ = PDBFile.writeFile(sim.topology, sim.context.getState(getPositions=True).getPositions(), file=f, keepIds=True)
    return None



if __name__=='__main__':
    direc = sys.argv[1]
    pdb_file = os.path.join(direc, 'input_top.pdb')
    xml_file = os.path.join(direc, 'input_sys.xml')
    assert os.path.isfile(pdb_file) and os.path.isfile(xml_file)
    print(direc)

    _ = simulate_from_filepair(pdb_file, xml_file,
                               temp=300, ts=0.001, n_steps=1e6, #1ns
                               dcd_fn='output.dcd', dcd_freq=1000,
                               stdout_fn='output.stdout', stdout_freq=1000,
                               working_dir=direc)