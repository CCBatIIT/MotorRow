import os
import sys
sys.path.append('/home/bxie4/cannabinoid_receptors/MotorRow') # Add MotorRow module path
import argparse
from datetime import datetime
import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
import MDAnalysis as mda
from MotorRow_utils import *
import MotorRow as m
#import Bridgeport.Minimizer.Minimizer as m

# Argument parser setup
parser = argparse.ArgumentParser(description="MotorRow simulation script.")
parser.add_argument('--name', required=True, help="Name of the system.")
parser.add_argument('--resume_org_dcd', help="Resume from existing DCD file.")
parser.add_argument('--work_dir', required=True, help="Working directory for simulation files.")
args = parser.parse_args()

name = args.name
RESUME = args.resume_org_dcd
abs_work_dir = args.work_dir

#--------------------------
stepnum = 5     # Stage number
nsteps= 250000  # Total simulation steps
#--------------------------
temp=300.0  # Temperature in Kelvin
press=1.0   # Pressure in bar
dt=2.0  # Time step in femtoseconds
nstdout=1000    # State data reporting interval
ndcd=5000   # DCD reporting interval
ncycles=50  # Number of cycles
steps_per_cycle = nsteps/ncycles
#--------------------------
start = datetime.now()

#abs_work_dir = "/home/bxie4/cannabinoid_receptors/CB2/MotorRow/"+name+"/systems"
input_system = os.path.join(abs_work_dir, name+".xml")
input_state = os.path.join(abs_work_dir, "Step_4.xml")
pdb_file = os.path.join(abs_work_dir,name+".pdb")

resumed_xml_file = os.path.join(abs_work_dir, "Step_5.xml")
resumed_dcd_file = os.path.join(abs_work_dir, "step5.dcd")

m_obj = m.MotorRow(input_system, pdb_file, abs_work_dir)

#--- set up the system -----#
system = XmlSerializer.deserialize(open(input_system, 'r').read())
print(f'Found input_system: {input_system}')
pdb = PDBFile(pdb_file)
topology = pdb.topology

#--- add barostat for stage 5 ---#
if stepnum == 5:
    system.addForce(MonteCarloBarostat(press*bar, temp*kelvin, 100))

#Set up integrator
integrator = LangevinMiddleIntegrator(temp*kelvin, 1/picosecond, dt*femtosecond)
    
# Select platfrom
try:
    platform = Platform.getPlatformByName('OpenCL')
    properties = {'OpenCLPrecision': 'mixed'}
    simulation = Simulation(topology, system, integrator, platform, properties)
except:
    print('Only found CPUs. The job will be terminated.')
    exit(1)
    #simulation = Simulation(topology, system, integrator)
#-----------------------------#
if RESUME:
    latest_state = XmlSerializer.deserialize(open(resumed_xml_file, 'r').read())
    simulation.context.setState(latest_state)
    simulation.context.setPositions(latest_state.getPositions())
    simulation.context.setVelocities(latest_state.getVelocities())
    
    tj = mda.Universe(resumed_dcd_file).trajectory
    frames = len(tj)
    completed_steps = frames * ndcd
    remaining_steps = nsteps - completed_steps
    curr_t = completed_steps*dt*femtosecond/picoseconds
    print(f'{completed_steps} steps have been finished, {remaining_steps} steps are left.')
    out_xml = f'Step_{stepnum}_cont.xml'
    fn_stdout = os.path.join(abs_work_dir, f'step{stepnum}_cont.stdout')
    nsteps = remaining_steps

else:
    simulation.loadState(input_state)
    context = simulation.context
    print(f'Found input_state: {input_state}')
    out_xml = f'Step_{stepnum}.xml'
    fn_stdout = os.path.join(abs_work_dir, f'step{stepnum}.stdout')
    curr_t = 0
fn_dcd = os.path.join(abs_work_dir, f'step{stepnum}.dcd')
SDR = app.StateDataReporter(fn_stdout, nstdout, step=True, time=True,
                            potentialEnergy=True, temperature=True, progress=False,
                            remainingTime=True, speed=True, volume=True,
                            totalSteps=nsteps, separator=' : ')

simulation.reporters.append(SDR)
DCDR = app.DCDReporter(file=fn_dcd, reportInterval=ndcd, append=RESUME)
simulation.reporters.append(DCDR)
print(f'Starting Step {stepnum} with forces {simulation.system.getForces()}')
print(f'Starting Step {stepnum} with box_vectors {simulation.system.getDefaultPeriodicBoxVectors()}')
state_xml_out = os.path.join(abs_work_dir, out_xml)
m_obj._describe_state(simulation, "initial")

# Determine no. of steps per cycle
#steps_per_cycle = int(nsteps / ncycles)
ncycles = int(nsteps/steps_per_cycle)

simulation.context.setTime(curr_t)
for cycle in range(1, ncycles+1):
    print('Cycle', cycle, 'to', ((cycle/ncycles) * ((nsteps * dt) / 1e6)), 'ns')
    simulation.step(steps_per_cycle)
    m_obj._describe_state(simulation, f'Step {stepnum}')
    m_obj._write_state(simulation, state_xml_out)

end = datetime.now() - start
print(f'Step {stepnum} completed after {end}')
print(f'Box Vectors after this step {simulation.system.getDefaultPeriodicBoxVectors()}')

pdb_out = os.path.join(abs_work_dir, f'Step_{stepnum}.pdb')
m_obj._write_structure(simulation, pdb_out)
