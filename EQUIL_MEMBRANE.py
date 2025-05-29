"""
Equilibrate usage
python EQUIL_MEMBRANE.py input_directory name output_directory
	input_directory: absolute path to directory where input OpenMM files are
	name: identifier for OpenMM.files without the extension 
	output_directory: absolute path to directory where subdirectory with equilibration steps will be saved
"""
from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np
import os, sys
sys.path.append('MotorRow')
from MotorRow import MotorRow
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input_dir', help='Directory with .pdb and .xml files for input')
parser.add_argument('name', help='Name to find .pdb and .xml. Should in format $INPUT_DIR/$NAME.pdb and $INPUT_DIR/$NAME.xml')
parser.add_argument('output_dir', help='Directory to store equilibration files')
parser.add_argument('--lig-resname', type=str, default=None, required=False, help='Resname of the ligand')
parser.add_argument('--lig-chain', type=str, default=None, required=False, help='Chain of the ligand (if peptide)')
parser.add_argument('--step-5-nsteps', type=int, default=1250000, required=False, help='Optional. Specify no. of 2 fs steps to perform in step 5 (NPT) no restraints')
args = parser.parse_args()

# Input files
input_dir = args.input_dir
input_name = args.name
input_xml = os.path.join(input_dir, input_name + '.xml')
input_pdb = os.path.join(input_dir, input_name + '.pdb')
lig_resname = args.lig_resname
lig_chain = args.lig_chain

# Output files
output_dir = args.output_dir
if not os.path.exists(output_dir):
	os.mkdir(output_dir)
output_dir = os.path.join(sys.argv[3], input_name)

# Run equilibration
final_xml, final_pdb = MotorRow(input_pdb, input_xml, output_dir, lig_resname=lig_resname, lig_chain=lig_chain).main(input_pdb, step_5_nsteps=args.step_5_nsteps)

# Print final files
print('Final pdb saved to', final_pdb)
print('Final xml saved to', final_xml)
