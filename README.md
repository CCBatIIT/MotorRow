# MotorRow - Equilibration of prepared membrane protein systems.
## Overview
This module executes a five stage equilibration protocol specific to membrane proteins.  These are:
1. Minimization (Typical OpenMM Minimization)​
2. 250 ps of NVT with restraints​
    1. Restraint 1 – All Protein Heavy Atoms​
    2. Restraint 2 – Z Coordinate of Membrane Heavy Atoms​
3. 250 ps of NPT with Monte Carlo Membrane Barostat and the same restraints as Step 2
4. 250 ps of NVT without restraints​
5. 2500 ps of NPT with Monte Carlo Membrane Barostat​
6. 2500 ps of NPT with Monte Carlo Barostat

Following these steps, the membrane system is ready for production simulation.

## MotorRow Usage 

```python
input_pdb, input_xml = 'built.pdb', 'built.xml' # The outputs of Bridgeport
MR = MotorRow(input_xml, input_pdb, os.path.getcwd())
production_xml, production_pdb = MR.main(input_pdb)
```

### Input files
An XML/PDB Pair which were generated with the main usage case of Bridgeport.

(https://github.com/CCBatIIT/Bridgeport/tree/main)

Note that this XML is an XML serialized OpenMM System, not a force field file

### Output files
An XML/PDB Pair which can now be used in a production simulation.\
Note that this XML is an XML serialized OpenMM System, not a force field file
