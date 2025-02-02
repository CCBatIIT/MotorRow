#!/bin/bash
#!/bin/bash
#
#SBATCH --job-name=four_fluoro_MDMB_BINACA
#SBATCH --output=/home/bxie4/cannabinoid_receptors/CB2/MotorRow/job/MR_four_fluoro_MDMB_BINACA.out
#SBATCH --error=/home/bxie4/cannabinoid_receptors/CB2/MotorRow/job/MR_four_fluoro_MDMB_BINACA.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodelist=compute-cpu-7
#SBATCH --mem=1GB
#SBATCH --time=7-00:00:00
#SBATCH --partition=normal
#SBATCH --requeue
#SBATCH --chdir=/home/bxie4/cannabinoid_receptors/CB2/MotorRow/four_fluoro_MDMB_BINACA/systems


python EncoreMotorRow.py --name four_fluoro_MDMB_BINACA --resume_org_dcd True --work_dir /home/bxie4/cannabinoid_receptors/CB2/MotorRow/four_fluoro_MDMB_BINACA/systems
