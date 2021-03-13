#!/bin/bash -l

# Tell SGE that we are using the bash shell
#$ -S /bin/bash

# Say which queue you want to submit to
#$ -q gpu.q

# Give the job a name
#$ -N Python_DI

source /etc/profile
shopt -s expand_aliases
module load mps/software/
module load hdf5
module load sge
module load gsl/gcc
module load python/3.5.1
source /home/di/di43/Python35/bin/activate

cd /lustre/scratch/astro/di43/Python/scripts
#python All_Data_58.py &
#python All_Data_56.py &
#python All_Data_55.py &
#wait
python SMF_Comp.py
#python SMF_Morph.py &
#python Galaxy_Fraction.py &
#python DiskMass_Vs_DiskSpin.py &
#python StellarMass_Decomposition.py &
#python BulgeMass_Vs_BlackHoleMass.py &
#python StellarMass_Vs_GalacticSpin.py &
#python BulgeMass_Vs_BlackHoleMass_Types.py &
#python Tully_Fisher.py &

#python DiskMass_Vs_PseudoBulgeMass.py &
#python CompMass_Vs_DiskScaleLength.py &
#python StellarMass_Vs_DiskScaleLength.py &

#python StellarMass_Vs_HalfMassRadius_LTG.py &
#python StellarMass_Vs_HalfMassRadius_ETG.py &
#wait
