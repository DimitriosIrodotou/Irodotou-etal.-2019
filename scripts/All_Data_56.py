# Path to the output folder #
OutputPath = '/Volumes/BAM-BLACK/output/'

# Define the number of tree files, their name and snap number #
firstFile = 0
lastFile = 511
maxFile = 512
filePostfix = '.h5'
filePrefix = 'SA_output_'

# Determine the size and declare arrays to hold the data #
nGal = 0
for iFile in range(firstFile, lastFile + 1):
    # The following line closes the file at the end of the loop
    with h5py.File(OutputPath + filePrefix + '%i' % iFile + filePostfix, 'r') as f:
        nGal += len(f[snap])
    print('nGal =', nGal)

# Declare numpy arrays to hold the data #
PBMass = np.empty(nGal)
CBMass = np.empty(nGal)
DiskMass = np.empty(nGal)
BulgeMass = np.empty(nGal)
DiskRadius = np.empty(nGal)
StellarMass = np.empty(nGal)
CBMassMinor = np.empty(nGal)

# Read in the data #
iGal = 0
for iFile in range(firstFile, lastFile + 1):
    # The following line closes the file at the end of the loop
    with h5py.File(OutputPath + filePrefix + '%i' % iFile + filePostfix, 'r') as f:
        nGalFile = len(f[snap])
        # Masses #
        PBMass[iGal:iGal + nGalFile] = f[snap]['PBMass']
        CBMass[iGal:iGal + nGalFile] = f[snap]['CBMass']
        DiskMass[iGal:iGal + nGalFile] = f[snap]['DiskMass']
        BulgeMass[iGal:iGal + nGalFile] = f[snap]['BulgeMass']
        StellarMass[iGal:iGal + nGalFile] = f[snap]['StellarMass']
        CBMassMinor[iGal:iGal + nGalFile] = f[snap]['CBMassMinor']

        # Sizes #
        DiskRadius[iGal:iGal + nGalFile] = f[snap]['DiskRadius']

        iGal += nGalFile

# Save the arrays so you can load them multiple times in different scripts #
np.save(SavePath + 'PBMass', PBMass)
np.save(SavePath + 'CBMass', CBMass)
np.save(SavePath + 'DiskMass', DiskMass)
np.save(SavePath + 'BulgeMass', BulgeMass)
np.save(SavePath + 'StellarMass', StellarMass)
np.save(SavePath + 'CBMassMinor', CBMassMinor)

np.save(SavePath + 'DiskRadius', DiskRadius)