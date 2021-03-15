# Path to the output folder #
OutputPath = '/Volumes/BAM-BLACK/output/output_ITH_Off/'
# OutputPath = '/Users/Bam/output/'

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
SHMR = np.empty(nGal)
Mvir = np.empty(nGal)
Sfr = np.empty(nGal)
Type = np.empty(nGal)
Vmax = np.empty(nGal)
CBMass = np.empty(nGal)
PBMass = np.empty(nGal)
ColdGas = np.empty(nGal)
DiskMass = np.empty(nGal)
BulgeMass = np.empty(nGal)
DiskRadius = np.empty(nGal)
InfallVmax = np.empty(nGal)
CBMassMajor = np.empty(nGal)
CBMassMinor = np.empty(nGal)
StellarMass = np.empty(nGal)
BlackHoleMass = np.empty(nGal)
HaloSpin = np.empty([nGal, 3])
DiskSpin = np.empty([nGal, 3])
BulgeSpin = np.empty([nGal, 3])

# Read in the data #
iGal = 0
for iFile in range(firstFile, lastFile + 1):
    # The following line closes the file at the end of the loop
    with h5py.File(OutputPath + filePrefix + '%i' % iFile + filePostfix, 'r') as f:
        nGalFile = len(f[snap])
        # Masses #
        Mvir[iGal:iGal + nGalFile] = f[snap]['Mvir']
        PBMass[iGal:iGal + nGalFile] = f[snap]['PBMass']
        CBMass[iGal:iGal + nGalFile] = f[snap]['CBMass']
        ColdGas[iGal:iGal + nGalFile] = f[snap]['ColdGas']
        DiskMass[iGal:iGal + nGalFile] = f[snap]['DiskMass']
        BulgeMass[iGal:iGal + nGalFile] = f[snap]['BulgeMass']
        StellarMass[iGal:iGal + nGalFile] = f[snap]['StellarMass']
        CBMassMajor[iGal:iGal + nGalFile] = f[snap]['CBMassMajor']
        CBMassMinor[iGal:iGal + nGalFile] = f[snap]['CBMassMinor']
        BlackHoleMass[iGal:iGal + nGalFile] = f[snap]['BlackHoleMass']

        # Sizes #
        DiskRadius[iGal:iGal + nGalFile] = f[snap]['DiskRadius']
        SHMR[iGal:iGal + nGalFile] = f[snap]['StellarHalfMassRadius']

        # Spins #
        HaloSpin[iGal:iGal + nGalFile] = f[snap]['HaloSpin']
        DiskSpin[iGal:iGal + nGalFile] = f[snap]['DiskSpin']
        BulgeSpin[iGal:iGal + nGalFile] = f[snap]['BulgeSpin']

        # Other properties #
        Sfr[iGal:iGal + nGalFile] = f[snap]['Sfr']
        Vmax[iGal:iGal + nGalFile] = f[snap]['Vmax']
        Type[iGal:iGal + nGalFile] = f[snap]['Type']
        InfallVmax[iGal:iGal + nGalFile] = f[snap]['InfallVmax']

        iGal += nGalFile

# Save the arrays so you can load them multiple times in different scripts #
np.save(SavePath + 'Mvir', Mvir)
np.save(SavePath + 'PBMass', PBMass)
np.save(SavePath + 'CBMass', CBMass)
np.save(SavePath + 'ColdGas', ColdGas)
np.save(SavePath + 'DiskMass', DiskMass)
np.save(SavePath + 'BulgeMass', BulgeMass)
np.save(SavePath + 'StellarMass', StellarMass)
np.save(SavePath + 'CBMassMajor', CBMassMajor)
np.save(SavePath + 'CBMassMinor', CBMassMinor)
np.save(SavePath + 'BlackHoleMass', BlackHoleMass)

np.save(SavePath + 'DiskRadius', DiskRadius)
np.save(SavePath + 'StellarHalfMassRadius', SHMR)

np.save(SavePath + 'HaloSpin', HaloSpin)
np.save(SavePath + 'DiskSpin', DiskSpin)
np.save(SavePath + 'BulgeSpin', BulgeSpin)

np.save(SavePath + 'Sfr', Sfr)
np.save(SavePath + 'Vmax', Vmax)
np.save(SavePath + 'Type', Type)
np.save(SavePath + 'InfallVmax', InfallVmax)