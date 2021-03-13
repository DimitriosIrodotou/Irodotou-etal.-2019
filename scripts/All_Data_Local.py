# Path to the output folder #
OutputPath = '/Users/Bam/Astronomy PhD/L-Galaxies/Development_Branch/output/'

# Define the number of tree files, their name and snap number #
firstFile = 5
lastFile = 5
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
Sfr = np.empty(nGal)
Type = np.empty(nGal)
Vmax = np.empty(nGal)
Vvir = np.empty(nGal)
Rvir = np.empty(nGal)
Mvir = np.empty(nGal)
SHMR = np.empty(nGal)
Fifty = np.empty(nGal)
Ninety = np.empty(nGal)
CBMass = np.empty(nGal)
PBMass = np.empty(nGal)
ColdGas = np.empty(nGal)
DiskMass = np.empty(nGal)
Mag = np.empty([nGal, 5])
BulgeMass = np.empty(nGal)
BulgeSize = np.empty(nGal)
InfallVmax = np.empty(nGal)
DiskRadius = np.empty(nGal)
CBMassMajor = np.empty(nGal)
CBMassMinor = np.empty(nGal)
StellarMass = np.empty(nGal)
ColdGasRadius = np.empty(nGal)
BlackHoleMass = np.empty(nGal)
NMinorMergers = np.empty(nGal)
NMajorMergers = np.empty(nGal)
DiskSpin = np.empty([nGal, 3])
HaloSpin = np.empty([nGal, 3])
BulgeSpin = np.empty([nGal, 3])
ColdGasSpin = np.empty([nGal, 3])

# Read in the data #
iGal = 0
for iFile in range(firstFile, lastFile + 1):
    # The following line closes the file at the end of the loop
    with h5py.File(OutputPath + filePrefix + '%i' % iFile + filePostfix, 'r') as f:
        nGalFile = len(f[snap])
        # Masses #
        PBMass[iGal:iGal + nGalFile] = f[snap]['PBMass']
        CBMass[iGal:iGal + nGalFile] = f[snap]['CBMass']
        ColdGas[iGal:iGal + nGalFile] = f[snap]['ColdGas']
        DiskMass[iGal:iGal + nGalFile] = f[snap]['DiskMass']
        BulgeMass[iGal:iGal + nGalFile] = f[snap]['BulgeMass']
        StellarMass[iGal:iGal + nGalFile] = f[snap]['StellarMass']
        CBMassMajor[iGal:iGal + nGalFile] = f[snap]['CBMassMajor']
        CBMassMinor[iGal:iGal + nGalFile] = f[snap]['CBMassMinor']
        BlackHoleMass[iGal:iGal + nGalFile] = f[snap]['BlackHoleMass']

        # Spins #
        DiskSpin[iGal:iGal + nGalFile] = f[snap]['DiskSpin']
        HaloSpin[iGal:iGal + nGalFile] = f[snap]['HaloSpin']
        BulgeSpin[iGal:iGal + nGalFile] = f[snap]['BulgeSpin']
        ColdGasSpin[iGal:iGal + nGalFile] = f[snap]['ColdGasSpin']

        # Sizes #
        BulgeSize[iGal:iGal + nGalFile] = f[snap]['BulgeSize']
        DiskRadius[iGal:iGal + nGalFile] = f[snap]['DiskRadius']
        SHMR[iGal:iGal + nGalFile] = f[snap]['StellarHalfMassRadius']
        ColdGasRadius[iGal:iGal + nGalFile] = f[snap]['ColdGasRadius']
        # Fifty[iGal:iGal + nGalFile] = f[snap]['StellarHalfLightRadius']
        # Ninety[iGal:iGal + nGalFile] = f[snap]['StellarNinetyLightRadius']

        # Other properties #
        Sfr[iGal:iGal + nGalFile] = f[snap]['Sfr']
        # Mag[iGal:iGal + nGalFile] = f[snap]['Mag']
        Type[iGal:iGal + nGalFile] = f[snap]['Type']
        Vmax[iGal:iGal + nGalFile] = f[snap]['Vmax']
        InfallVmax[iGal:iGal + nGalFile] = f[snap]['InfallVmax']
        # NMajorMergers[iGal:iGal + nGalFile] = f[snap]['NMajorMergers']
        # NMinorMergers[iGal:iGal + nGalFile] = f[snap]['NMinorMergers']

        iGal += nGalFile

# Save the arrays so you can load them multiple times in different scripts #
np.save(SavePath + 'PBMass', PBMass)
np.save(SavePath + 'CBMass', CBMass)
np.save(SavePath + 'ColdGas', ColdGas)
np.save(SavePath + 'DiskMass', DiskMass)
np.save(SavePath + 'BulgeMass', BulgeMass)
np.save(SavePath + 'StellarMass', StellarMass)
np.save(SavePath + 'CBMassMajor', CBMassMajor)
np.save(SavePath + 'CBMassMinor', CBMassMinor)
np.save(SavePath + 'BlackHoleMass', BlackHoleMass)

np.save(SavePath + 'DiskSpin', DiskSpin)
np.save(SavePath + 'HaloSpin', HaloSpin)
np.save(SavePath + 'BulgeSpin', BulgeSpin)
np.save(SavePath + 'ColdGasSpin', ColdGasSpin)

np.save(SavePath + 'BulgeSize', BulgeSize)
np.save(SavePath + 'DiskRadius', DiskRadius)
np.save(SavePath + 'StellarHalfMassRadius', SHMR)
np.save(SavePath + 'ColdGasRadius', ColdGasRadius)
# np.save(SavePath + 'StellarHalfLightRadius', Fifty)
# np.save(SavePath + 'StellarNinetyLightRadius', Ninety)

np.save(SavePath + 'Sfr', Sfr)
# np.save(SavePath + 'Mag', Mag)
np.save(SavePath + 'Type', Type)
np.save(SavePath + 'Vmax', Vmax)
np.save(SavePath + 'InfallVmax', InfallVmax)  # np.save(SavePath + 'NMajorMergers', NMajorMergers)
# np.save(SavePath + 'NMinorMergers', NMinorMergers)