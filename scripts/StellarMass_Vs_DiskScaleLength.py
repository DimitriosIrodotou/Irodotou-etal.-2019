# Load data arrays #
BulgeMass = np.load(SavePath + 'BulgeMass.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')

DiskRadius = np.load(SavePath + 'DiskRadius.npy')

dlog10 = 0.25

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure = plt.figure(0, figsize=(10, 7.5))

# Figure parameters #
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e9, 1e12)
plt.ylim(6e-1, 1e1)
plt.xlabel(r'$\mathrm{M_{\bigstar} / M_{\odot}}$',size=25)
plt.ylabel(r'$\mathrm{R_{disc} / kpc}$',size=25)
plt.tick_params(direction='in', which='both', top='on', right='on', labelsize=25)

######################################################################################################################################################

# Read observational data from G09 #
G09 = np.genfromtxt('./Obs_Data/G09.txt', delimiter=';', skip_header=1,
                    dtype=[('Type', 'U15'), ('Type2', 'U15'), ('h', 'f8'), ('re', 'f8'), ('BtoT', 'f8'), ('MtoLb', 'f8'), ('MtoLd', 'f8'),
                           ('Mb', 'f8'), ('Md', 'f8')])

# Plot observational data from G09 #
scatter = plt.scatter(G09['Md'] + G09['Mb'], G09['h'], c=G09['BtoT'], edgecolor='black', cmap='RdYlBu_r', s=50, zorder=2)

# Trim the data #
index = np.where((BulgeMass >= 0.0 * StellarMass) & (BulgeMass <= 0.2 * StellarMass))
Stellar_Mass = StellarMass[index] * MassUnits
Disk_Scale_Length = DiskRadius[index] * LengthUnits / 3.0

# Calculate median and 1-sigma #
log10X = np.log10(Stellar_Mass)
log10XMax = 11.7
log10XMin = 9.4
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Stellar_Mass[index])
    if len(index) > 0:
        median[i] = np.median(Disk_Scale_Length[index])
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median1, = plt.plot(X, median, color='blue', lw=lw)

# Trim the data #
index = np.where((BulgeMass > 0.2 * StellarMass) & (BulgeMass <= 0.4 * StellarMass))

Stellar_Mass = StellarMass[index] * MassUnits
Disk_Scale_Length = DiskRadius[index] * LengthUnits / 3.0

# Calculate median and 1-sigma #
log10X = np.log10(Stellar_Mass)
log10XMax = 11.6
log10XMin = 9.5
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Stellar_Mass[index])
    if len(index) > 0:
        median[i] = np.median(Disk_Scale_Length[index])
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median2, = plt.plot(X, median, color='lightblue', lw=lw)

# Trim the data #
index = np.where((BulgeMass > 0.4 * StellarMass) & (BulgeMass <= 0.6 * StellarMass))

Stellar_Mass = StellarMass[index] * MassUnits
Disk_Scale_Length = DiskRadius[index] * LengthUnits / 3.0

# Calculate median and 1-sigma #
log10X = np.log10(Stellar_Mass)
log10XMax = 11.55
log10XMin = 9.6
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Stellar_Mass[index])
    if len(index) > 0:
        median[i] = np.median(Disk_Scale_Length[index])
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median3, = plt.plot(X, median, color='gold', lw=lw)

# Trim the data #
index = np.where((BulgeMass > 0.6 * StellarMass) & (BulgeMass <= 0.8 * StellarMass))

Stellar_Mass = StellarMass[index] * MassUnits
Disk_Scale_Length = DiskRadius[index] * LengthUnits / 3.0

# Calculate median and 1-sigma #
log10X = np.log10(Stellar_Mass)
log10XMax = 11.5
log10XMin = 9.8
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Stellar_Mass[index])
    if len(index) > 0:
        median[i] = np.median(Disk_Scale_Length[index])
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median4, = plt.plot(X, median, color='tomato', lw=lw)

# Create the legends #
legend = plt.legend([median1, median2, median3, median4],
                    [r'$\mathrm{This\; work: 0.0 \leq B/T \leq 0.2}$', r'$\mathrm{This\; work: 0.2 < B/T \leq 0.4}$',
                     r'$\mathrm{This\; work: 0.4 < B/T \leq 0.6}$', r'$\mathrm{This\; work: 0.6 < B/T \leq 0.8}$'], frameon=False, loc=4)

colors = ['darkred', 'gold', 'darkblue']
circles = collections.CircleCollection([10] * 3, facecolor=colors)
legend1 = plt.legend([circles], [r'$\mathrm{Gadotti\, 09}$'], scatterpoints=len(colors), scatteryoffsets=[.5], handlelength=len(colors),
                     markerscale=2, frameon=False, loc=2)

plt.gca().add_artist(legend)
plt.gca().add_artist(legend1)

cbaxes = figure.add_axes([0.9, 0.11, 0.02, 0.77])
cb = plt.colorbar(scatter, cax=cbaxes)
cbaxes.tick_params(direction='out', which='both', right='on', labelsize=25)
cb.set_label(r'$\mathrm{B/T}$')

######################################################################################################################################################

# Save the figure #
plt.savefig('SM_Vs_DSL_56-' + date + '.png', bbox_inches='tight')