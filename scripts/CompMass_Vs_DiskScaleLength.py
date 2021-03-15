# Load data arrays #
PBMass = np.load(SavePath + 'PBMass.npy')
BulgeMass = np.load(SavePath + 'BulgeMass.npy')
CBMassMinor = np.load(SavePath + 'CBMassMinor.npy')

DiskRadius = np.load(SavePath + 'DiskRadius.npy')

# Trim the data #
indexPB = np.where(PBMass > 0.0)
indexCB = np.where(CBMassMinor > 0.0)

PB_Mass = PBMass[indexPB] * MassUnits
CB_Mass_Minor = CBMassMinor[indexCB] * MassUnits
Disk_Scale_Length_PB = DiskRadius[indexPB] * LengthUnits / 3.
Disk_Scale_Length_CB = DiskRadius[indexCB] * LengthUnits / 3.

dlog10 = 0.2

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure, ax = plt.subplots()
figure, (ax1, ax2) = plt.subplots(nrows=2, sharey=True, figsize=(10, 7.5))
figure.subplots_adjust(hspace=0, wspace=0)

# Figure parameters #
plt.yscale('log')
ax1.set_xscale('log')
ax2.set_xscale('log')

plt.ylim(1e-1, 9e1)
ax1.set_xlim(1e9, 1e12)
ax2.set_xlim(1e9, 1e12)

ax1.set_ylabel(r'$\mathrm{R_{disc} / kpc}$', size=25)
ax2.set_ylabel(r'$\mathrm{R_{disc} / kpc}$', size=25)
ax2.set_xlabel(r'$\mathrm{M_{comp.} / M_{\odot}}$', size=25)

ax1.set_xticklabels([])
ax1.tick_params(direction='in', which='both', top='on', right='on', labelsize=25)
ax2.tick_params(direction='in', which='both', top='on', right='on', labelsize=25)

######################################################################################################################################################

# Read observational data from G09 #
G09 = np.genfromtxt('./Obs_Data/G09.txt', delimiter=';', skip_header=1,
                    dtype=[('Type', 'U15'), ('Type2', 'U15'), ('h', 'f8'), ('re', 'f8'), ('BtoT', 'f8'),
                           ('MtoLb', 'f8'), ('MtoLd', 'f8'), ('Mb', 'f8'), ('Md', 'f8')])

# Trim observational data from G09 #
index = np.where((G09['Type'] == 'classical   ') | (G09['Type2'] == 'classical   '))

# Plot observational data from G09 #
scatter = ax1.scatter(G09['Mb'][index], G09['h'][index], c='magenta', edgecolor='black', s=50, zorder=2,
                      label=r'$\mathrm{Gadotti\, 09: CBs}$')

# Calculate median and 1-sigma #
log10X = np.log10(CB_Mass_Minor)
log10XMax = np.log10(max(CB_Mass_Minor))
log10XMin = np.log10(min(CB_Mass_Minor))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(np.absolute(CB_Mass_Minor)[index])
    if len(index) > 0:
        median[i] = np.nanmedian(Disk_Scale_Length_CB[index])
        slow[i] = np.nanpercentile(Disk_Scale_Length_CB[index], 15.87)
        shigh[i] = np.nanpercentile(Disk_Scale_Length_CB[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = ax1.plot(X, median, color='magenta', lw=lw)
ax1.fill_between(X, shigh, slow, color='magenta', alpha='0.5', zorder=2)
fill, = plt.fill(np.NaN, np.NaN, color='magenta', alpha=0.5)

# Create the legends #
legend1 = ax1.legend([median, fill],
                     [r'$\mathrm{This\; work: Median}$', r'$\mathrm{This\; work:16^{th} - 84^{th}\,\%ile}$'],
                     frameon=False, loc=2)

ax1.add_artist(legend1)
ax1.legend(frameon=False, loc=4, scatterpoints=3)

######################################################################################################################################################

# Trim observational data from G09 #
index = np.where((G09['Type'] == 'pseudo-bulge') | (G09['Type2'] == 'pseudo-bulge'))

# Plot observational data from G09 #
scatter = ax2.scatter(G09['Mb'][index], G09['h'][index], c='green', edgecolor='black', s=50, zorder=2,
                      label=r'$\mathrm{Gadotti\, 09: PBs}$')

# Calculate median and 1-sigma #
log10X = np.log10(PB_Mass)
log10XMax = np.log10(max(PB_Mass))
log10XMin = np.log10(min(PB_Mass))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(np.absolute(PB_Mass)[index])
    if len(index) > 0:
        median[i] = np.nanmedian(Disk_Scale_Length_PB[index])
        slow[i] = np.nanpercentile(Disk_Scale_Length_PB[index], 15.87)
        shigh[i] = np.nanpercentile(Disk_Scale_Length_PB[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = ax2.plot(X, median, color='green', lw=lw)
ax2.fill_between(X, shigh, slow, color='green', alpha='0.5', zorder=2)
fill, = plt.fill(np.NaN, np.NaN, color='green', alpha=0.5)

# Create the legends #
legend3 = ax2.legend([median, fill],
                     [r'$\mathrm{This\; work: Median}$', r'$\mathrm{This\; work:16^{th}-84^{th}\,\%ile}$'],
                     frameon=False, loc=2)

ax2.add_artist(legend3)
ax2.legend(frameon=False, loc=4, scatterpoints=3)

######################################################################################################################################################

# Save the figure #
plt.savefig('CM_Vs_DSL_56-' + date + '.png', bbox_inches='tight')
