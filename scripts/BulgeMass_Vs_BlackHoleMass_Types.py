# Load data arrays #
PBMass = np.load(SavePath + 'PBMass.npy')
CBMass = np.load(SavePath + 'CBMass.npy')
BulgeMass = np.load(SavePath + 'BulgeMass.npy')
CBMassMajor = np.load(SavePath + 'CBMassMajor.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')
BlackHoleMass = np.load(SavePath + 'BlackHoleMass.npy')

# Trim the data #
indexID = np.where((PBMass > 0.0) & (BlackHoleMass > 0))
indexMD = np.where((CBMass > 0.0) & (BlackHoleMass > 0))
indexEll = np.where((CBMassMajor > 0.0) & (BlackHoleMass > 0))

Bulge_Mass3 = BulgeMass[indexID] * MassUnits
Bulge_Mass2 = BulgeMass[indexMD] * MassUnits
Bulge_Mass1 = BulgeMass[indexEll] * MassUnits
BH_Mass3 = BlackHoleMass[indexID] * MassUnits
BH_Mass2 = BlackHoleMass[indexMD] * MassUnits
BH_Mass1 = BlackHoleMass[indexEll] * MassUnits

dlog10 = 0.3

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure, ax = plt.subplots()
figure, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True, figsize=(10, 15))
figure.subplots_adjust(hspace=0)

# Figure parameters #
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e9, 1e12)
plt.ylim(1.1e6, 1e11)

ax3.set_xlabel(r'$\mathrm{M_{b} / M_{\odot}}$')
ax2.set_ylabel(r'$\mathrm{M_{\bullet} / M_{\odot}}$')
ax1.tick_params(direction='in', which='both', top='on', right='on')
ax2.tick_params(direction='in', which='both', top='on', right='on')
ax3.tick_params(direction='in', which='both', top='on', right='on')

######################################################################################################################################################

# Read observational data from KH13_EL, KH13_CB and KH13_PB #
KH13_EL = np.genfromtxt('./Obs_Data/KH13_EL.csv', delimiter=',', names=['Mb', 'Mbh'])
KH13_CB = np.genfromtxt('./Obs_Data/KH13_CB.csv', delimiter=',', names=['Mb', 'Mbh'])
KH13_PB = np.genfromtxt('./Obs_Data/KH13_PB.csv', delimiter=',', names=['Mb', 'Mbh'])

# Plot observational data from KH13_EL, KH13_CB and KH13_PB #
ax1.scatter(KH13_EL['Mb'], KH13_EL['Mbh'], color='red', edgecolor='black', s=3 * size, marker='*', label=r'$\mathrm{Kormendy+13: Els}$', zorder=2)

# Plot L-Galaxies data - 2D histogram #
h = ax1.hexbin(Bulge_Mass1, BH_Mass1, xscale='log', yscale='log', bins='log', cmap='Reds', mincnt=1, gridsize=50)

# Adjust the color bar #
cbaxes = figure.add_axes([0.9, 0.623, 0.02, 0.257])
cb = plt.colorbar(h, cax=cbaxes)
cb.set_label('$\mathrm{Counts\; per\; hexbin}$')

# Create the legends #
colors = ['darkred', 'red', 'lightgrey']
squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
legend1 = ax1.legend([squares], [r'$\mathrm{This\;work:M_{cb(ma)} > 0.0}$'], scatterpoints=len(colors), scatteryoffsets=[.5],
                     handlelength=len(colors), markerscale=2, frameon=False, loc=2)

ax1.add_artist(legend1)
ax1.legend(frameon=False, loc=1, scatterpoints=sp)

######################################################################################################################################################

ax2.scatter(KH13_CB['Mb'], KH13_CB['Mbh'], color='orange', edgecolor='black', s=3 * 50, marker='*', label=r'$\mathrm{Kormendy+13: CBs}$', zorder=2)

# Plot L-Galaxies data - 2D histogram #
h = ax2.hexbin(Bulge_Mass2, BH_Mass2, xscale='log', yscale='log', bins='log', cmap='Oranges', mincnt=2, gridsize=100)

# Adjust the color bar #
cbaxes = figure.add_axes([0.9, 0.366, 0.02, 0.257])
cb = plt.colorbar(h, cax=cbaxes)
cb.set_label('$\mathrm{Counts\; per\; hexbin}$')

# Create the legends #
colors = ['darkorange', 'orange', 'lightgrey']
squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
legend2 = ax2.legend([squares], [r'$\mathrm{This\;work:M_{cb} > 0.0}$'], scatterpoints=len(colors), scatteryoffsets=[.5], handlelength=len(colors),
                     markerscale=2, frameon=False, loc=2)

ax2.add_artist(legend2)
ax2.legend(frameon=False, loc=1, scatterpoints=3)

########################################################################################################################

ax3.scatter(KH13_PB['Mb'], KH13_PB['Mbh'], color='g', edgecolor='black', s=3 * 50, marker='*', label=r'$\mathrm{Kormendy+13: PBs}$', zorder=2)

# Plot L-Galaxies data - 2D histogram #
h = ax3.hexbin(Bulge_Mass3, BH_Mass3, xscale='log', yscale='log', bins='log', cmap='Greens', mincnt=2, gridsize=50)

# Adjust the color bar #
cbaxes = figure.add_axes([0.9, 0.11, 0.02, 0.257])
cb = plt.colorbar(h, cax=cbaxes)
cb.set_label('$\mathrm{Counts\; per\; hexbin}$')

# Create the legends #
colors = ['darkgreen', 'green', 'lightgrey']
squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
legend3 = ax3.legend([squares], [r'$\mathrm{This\;work:M_{pb} > 0.0}$'], scatterpoints=len(colors), scatteryoffsets=[.5], handlelength=len(colors),
                     markerscale=2, frameon=False, loc=2)

ax3.add_artist(legend3)
ax3.legend(frameon=False, loc=1, scatterpoints=3)

######################################################################################################################################################

# Save the figure #
plt.savefig('BM_Vs_BHM_Types_58-' + date + '.pdf', bbox_inches='tight')