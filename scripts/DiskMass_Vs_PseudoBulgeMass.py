# Load data arrays #
PBMass = np.load(SavePath + 'PBMass.npy')
DiskMass = np.load(SavePath + 'DiskMass.npy')

# Trim the data #
index = np.where(PBMass > 0.0)

PBMass = PBMass[index] * MassUnits
DiskMass = DiskMass[index] * MassUnits

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure = plt.figure(0, figsize=(10, 7.5))

# Figure parameters #
plt.xscale('log')
plt.yscale('log')

plt.xlim(1e7, 1e12)
plt.ylim(1e8, 1e12)

plt.ylabel(r'$\mathrm{M_{pb} / M_{\odot}}$')
plt.xlabel(r'$\mathrm{M_{d,\bigstar} / M_{\odot}}$')

plt.tick_params(direction='in', which='both', top='on', right='on')

######################################################################################################################################################

# Read observational data from G09 #
G09 = np.genfromtxt('./Obs_Data/G09.txt', delimiter=';', skip_header=1,
                    dtype=[('Type', 'U15'), ('Type2', 'U15'), ('h', 'f8'), ('re', 'f8'), ('BtoT', 'f8'), ('MtoLb', 'f8'), ('MtoLd', 'f8'),
                           ('Mb', 'f8'), ('Md', 'f8')])

# Trim observational data from G09 #
index = np.where((G09['Type'] == 'pseudo-bulge'))

# Plot observational data from G09 #
scatter = plt.scatter(G09['Md'][index], G09['Mb'][index], c='lime', edgecolor='black', s=50, zorder=2, label=r'$\mathrm{Gadotti\, 09: PBs}$')

# Create the legends #
colors = ['black', 'grey', 'lightgrey']
squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
legend = plt.legend([squares], [r'$\mathrm{This\;work}$'], scatterpoints=len(colors), scatteryoffsets=[.5], handlelength=len(colors), markerscale=2,
                    frameon=False, loc=1)

plt.gca().add_artist(legend)
plt.legend(frameon=False, scatterpoints=3, loc=2)

# Plot L-Galaxies data - 2D histogram #
hexbin = plt.hexbin(DiskMass, PBMass, xscale='log', yscale='log', bins='log', cmap='Greys', mincnt=2, gridsize=150)

# Adjust the color bar #
cbaxes = figure.add_axes([0.9, 0.11, 0.02, 0.77])
cb = plt.colorbar(hexbin, cax=cbaxes)
cb.set_label('$\mathrm{Counts\; per\; hexbin}$')

######################################################################################################################################################

# Save the figure #
plt.savefig('DM_Vs_PBM_56-' + date + '.pdf', bbox_inches='tight')