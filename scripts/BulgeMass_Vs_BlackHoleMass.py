# Load data arrays #
BulgeMass = np.load(SavePath + 'BulgeMass.npy')
BlackHoleMass = np.load(SavePath + 'BlackHoleMass.npy')

# Trim the data #
index = np.where((BulgeMass > 0.0) & (BlackHoleMass > 0.0))

Bulge_Mass = np.log10(BulgeMass[index] * MassUnits * hubble ** 2)
Black_Hole_Mass = np.log10(BlackHoleMass[index] * MassUnits * hubble)

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure = plt.figure(0, figsize=(10, 7.5))

# Figure parameters #
plt.ylim(5.8, 10.5)
plt.xlim(8.5, 12.5)
plt.xlabel(r'$\mathrm{log_{10}(M_{b} / h^{-2}M_{\odot})}$')
plt.ylabel(r'$\mathrm{log_{10}(M_{\bullet} / h^{-1}M_{\odot})}$')
plt.tick_params(direction='in', which='both', top='on', right='on')

######################################################################################################################################################

# Read observational data from McConnell12 and HR04 #
MM04 = Table.read('./Obs_Data/McConnell12_BHBM.dat', format='ascii', data_start=20)
HR04 = np.genfromtxt('./Obs_Data/HR04.csv', delimiter=',', names=['Mb', 'Mbh', 'yplus', 'yminus'])
BN18 = np.genfromtxt('./Obs_Data/BN18.csv', delimiter=',', names=['Mb', 'xerr', 'Mbh', 'yplus', 'yminus'])

# Plot observational data from McConnell12, HR04 and BM18 #
yerr = [np.log10(HR04['Mbh'] / HR04['yminus']), np.log10(HR04['yplus'] / HR04['Mbh'])]

plt.errorbar(np.log10(HR04['Mb'] * 1e10 * hubble ** 2), np.log10(HR04['Mbh'] * 1e10 * hubble), yerr=yerr, fmt='s', ecolor='blue', color='blue',
             markersize=3, label=r'$\mathrm{H\ddot{a}ring\; &\; Rix\, 04}$', zorder=3)

xerr = np.zeros(len(np.log10(MM04['col14'] * hubble ** 2)), dtype=np.float64) + 0.24
yerr = [np.log10(MM04['col3'] / MM04['col4']), np.log10(MM04['col5'] / MM04['col3'])]

plt.errorbar(np.log10(MM04['col14'] * hubble ** 2), np.log10(MM04['col3'] * hubble), xerr=xerr, yerr=yerr, fmt='o', ecolor='lime', color='lime',
             markersize=3, label=r'$\mathrm{McConnell\; &\; Ma\, 12}$', zorder=3)

xerr = [BN18['xerr'], BN18['xerr']]
yerr = [BN18['yminus'], BN18['yplus']]
plt.errorbar(BN18['Mb'], BN18['Mbh'], xerr=xerr, yerr=yerr, fmt='s', ecolor='red', color='red', markersize=3,
             label=r'$\mathrm{Bentz\; &\; Manne\operatorname{-}Nicholas\, 18}$', zorder=3)

# Create the legends #
colors = ['black', 'grey', 'lightgrey']
squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
legend = plt.legend([squares], [r'$\mathrm{This\;work}$'], scatterpoints=len(colors), scatteryoffsets=[.5], handlelength=len(colors), markerscale=2,
                    frameon=False, loc=4)

plt.gca().add_artist(legend)
plt.legend(frameon=False, loc=2)

# Plot L-Galaxies data - 2D histogram #
hexbin = plt.hexbin(Bulge_Mass, Black_Hole_Mass, bins='log', cmap="Greys", mincnt=2)

# Adjust the color bar #
cbaxes = figure.add_axes([0.9, 0.11, 0.02, 0.77])
cb = plt.colorbar(hexbin, cax=cbaxes)
cb.set_label('$\mathrm{Counts\; per\; hexbin}$')

######################################################################################################################################################

# Save the figure #
plt.savefig('BM_Vs_BHM_58-' + date + '.pdf', bbox_inches='tight')

# Plot L-Galaxies data - contour #
# xlim = [8.5, 12.5]
# ylim = [5.0, 10.5]
# bin = [0.1, 0.05]
# NGals = len(index[0])
# Nbins = [int((xlim[1] - xlim[0]) / bin[0]), int((ylim[1] - ylim[0]) / bin[1])]

# H, xedges, yedges = np.histogram2d(Bulge_Mass, Black_Hole_Mass, bins=Nbins)
# extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
# mylevels = np.linspace(1., Nbins[1], Nbins[1]) * NGals / (Nbins[1] ** 2 / 0.7)
# H = zoom(H, 20)
# plt.contourf(H.transpose()[::], origin='lower', cmap='inferno', levels=mylevels, extent=extent)
# plt.colorbar(format='%d')