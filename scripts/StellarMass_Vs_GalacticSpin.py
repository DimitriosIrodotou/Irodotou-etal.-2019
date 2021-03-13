# Load data arrays #
DiskMass = np.load(SavePath + 'DiskMass.npy')
DiskSpin = np.load(SavePath + 'DiskSpin.npy')
BulgeMass = np.load(SavePath + 'BulgeMass.npy')
BulgeSpin = np.load(SavePath + 'BulgeSpin.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')

# Trim the data #
index = np.where((np.linalg.norm(DiskSpin, axis=1) > 0.0) & (BulgeMass < 0.8 * StellarMass))

Disk_Mass = DiskMass[index] * MassUnits
Bulge_Mass = BulgeMass[index] * MassUnits
Stellar_Mass = StellarMass[index] * MassUnits
Disk_Spin = np.linalg.norm(DiskSpin[index], axis=1) * SpinUnits
Bulge_Spin = np.linalg.norm(BulgeSpin[index], axis=1) * SpinUnits

Mass_Ratio = np.divide(Bulge_Mass, Stellar_Mass)
Total_Spin = np.divide((Disk_Spin * Disk_Mass + Bulge_Spin * Bulge_Mass), Stellar_Mass)

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure = plt.figure(0, figsize=(10, 7.5))

# Figure parameters #
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e0, 1e5)
plt.xlim(1e9, 1e12)
plt.xlabel(r'$\mathrm{M_{\bigstar} / M_{\odot}}$')
plt.ylabel(r'$\mathrm{|j_{\bigstar}| / (km \cdot s^{-1} \cdot kpc)}$')
plt.tick_params(direction='in', which='both', top='on', right='on')

######################################################################################################################################################

# Read observational data from FR18 #
FR18 = np.genfromtxt('./Obs_Data/FR18.csv', delimiter=',', names=['Mstar', 'jstar'])

# Plot observational data from FR13 and OG14 #
plt.plot(np.power(10, FR18['Mstar'][0:2]), np.power(10, FR18['jstar'][0:2]), color='blue', lw=lw, linestyle='dashed',
         label=r'$\mathrm{Fall\; &\; Romanowsky\, 18:Discs}$', zorder=4)
plt.plot(np.power(10, FR18['Mstar'][2:4]), np.power(10, FR18['jstar'][2:4]), color='red', lw=lw, linestyle='dashed',
         label=r'$\mathrm{Fall\; &\; Romanowsky\, 18:Bulges}$', zorder=4)

# Create the legends #
colors = ['darkred', 'gold', 'darkblue']
circles = collections.CircleCollection([10] * 3, facecolor=colors)
legend = plt.legend([circles], [r'$\mathrm{This\; work:M_{b} / M_{\bigstar} < 0.8}$'], scatterpoints=len(colors), scatteryoffsets=[.5],
                    handlelength=len(colors), markerscale=2, frameon=False, loc=4)

plt.gca().add_artist(legend)
plt.legend(frameon=False, ncol=1, loc=2)

# Plot L-Galaxies data - color coded #
scatter = plt.scatter(Stellar_Mass, Total_Spin, c=Mass_Ratio, s=10, cmap='RdYlBu_r')
cbaxes = figure.add_axes([0.9, 0.11, 0.02, 0.77])
cb = plt.colorbar(scatter, cax=cbaxes)
cb.set_label(r'$\mathrm{M_{b} / M_{\bigstar}}$')

######################################################################################################################################################

# Save the figure #
plt.savefig('SM_Vs_GS_58-' + date + '.png', bbox_inches='tight')