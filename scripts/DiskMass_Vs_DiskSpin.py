# Load data arrays #
DiskMass = np.load(SavePath + 'DiskMass.npy')
DiskSpin = np.load(SavePath + 'DiskSpin.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')

SavePath = '/Volumes/BAM-BLACK/output/output_ITH_Off/58/'
DiskMassHWT15 = np.load(SavePath + 'DiskMass.npy')
DiskSpinHWT15 = np.load(SavePath + 'DiskSpin.npy')
StellarMassHWT15 = np.load(SavePath + 'StellarMass.npy')

# Trim the data #
index = np.where(DiskMass > 0.7 * StellarMass)
indexHWT15 = np.where(DiskMassHWT15 > 0.7 * StellarMassHWT15)

Disk_Mass = DiskMass[index] * MassUnits
Disk_Spin = np.linalg.norm(DiskSpin[index], axis=1) * SpinUnits
Disk_MassHWT15 = DiskMassHWT15[indexHWT15] * MassUnits
Disk_SpinHWT15 = np.linalg.norm(DiskSpinHWT15[indexHWT15], axis=1) * SpinUnits

dlog10 = 0.1

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20, 7.5))
figure.subplots_adjust(hspace=0, wspace=0.37)

# Figure parameters #
ax1.set_xscale('log')
ax2.set_xscale('log')
ax1.set_yscale('log')
ax2.set_yscale('log')

ax1.set_ylim(1e1, 1e5)
ax2.set_ylim(1e1, 1e5)
ax1.set_xlim(1e9, 1e12)
ax2.set_xlim(1e9, 1e12)

ax1.set_xlabel(r'$\mathrm{M_{d,\bigstar} / M_{\odot}}$', size=25)
ax2.set_xlabel(r'$\mathrm{M_{d,\bigstar} / M_{\odot}}$', size=25)
ax1.set_ylabel(r'$\mathrm{|j_{d,\bigstar}| / (km \cdot s^{-1} \cdot kpc)}$', size=25)
ax2.set_ylabel(r'$\mathrm{|j_{d,\bigstar}| / (km \cdot s^{-1} \cdot kpc)}$', size=25)

ax1.set_xticklabels([])
ax1.tick_params(direction='in', which='both', top='on', right='on', labelsize=25)
ax2.tick_params(direction='in', which='both', top='on', right='on', labelsize=25)

######################################################################################################################################################

# Read observational data from FR13 and OG14 #
FR13 = np.genfromtxt('./Obs_Data/FR13.csv', delimiter=',', names=['Md', 'jd'])
x = [10 ** 8.407256236, 10 ** 12.00269017]
y = [10 ** 1.995464853, 10 ** 3.986394558]
OG14 = np.genfromtxt('./Obs_Data/OG14.csv', delimiter=',', names=['Md', 'jd'])

# Plot observational data from FR13 and OG14 #
ax1.scatter(np.power(10, FR13['Md']), np.power(10, FR13['jd']), edgecolor='black', color='blue', s=50, marker='s',
            label=r'$\mathrm{Fall\; &\; Romanowsky\, 13}$', zorder=4)
ax1.scatter(np.power(10, OG14['Md']), np.power(10, OG14['jd']), edgecolor='black', color='lime', s=50, marker='^',
            label=r'$\mathrm{Obreschkow\; &\; Glazebrook\, 14}$', zorder=4)

# Plot L-Galaxies data - 2D histogram #
hexbin = ax1.hexbin(Disk_Mass, Disk_Spin, xscale='log', yscale='log', bins='log', cmap='Greys', mincnt=2)

# Adjust the color bar #
cbaxes = figure.add_axes([0.452, 0.11, 0.01, 0.77])
cb = plt.colorbar(hexbin, cax=cbaxes)
cbaxes.tick_params(direction='out', which='both', right='on', labelsize=25)
cb.set_label('$\mathrm{Counts\; per\; hexbin}$', size=25)

# Create the legends #
colors = ['black', 'grey', 'lightgrey']
squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
legend1 = ax1.legend([squares], [r'$\mathrm{This\;work:M_{d,\bigstar} / M_{\bigstar}> 0.7}$'], scatterpoints=len(colors), scatteryoffsets=[.5],
                     handlelength=len(colors), markerscale=2, frameon=False, loc=2)

ax1.add_artist(legend1)
ax1.legend(frameon=False, ncol=1, loc=4, scatterpoints=3)

######################################################################################################################################################

# Calculate median and 1-sigma #
log10X = np.log10(Disk_Mass)
log10XMax = np.log10(max(Disk_Mass))
log10XMin = np.log10(min(Disk_Mass))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Disk_Mass[index])
    if len(index) > 0:
        median[i] = np.median(Disk_Spin[index])
        slow[i] = np.percentile(Disk_Spin[index], 15.87)
        shigh[i] = np.percentile(Disk_Spin[index], 84.13)
    log10XLow += dlog10
# Plot median and 1-sigma lines #
median, = ax2.plot(X, median, color='black', lw=lw)
ax2.fill_between(X, shigh, slow, color='black', alpha='0.5', zorder=2)
fill, = plt.fill(np.NaN, np.NaN, color='black', alpha=0.5)

legend2 = ax2.legend([median, fill], [r'$\mathrm{This\; work: Median}$', r'$\mathrm{This\; work:16^{th}-84^{th}\,\%ile}$'], scatterpoints=len(colors),
                     scatteryoffsets=[.5], markerscale=2, frameon=False, loc=2)

log10X = np.log10(Disk_MassHWT15)
log10XMax = np.log10(max(Disk_MassHWT15))
log10XMin = np.log10(min(Disk_MassHWT15))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Disk_MassHWT15[index])
    if len(index) > 0:
        median[i] = np.median(Disk_SpinHWT15[index])
        slow[i] = np.percentile(Disk_SpinHWT15[index], 15.87)
        shigh[i] = np.percentile(Disk_SpinHWT15[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = ax2.plot(X, median, color='black', lw=lw, linestyle='dotted')

# Plot median the theoretical relation #
line1, = ax2.plot(x, y, color='blue', lw=lw, linestyle='dashed', label=r'$\mathrm{j_{\bigstar} \propto M_{\bigstar}^{0.6}}$')

# Create the legends #
legend3 = ax2.legend([median], [r'$\mathrm{HWT15: Median}$'], frameon=False, loc=1)

ax2.add_artist(legend2)
ax2.add_artist(legend3)
ax2.legend(frameon=False, ncol=1, loc=4, scatterpoints=3)

######################################################################################################################################################

# Save the figure #
plt.savefig('DM_Vs_DS_58-' + date + '.pdf', bbox_inches='tight')