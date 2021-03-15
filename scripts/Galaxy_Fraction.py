# Load data arrays #
DiskMass = np.load(SavePath + 'DiskMass.npy')
BulgeMass = np.load(SavePath + 'BulgeMass.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')

SavePath = '/Volumes/BAM-BLACK/output/output_ITH_Off/58/'
DiskMassHWT15 = np.load(SavePath + 'DiskMass.npy')
BulgeMassHWT15 = np.load(SavePath + 'BulgeMass.npy')
StellarMassHWT15 = np.load(SavePath + 'StellarMass.npy')

# Put into observer units and add scatter to stellar mass estimate #
offset = 10 + np.log10(hubble)

log10Disk = np.log10(DiskMass) + offset
log10Bulge = np.log10(BulgeMass) + offset
log10StellarMass = np.log10(StellarMass) + offset
# random.seed(datetime.now())
log10StellarMassObs = log10StellarMass + np.random.randn(len(StellarMass)) * 0.08 * (1 + redshift)

log10DiskHWT15 = np.log10(DiskMassHWT15) + offset
log10BulgeHWT15 = np.log10(BulgeMassHWT15) + offset
log10StellarMassHWT15 = np.log10(StellarMassHWT15) + offset
log10StellarMassObsHWT15 = log10StellarMassHWT15 + np.random.randn(len(StellarMassHWT15)) * 0.08 * (1 + redshift)

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure, ax = plt.subplots()
figure, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20, 7.5))
figure.subplots_adjust(hspace=0, wspace=0)

# Figure parameters #
ax1.set_xlim(7.8, 11.9)
ax2.set_xlim(7.8, 11.9)
ax1.set_ylim(-0.1, 1.1)
ax2.set_ylim(-0.1, 1.1)

ax1.set_ylabel(r'$\mathrm{Fraction}$', size=22)
ax2.set_ylabel(r'$\mathrm{Fraction}$', size=22)
ax1.set_xlabel(r'$\mathrm{log_{10}(M_{\bigstar} / h^{-2}M_\odot)}$', size=22)
ax2.set_xlabel(r'$\mathrm{log_{10}(M_{\bigstar} / h^{-2}M_\odot)}$', size=22)

ax2.yaxis.set_label_position("right")

ax2.yaxis.tick_right()
ax1.tick_params(direction='in', which='both', top='on', right='on', labelsize=22)
ax2.tick_params(direction='in', which='both', top='on', left='on', labelsize=22)

######################################################################################################################################################

# Read observational data from C06 and ﻿KDR14 #
obsBulge = np.loadtxt('./Obs_Data/Conselice06_Bulge_Frac.txt')
obsDisk = np.loadtxt('./Obs_Data/Conselice06_Disk_Frac.txt')
obsIrr = np.loadtxt('./Obs_Data/Conselice06_Irr_Frac.txt')

x = [9, 9.5, 10, 10.5, 11, 11.5]
yE = [0.23, 0.25, 0.24, 0.31, 0.72, 1.00]
yS0Scd = [0.46, 0.64, 0.74, 0.68, 0.28, 0.00]
ySd = [0.31, 0.11, 0.02, 0.00, 0.00, 0.00]

# Plot observational data from C06 and ﻿KDR14 #
obsMass = obsBulge[:, 0] + 2 * np.log10(obsHubble)
ax1.errorbar(obsMass, obsBulge[:, 1], yerr=obsBulge[:, 2], marker='o', color='red', linestyle='None', elinewidth=1,
             capsize=3, capthick=1)

obsMass = obsDisk[:, 0] + 2 * np.log10(obsHubble)
ax1.errorbar(obsMass, obsDisk[:, 1], yerr=obsDisk[:, 2], marker='o', color='green', linestyle='None', elinewidth=1,
             capsize=3, capthick=1)

obsMass = obsIrr[:, 0] + 2 * np.log10(obsHubble)
ax1.errorbar(obsMass, obsIrr[:, 1], yerr=obsIrr[:, 2], marker='o', color='blue', linestyle='None', elinewidth=1,
             capsize=3, capthick=1)

ax1.errorbar(x, yE, xerr=0.1, color='red', marker='s', linestyle='None', elinewidth=1, capsize=3, capthick=1)
ax1.errorbar(x, yS0Scd, xerr=0.1, color='green', marker='s', linestyle='None', elinewidth=1, capsize=3, capthick=1)
ax1.errorbar(x, ySd, xerr=0.1, color='blue', marker='s', linestyle='None', elinewidth=1, capsize=3, capthick=1)

# Divisions between bulge classes #
logRatio1 = -0.154902
logRatio2 = -2

# Bins for histogram and plotting #
binwidth = 0.25
xrange = np.array([8.67, 11.8])
bins = np.arange(xrange[0], xrange[1] + 0.001, binwidth)

# Put galaxies into bins #
indBin = np.digitize(log10StellarMassObs, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yIrr = np.empty(nBin)
yDisk = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((log10Bulge[indThisBin] - log10StellarMass[indThisBin]) > logRatio1)[0]) / float(allBin)
    # Disks
    yIrr[iBin] = len(np.where((log10Bulge[indThisBin] - log10StellarMass[indThisBin]) < logRatio2)[0]) / float(allBin)
    # Intermediates
    yDisk[iBin] = 1. - yBulge[iBin] - yIrr[iBin]

# Plot L-Galaxies data #
ax1.plot(x, yIrr, color='blue', lw=2)
ax1.plot(x, yBulge, color='red', lw=2)
ax1.plot(x, yDisk, color='green', lw=2)

# Put HWT15 galaxies into bins #
indBin = np.digitize(log10StellarMassObsHWT15, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yIrr = np.empty(nBin)
yDisk = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(
        np.where((log10BulgeHWT15[indThisBin] - log10StellarMassHWT15[indThisBin]) > logRatio1)[0]) / float(allBin)
    # Disks
    yIrr[iBin] = len(
        np.where((log10BulgeHWT15[indThisBin] - log10StellarMassHWT15[indThisBin]) < logRatio2)[0]) / float(allBin)
    # Intermediates
    yDisk[iBin] = 1. - yBulge[iBin] - yIrr[iBin]

# Plot L-Galaxies data #
ax1.plot(x, yIrr, color='blue', lw=2, linestyle='dotted')
ax1.plot(x, yBulge, color='red', lw=2, linestyle='dotted')
ax1.plot(x, yDisk, color='green', lw=2, linestyle='dotted')

# Bins for histogram and plotting #
binwidth = 0.25
xrange = np.array([7.8, 9.05])
bins = np.arange(xrange[0], xrange[1] + 0.001, binwidth)

# Put galaxies into bins #
indBin = np.digitize(log10StellarMassObs, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yIrr = np.empty(nBin)
yDisk = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((log10Bulge[indThisBin] - log10StellarMass[indThisBin]) > logRatio1)[0]) / float(allBin)
    # Disks
    yIrr[iBin] = len(np.where((log10Bulge[indThisBin] - log10StellarMass[indThisBin]) < logRatio2)[0]) / float(allBin)
    # Intermediates
    yDisk[iBin] = 1. - yBulge[iBin] - yIrr[iBin]

# Plot L-Galaxies data #
# ax1.plot(x, yIrr, color='blue', alpha=0.3, lw=2)
# ax1.plot(x, yBulge, color='red', alpha=0.3, lw=2)
# ax1.plot(x, yDisk, color='green', alpha=0.3, lw=2)

# Put HWT15 galaxies into bins #
indBin = np.digitize(log10StellarMassObsHWT15, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yIrr = np.empty(nBin)
yDisk = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(
        np.where((log10BulgeHWT15[indThisBin] - log10StellarMassHWT15[indThisBin]) > logRatio1)[0]) / float(allBin)
    # Disks
    yIrr[iBin] = len(
        np.where((log10BulgeHWT15[indThisBin] - log10StellarMassHWT15[indThisBin]) < logRatio2)[0]) / float(allBin)
    # Intermediates
    yDisk[iBin] = 1. - yBulge[iBin] - yIrr[iBin]


# Plot L-Galaxies data #
# ax1.plot(x, yIrr, color='blue', lw=2, alpha=0.3, linestyle='dotted')
# ax1.plot(x, yBulge, color='red', lw=2, alpha=0.3, linestyle='dotted')
# ax1.plot(x, yDisk, color='green', lw=2, alpha=0.3, linestyle='dotted')


# Create the legends #
# Create an object to combine the HWT15 lines in the legend #
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [1.0 * height, 1.0 * height], color='blue')
        l2 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], color='green')
        l3 = plt.Line2D([x0, y0 + width], [0.0 * height, 0.0 * height], color='red')
        return [l1, l2, l3]


legend1 = ax1.legend([object], [r'$\mathrm{This\; work}$'], handler_map={object:AnyObjectHandler()}, frameon=False,
                     loc=2)


class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [1.0 * height, 1.0 * height], linestyle='dotted', color='blue')
        l2 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], linestyle='dotted', color='green')
        l3 = plt.Line2D([x0, y0 + width], [0.0 * height, 0.0 * height], linestyle='dotted', color='red')
        return [l1, l2, l3]


legend2 = ax1.legend([object], [r'$\mathrm{HWT15}$'], handler_map={object:AnyObjectHandler()}, frameon=False, loc=1)

colors = ['blue', 'green', 'red']
circles = collections.CircleCollection([10] * 3, facecolor=colors)
colors = ['blue', 'green', 'red']
squares = collections.RegularPolyCollection(numsides=4, rotation=0.785, sizes=(20,), facecolors=colors)

legend3 = ax1.legend([squares, circles], [r'$\mathrm{Kelvin+14}$', r'$\mathrm{Conselice\, 06}$'],
                     scatterpoints=len(colors), scatteryoffsets=[.5], handlelength=len(colors), markerscale=2,
                     frameon=False, loc=9)

ax1.add_artist(legend1)
ax1.add_artist(legend2)
ax1.add_artist(legend3)

# Add text annotation #
ax1.annotate(r'$\mathrm{0.01 < M_{b} / M_{\bigstar} < 0.7}$', xy=(0.01, 0.02), xycoords='axes fraction', color='green',
             size=15)
ax1.annotate(r'$\mathrm{M_{b} / M_{\bigstar} > 0.7}$', xy=(0.46, 0.02), xycoords='axes fraction', color='red', size=15)
ax1.annotate(r'$\mathrm{M_{b} / M_{\bigstar} < 0.01}$', xy=(0.81, 0.02), xycoords='axes fraction', color='blue',
             size=15)

######################################################################################################################################################

# Read observational data from MID16 #
MID16SphDom = np.genfromtxt('./Obs_Data/MID16SphDom.csv', delimiter=',', names=['StellarMass', 'Fraction'])
MID16SphDomTop = np.genfromtxt('./Obs_Data/MID16SphDomTop.csv', delimiter=',', names=['StellarMass', 'Fraction'])
MID16SphDomBot = np.genfromtxt('./Obs_Data/MID16SphDomBot.csv', delimiter=',', names=['StellarMass', 'Fraction'])

MID16DiskDom = np.genfromtxt('./Obs_Data/MID16DiskDom.csv', delimiter=',', names=['StellarMass', 'Fraction'])
MID16DiskDomTop = np.genfromtxt('./Obs_Data/MID16DiskDomTop.csv', delimiter=',', names=['StellarMass', 'Fraction'])
MID16DiskDomBot = np.genfromtxt('./Obs_Data/MID16DiskDomBot.csv', delimiter=',', names=['StellarMass', 'Fraction'])

# Plot observational data from MID16 #
median, = ax2.plot(MID16SphDom['StellarMass'], MID16SphDom['Fraction'], linestyle="dashed", color='red')
ax2.fill_between(MID16SphDomTop['StellarMass'], MID16SphDomTop['Fraction'], MID16SphDomBot['Fraction'], color='red',
                 alpha='0.5', zorder=2)
fill, = plt.fill(np.NaN, np.NaN, color='red', alpha=0.5)

median1, = ax2.plot(MID16DiskDom['StellarMass'], MID16DiskDom['Fraction'], linestyle="dashed", color='blue')
ax2.fill_between(MID16DiskDomTop['StellarMass'], MID16DiskDomTop['Fraction'], MID16DiskDomBot['Fraction'], color='blue',
                 alpha='0.5', zorder=2)
fill1, = plt.fill(np.NaN, np.NaN, color='blue', alpha=0.5)

# Divisions between bulge classes #
logRatio1 = -0.523
logRatio2 = -0.523

# Bins for histogram and plotting #
binwidth = 0.25
xrange = np.array([8.67, 11.8])
bins = np.arange(xrange[0], xrange[1] + 0.001, binwidth)

# Put galaxies into bins #
indBin = np.digitize(log10StellarMassObs, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yIrr = np.empty(nBin)
yDisk = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((log10Bulge[indThisBin] - log10StellarMass[indThisBin]) > logRatio1)[0]) / float(allBin)
    # Disks
    yDisk[iBin] = len(np.where((log10Bulge[indThisBin] - log10StellarMass[indThisBin]) < logRatio2)[0]) / float(allBin)

# Plot L-Galaxies data #
ax2.plot(x, yDisk, color='blue', lw=2)
ax2.plot(x, yBulge, color='red', lw=2)

# Put HWT15 galaxies into bins #
indBin = np.digitize(log10StellarMassObsHWT15, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yIrr = np.empty(nBin)
yDisk = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(
        np.where((log10BulgeHWT15[indThisBin] - log10StellarMassHWT15[indThisBin]) > logRatio1)[0]) / float(allBin)
    # Disks
    yDisk[iBin] = len(
        np.where((log10BulgeHWT15[indThisBin] - log10StellarMassHWT15[indThisBin]) < logRatio2)[0]) / float(allBin)

# Plot L-Galaxies data #
ax2.plot(x, yDisk, color='blue', lw=2, linestyle='dotted')
ax2.plot(x, yBulge, color='red', lw=2, linestyle='dotted')

# Bins for histogram and plotting #
binwidth = 0.25
xrange = np.array([7.8, 9.05])
bins = np.arange(xrange[0], xrange[1] + 0.001, binwidth)

# Put galaxies into bins #
indBin = np.digitize(log10StellarMassObs, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yIrr = np.empty(nBin)
yDisk = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((log10Bulge[indThisBin] - log10StellarMass[indThisBin]) > logRatio1)[0]) / float(allBin)
    # Disks
    yDisk[iBin] = len(np.where((log10Bulge[indThisBin] - log10StellarMass[indThisBin]) < logRatio2)[0]) / float(allBin)

# Plot L-Galaxies data #
# ax2.plot(x, yDisk, color='blue', alpha=0.3, lw=2)
# ax2.plot(x, yBulge, color='red', alpha=0.3, lw=2)

# Put HWT15 galaxies into bins #
indBin = np.digitize(log10StellarMassObsHWT15, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yIrr = np.empty(nBin)
yDisk = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(
        np.where((log10BulgeHWT15[indThisBin] - log10StellarMassHWT15[indThisBin]) > logRatio1)[0]) / float(allBin)
    # Disks
    yDisk[iBin] = len(
        np.where((log10BulgeHWT15[indThisBin] - log10StellarMassHWT15[indThisBin]) < logRatio2)[0]) / float(allBin)


# Plot L-Galaxies data #
# ax2.plot(x, yDisk, color='blue', lw=2, alpha=0.3, linestyle='dotted')
# ax2.plot(x, yBulge, color='red', lw=2, alpha=0.3, linestyle='dotted')


# Create the legends #
# Create an object to combine the HWT15 lines in the legend #
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.9 * height, 0.9 * height], color='blue')
        l2 = plt.Line2D([x0, y0 + width], [0.1 * height, 0.1 * height], color='red')
        return [l1, l2]


legend1 = ax2.legend([object, object], [r'$\mathrm{This\; work}$'], handler_map={object:AnyObjectHandler()},
                     frameon=False, loc=2)


class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.9 * height, 0.9 * height], linestyle='dotted', color='blue')
        l2 = plt.Line2D([x0, y0 + width], [0.1 * height, 0.1 * height], linestyle='dotted', color='red')
        return [l1, l2]


legend2 = ax2.legend([object, object], [r'$\mathrm{HWT15}$'], handler_map={object:AnyObjectHandler()}, frameon=False,
                     loc=1)

legend3 = plt.legend([median1, fill1, median, fill],
                     [r'$\mathrm{Disc-Dom}$', r'$\mathrm{Error\,ranges}$', r'$\mathrm{Sph-Dom}$',
                      r'$\mathrm{Error\,ranges}$'], loc=6)

ax2.add_artist(legend1)
ax2.add_artist(legend2)
ax2.add_artist(legend3)

# Add text annotation #
ax2.annotate(r'$\mathrm{Moffett+16:}$', xy=(0.022, 0.645), xycoords='axes fraction', size=15)
ax2.annotate(r'$\mathrm{M_{b} / M_{\bigstar} > 0.3}$', xy=(0.01, 0.02), xycoords='axes fraction', color='red', size=15)
ax2.annotate(r'$\mathrm{M_{b} / M_{\bigstar} < 0.3}$', xy=(0.83, 0.02), xycoords='axes fraction', color='blue', size=15)

######################################################################################################################################################

# Save the figure #
plt.savefig('Gal_Frac_58-' + date + '.pdf', bbox_inches='tight')
