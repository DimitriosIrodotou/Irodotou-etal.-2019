# Load data arrays #
BulgeMass = np.load(SavePath + 'BulgeMass.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')

SavePath = '/Volumes/BAM-BLACK/output/output_ITH_Off/58/'
BulgeMassHWT15 = np.load(SavePath + 'BulgeMass.npy')
StellarMassHWT15 = np.load(SavePath + 'StellarMass.npy')

indexDD = np.where(BulgeMass < 0.3 * StellarMass)
indexBD = np.where(BulgeMass > 0.3 * StellarMass)
indexIrr = np.where(BulgeMass == 0.0)

indexDDHWT15 = np.where(BulgeMassHWT15 < 0.3 * StellarMassHWT15)
indexBDHWT15 = np.where(BulgeMassHWT15 > 0.3 * StellarMassHWT15)
indexIrrHWT15 = np.where(BulgeMassHWT15 == 0.0)

# Bins for histogram and plotting
firstFile = 0
lastFile = 511
maxFile = 512
binperdex = 5
xrange = np.array([7.8, 12.0])
nbin = (xrange[1] - xrange[0]) * binperdex

######################################################################################################################################################

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(StellarMass[indexDD] * 1e10 * hubble), bins=np.int(nbin), range=xrange, log=True)
yDD = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xDD = 0.5 * (bins[:-1] + bins[1:])

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(StellarMass[indexBD] * 1e10 * hubble), bins=np.int(nbin), range=xrange, log=True)
yBD = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xBD = 0.5 * (bins[:-1] + bins[1:])

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(StellarMass[indexIrr] * 1e10 * hubble), bins=np.int(nbin), range=xrange, log=True)
yIrr = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xIrr = 0.5 * (bins[:-1] + bins[1:])

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(StellarMassHWT15[indexDDHWT15] * 1e10 * hubble), bins=np.int(nbin), range=xrange,
                            log=True)
yDDHen = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xDDHen = 0.5 * (bins[:-1] + bins[1:])

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(StellarMassHWT15[indexBDHWT15] * 1e10 * hubble), bins=np.int(nbin), range=xrange,
                            log=True)
yBDHen = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xBDHen = 0.5 * (bins[:-1] + bins[1:])

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(StellarMassHWT15[indexIrrHWT15] * 1e10 * hubble), bins=np.int(nbin), range=xrange,
                            log=True)
yIrrHen = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xIrrHen = 0.5 * (bins[:-1] + bins[1:])

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure = plt.figure(0, figsize=(10, 7.5))

# Figure parameters #
plt.ylim(-5.2, -0.2)
plt.xlim(8.8, 12.2)

plt.ylabel(r'$\mathrm{log_{10}(N / (dex / (h^{-1}Mpc)^{3}))}$')
plt.xlabel(r'$\mathrm{log_{10}(M_{\bigstar} / h^{-2}M_\odot)}$')

plt.tick_params(direction='in', which='both', top='on', right='on')

######################################################################################################################################################

# Read MCMC data #
HWT15 = np.genfromtxt('./Obs_Data/Henriques2015a_smf_z01.csv', delimiter=',', names=['x1', 'x2', 'y', 'err'])

# Plot MCMC data #
ytop = np.log10(HWT15['y'] + HWT15['err']) - np.log10(HWT15['y'])
ybot = np.log10(HWT15['y']) - np.log10(HWT15['y'] - HWT15['err'])
plt.errorbar((HWT15['x2'] + HWT15['x1']) / 2, np.log10(HWT15['y']), color='black', yerr=(ybot, ytop), marker='o',
             linestyle="None", elinewidth=1, capsize=4, capthick=1, zorder=4,
             label=r'$\mathrm{Observations\, used\, in\, MCMC}$')

# Read observational data from MID16 #
MID16D = np.genfromtxt('./Obs_Data/MID16D.csv', delimiter=',', names=['x', 'y'])
MID16Sph = np.genfromtxt('./Obs_Data/MID16Sph.csv', delimiter=',', names=['x', 'y'])

# Plot observational data from MID16 #
ytop = np.array(np.log10(MID16Sph['y'][0::3].copy()) - 3 * np.log10(obsHubble)) - np.array(
    np.log10(MID16Sph['y'][1::3].copy()) - 3 * np.log10(obsHubble))

ybot = np.array(np.log10(MID16Sph['y'][1::3].copy()) - 3 * np.log10(obsHubble)) - np.array(
    np.log10(MID16Sph['y'][2::3].copy()) - 3 * np.log10(obsHubble))

plt.errorbar(MID16Sph['x'][1::3] + 2 * np.log10(obsHubble), np.log10(MID16Sph['y'][1::3]) - 3 * np.log10(obsHubble),
             color='red', yerr=(ybot, ytop), marker='o', linestyle="None", elinewidth=1, capsize=4, capthick=1,
             zorder=3)

ytop = np.array(np.log10(MID16D['y'][0::3].copy()) - 3 * np.log10(obsHubble)) - np.array(
    np.log10(MID16D['y'][1::3].copy()) - 3 * np.log10(obsHubble))

ybot = np.array(np.log10(MID16D['y'][1::3].copy()) - 3 * np.log10(obsHubble)) - np.array(
    np.log10(MID16D['y'][2::3].copy()) - 3 * np.log10(obsHubble))

plt.errorbar(MID16D['x'][1::3] + 2 * np.log10(obsHubble), np.log10(MID16D['y'][1::3]) - 3 * np.log10(obsHubble),
             color='blue', yerr=(ybot, ytop), marker='o', linestyle="None", elinewidth=1, capsize=4, capthick=1,
             zorder=3)

# Read observational data from MLD16 #
MLD16 = np.loadtxt('./Obs_Data/MLD16.txt', skiprows=2)

# Plot observational data from MLD16 #
ytop = np.array(np.log10(MLD16[0:22, 21]) - 3 * np.log10(obsHubble)) - np.array(
    np.log10(MLD16[0:22, 20]) - 3 * np.log10(obsHubble))
ybot = np.array(np.log10(MLD16[0:22, 20]) - 3 * np.log10(obsHubble)) - np.array(
    np.log10(MLD16[0:22, 19]) - 3 * np.log10(obsHubble))
plt.errorbar(MLD16[0:22, 0] + 2 * np.log10(obsHubble), np.log10(MLD16[0:22, 20]) - 3 * np.log10(obsHubble),
             color='darkblue', yerr=(ybot, ytop), marker='o', linestyle="None", elinewidth=1, capsize=4, capthick=1,
             zorder=3)

plt.plot(xBD, np.log10(yBD), color='red')
plt.plot(xDD, np.log10(yDD), color='blue')
plt.plot(xIrr, np.log10(yIrr), color='darkblue')

plt.plot(xBDHen, np.log10(yBDHen), linestyle='dotted', color='red')
plt.plot(xDDHen, np.log10(yDDHen), linestyle='dotted', color='blue')
plt.plot(xIrrHen, np.log10(yIrrHen), linestyle='dotted', color='darkblue')


# Create the legends #
# Create an object to combine the HWT15 lines in the legend #
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [1.0 * height, 1.0 * height], color='red')
        l2 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], color='blue')
        l3 = plt.Line2D([x0, y0 + width], [0.0 * height, 0.0 * height], color='darkblue')
        return [l1, l2, l3]


legend1 = plt.legend([object], [r'$\mathrm{This\; work}$'], handler_map={object:AnyObjectHandler()}, frameon=False,
                     loc=2)


class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [1.0 * height, 1.0 * height], linestyle='dotted', color='red')
        l2 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], linestyle='dotted', color='blue')
        l3 = plt.Line2D([x0, y0 + width], [0.0 * height, 0.0 * height], linestyle='dotted', color='darkblue')
        return [l1, l2, l3]


legend2 = plt.legend([object], [r'$\mathrm{HWT15}$'], handler_map={object:AnyObjectHandler()}, frameon=False, loc=1)

colors = ['red', 'blue', 'darkblue']
circles = collections.CircleCollection([10] * 3, facecolor=colors)

legend3 = plt.legend([circles], [r'$\mathrm{Moffett+16\, a, b}$'], scatterpoints=len(colors), scatteryoffsets=[.5],
                     handlelength=len(colors), markerscale=2, frameon=False, loc=9)

plt.gca().add_artist(legend1)
plt.gca().add_artist(legend2)
plt.gca().add_artist(legend3)
plt.legend(frameon=False, loc=(0.015, 0.815))

# Add text annotation #
plt.annotate(r'$\mathrm{M_{b} / M_{\bigstar} < 0.3}$', xy=(0.5, 0.02), xycoords='axes fraction', color='blue', size=15)
plt.annotate(r'$\mathrm{M_{b} / M_{\bigstar} = 0.0}$', xy=(0.2, 0.2), xycoords='axes fraction', color='darkblue',
             size=15)
plt.annotate(r'$\mathrm{M_{b} / M_{\bigstar} > 0.3}$', xy=(0.81, 0.12), xycoords='axes fraction', color='red', size=15)

######################################################################################################################################################

# Save the figure #
plt.savefig('SMF_Morph_58-' + date + '.pdf', bbox_inches='tight')
