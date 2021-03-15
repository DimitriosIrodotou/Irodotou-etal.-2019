# Load data arrays #
DiskMass = np.load(SavePath + 'DiskMass.npy')
BulgeMass = np.load(SavePath + 'BulgeMass.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')

SavePath = '/Volumes/BAM-BLACK/output/output_ITH_Off/58/'
DiskMassHWT15 = np.load(SavePath + 'DiskMass.npy')
BulgeMassHWT15 = np.load(SavePath + 'BulgeMass.npy')
StellarMassHWT15 = np.load(SavePath + 'StellarMass.npy')

######################################################################################################################################################

# Bins for histogram and plotting
firstFile = 0
lastFile = 511
maxFile = 512
binperdex = 5
xrange = np.array([7.8, 12.0])
nbin = (xrange[1] - xrange[0]) * binperdex

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(StellarMass * 1e10 * hubble), bins=np.int(nbin), range=xrange, log=True)
y = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
x = 0.5 * (bins[:-1] + bins[1:])

# Do the same for the HWT15 data #
# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(StellarMassHWT15 * 1e10 * hubble), bins=np.int(nbin), range=xrange, log=True)
yHen = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xHen = 0.5 * (bins[:-1] + bins[1:])

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

plt.plot(x, np.log10(y), color='black')
plt.plot(xHen, np.log10(yHen), color='black', linestyle='dotted')


# Create the legends #
# Create an object to combine the HWT15 lines in the legend #
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], color='black')
        return [l1]


legend1 = plt.legend([object], [r'$\mathrm{This\; work}$'], handler_map={object:AnyObjectHandler()}, frameon=False,
                     loc=2)


class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], linestyle='dotted', color='black')
        return [l1]


legend2 = plt.legend([object], [r'$\mathrm{HWT15}$'], handler_map={object:AnyObjectHandler()}, frameon=False, loc=1)

plt.gca().add_artist(legend1)
plt.gca().add_artist(legend2)
plt.legend(frameon=False, loc=3)

######################################################################################################################################################

# Save the figure #
plt.savefig('SMF_58-' + date + '.pdf', bbox_inches='tight')
