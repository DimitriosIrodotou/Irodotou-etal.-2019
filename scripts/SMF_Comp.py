# Import required python libraries #
import matplotlib
matplotlib.use('Agg')
import time
import h5py
import random
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.optimize as so
from datetime import datetime
from scipy.ndimage import zoom
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors, ticker, cm
import matplotlib.collections as collections
from matplotlib.legend_handler import HandlerBase
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

# Declare the plotting style #
sns.set()
sns.set_style('ticks', {'axes.grid': True})
sns.set_context('notebook', font_scale=1.6)
lw = 3  # linewidth
mc = 2  # mincnt
sp = 3  # scatterpoints
ms = 5  # markersize
gs = 150  # gridsize
size = 50  # size

# Physical units and simulation parameters #
redshift = 0.0
hubble = 0.673  # [dimensionless]
obsHubble = 0.70  # [dimensionless]
boxside = 480.28  # [Mpc/h]
VelocityUnits = 1  # [Km/s]
LengthUnits = 1e3 / hubble  # [Kpc]
MassUnits = 1e10 / hubble  # [Msun]
SpinUnits = 1e3 / hubble  # [Km/s * Kpc]
date = time.strftime("%d\%m\%y\%H%M")

snap='58'
SavePath = '/lustre/scratch/astro/di43/Python/output/' + snap + '/'

# Load data arrays #
DiskMass = np.load(SavePath + 'DiskMass.npy')
BulgeMass = np.load(SavePath + 'BulgeMass.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')

DiskMassHen15 = np.load(SavePath + 'Output_All_Off/' + 'DiskMass.npy')
BulgeMassHen15 = np.load(SavePath + 'Output_All_Off/' + 'BulgeMass.npy')
StellarMassHen15 = np.load(SavePath + 'Output_All_Off/' + 'StellarMass.npy')

# Bins for histogram and plotting
firstFile = 0
lastFile = 511
maxFile = 512
binperdex = 5
xrange = np.array([7.8, 12.0])
nbin = (xrange[1] - xrange[0]) * binperdex

########################################################################################################################

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(StellarMass * 1e10 * hubble), bins=np.int(nbin), range=xrange, log=True)
y = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
x = 0.5 * (bins[:-1] + bins[1:])

# Do the same for the Hen15 data #
# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(StellarMassHen15 * 1e10 * hubble), bins=np.int(nbin), range=xrange, log=True)
yHen = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xHen = 0.5 * (bins[:-1] + bins[1:])

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(BulgeMass * 1e10 * hubble), bins=np.int(nbin), range=xrange, log=True)
yBM = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xBM = 0.5 * (bins[:-1] + bins[1:])

nobj, bins, junk = plt.hist(np.log10(BulgeMassHen15 * 1e10 * hubble), bins=np.int(nbin), range=xrange, log=True)
yBMHen = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xBMHen = 0.5 * (bins[:-1] + bins[1:])

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(DiskMass * 1e10 * hubble), bins=np.int(nbin), range=xrange, log=True)
yDM = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xDM = 0.5 * (bins[:-1] + bins[1:])

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj, bins, junk = plt.hist(np.log10(DiskMassHen15 * 1e10 * hubble), bins=np.int(nbin), range=xrange, log=True)
yDMHen = nobj * maxFile / ((lastFile - firstFile + 1) * boxside ** 3) * binperdex

# Plot at centre of bins
xDMHen = 0.5 * (bins[:-1] + bins[1:])

########################################################################################################################

# Generate initial figure #
plt.close()
figure, ax = plt.subplots()
figure, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(20, 15))
figure.subplots_adjust(hspace=0, wspace=0)

# Plot parameters #
ax2.set_ylim(-5.2, -0.2)
ax4.set_ylim(-5.2, -0.2)
ax6.set_ylim(-5.2, -0.2)
ax2.set_xlim(8.8, 12.2)
ax4.set_xlim(8.8, 12.2)
ax6.set_xlim(8.8, 12.2)

ax2.set_ylabel(r'$\mathrm{log_{10}(N / (dex / (h^{-1}Mpc)^{3})}$')
ax4.set_ylabel(r'$\mathrm{log_{10}(N / (dex / (h^{-1}Mpc)^{3})}$')
ax6.set_ylabel(r'$\mathrm{log_{10}(N / (dex / (h^{-1}Mpc)^{3})}$')
ax2.set_xlabel(r'$\mathrm{log_{10}(M_{\bigstar} / h^{-2}M_\odot)}$')
ax4.set_xlabel(r'$\mathrm{log_{10}(M_{b} / h^{-2}M_\odot)}$')
ax6.set_xlabel(r'$\mathrm{log_{10}(M_{\bigstar,d} / h^{-2}M_\odot)}$')

# Remove unwanted subplots #
ax1.axis('off')
ax3.axis('off')
ax5.axis('off')

ax2.tick_params(direction='in', which='both', top='on', right='on')
ax4.tick_params(direction='in', which='both', top='on', right='on')
ax6.tick_params(direction='in', which='both', top='on', right='on')

########################################################################################################################

# Read MCMC data #
HWT15 = np.genfromtxt('./Obs_Data/Henriques2015a_smf_z01.csv', delimiter=',', names=['x1', 'x2', 'y', 'err'])

# Plot MCMC data #
ytop = np.log10(HWT15['y'] + HWT15['err']) - np.log10(HWT15['y'])
ybot = np.log10(HWT15['y']) - np.log10(HWT15['y'] - HWT15['err'])
ax2.errorbar((HWT15['x2'] + HWT15['x1']) / 2, np.log10(HWT15['y']), color='black', yerr=(ybot, ytop), marker='o', linestyle="None", elinewidth=1,
             capsize=4, capthick=1, zorder=4, label=r'$\mathrm{Observations\, used\, in\, MCMC}$')

ax2.plot(x, np.log10(y), color='black')
ax2.plot(xHen, np.log10(yHen), color='black', linestyle='dotted')


# Create the legends #
# Create an object to combine the Hen15 lines in the legend #
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], color='black')
        return [l1]


legend1 = ax2.legend([object], [r'$\mathrm{This\; work}$'], handler_map={object: AnyObjectHandler()}, frameon=False, loc=2)


class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], linestyle='dotted', color='black')
        return [l1]


legend2 = ax2.legend([object], [r'$\mathrm{HWT15}$'], handler_map={object: AnyObjectHandler()}, frameon=False, loc=1)

ax2.add_artist(legend1)
ax2.add_artist(legend2)
ax2.legend(frameon=False, loc=3)

########################################################################################################################

# Read observational data from MLD16 #
MLD16 = np.loadtxt('./Obs_Data/MLD16.txt', skiprows=2)

# Plot L-Galaxies and MCMC data #
ax4.plot(xBM, np.log10(yBM), color='red')
ax4.plot(xBMHen, np.log10(yBMHen), linestyle='dotted', color='red')

# Plot observational data from #
ytop = np.array(np.log10(MLD16[0:34, 27]) - 3 * np.log10(obsHubble)) - np.array(np.log10(MLD16[0:34, 26]) - 3 * np.log10(obsHubble))
ybot = np.array(np.log10(MLD16[0:34, 26]) - 3 * np.log10(obsHubble)) - np.array(np.log10(MLD16[0:34, 25]) - 3 * np.log10(obsHubble))
ax4.errorbar(MLD16[0:34, 0] + 2 * np.log10(obsHubble), np.log10(MLD16[0:34, 26]) - 3 * np.log10(obsHubble), color='red', yerr=(ybot, ytop),
             marker='o', linestyle="None", elinewidth=1, capsize=4, capthick=1, zorder=3, label=r'$\mathrm{Moffett+16\,b}$')


# Create the legends #
# Create an object to combine the Hen15 lines in the legend #
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], color='red')
        return [l1]


legend3 = ax4.legend([object], [r'$\mathrm{This\; work}$'], handler_map={object: AnyObjectHandler()}, frameon=False, loc=2)


class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], linestyle='dotted', color='red')
        return [l1]


legend4 = ax4.legend([object], [r'$\mathrm{HWT15}$'], handler_map={object: AnyObjectHandler()}, frameon=False, loc=1)

ax4.add_artist(legend3)
ax4.add_artist(legend4)
ax4.legend(frameon=False, loc=3)

########################################################################################################################

# Read observational data from MLD16 #
MLD16 = np.loadtxt('./Obs_Data/MLD16.txt', skiprows=2)

# Plot L-Galaxies and MCMC data #
ax6.plot(xDM, np.log10(yDM), color='blue')
ax6.plot(xDMHen, np.log10(yDMHen), linestyle='dotted', color='blue')

# Plot observational data from #
ytop = np.array(np.log10(MLD16[0:34, 30]) - 3 * np.log10(obsHubble)) - np.array(np.log10(MLD16[0:34, 29]) - 3 * np.log10(obsHubble))
ybot = np.array(np.log10(MLD16[0:34, 29]) - 3 * np.log10(obsHubble)) - np.array(np.log10(MLD16[0:34, 28]) - 3 * np.log10(obsHubble))
ax6.errorbar(MLD16[0:34, 0] + 2 * np.log10(obsHubble), np.log10(MLD16[0:34, 28]) - 3 * np.log10(obsHubble), color='blue', yerr=(ybot, ytop),
             marker='o', linestyle="None", elinewidth=1, capsize=4, capthick=1, zorder=3, label=r'$\mathrm{Moffett+16\,b}$')


# Create the legends #
# Create an object to combine the Hen15 lines in the legend #
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], color='blue')
        return [l1]


legend5 = ax6.legend([object], [r'$\mathrm{This\; work}$'], handler_map={object: AnyObjectHandler()}, frameon=False, loc=2)


class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], linestyle='dotted', color='blue')
        return [l1]


legend6 = ax6.legend([object], [r'$\mathrm{HWT15}$'], handler_map={object: AnyObjectHandler()}, frameon=False, loc=1)

ax6.add_artist(legend5)
ax6.add_artist(legend6)
ax6.legend(frameon=False, loc=3)

# Create arrows to show how we split each component #
ax1.annotate("", xy=(0.5, 0.07), xycoords='axes fraction', xytext=(0.8, 0.5),
             arrowprops=dict(arrowstyle="-", color='black', connectionstyle="arc3,rad=0"))

ax3.annotate("", xy=(0.5, 0.07), xycoords='axes fraction', xytext=(0.2, 0.5),
             arrowprops=dict(arrowstyle="-", color='black', connectionstyle="arc3,rad=0"))

ax2.annotate(r'$\mathrm{(1,2)}$', xy=(0.82, 0.78), xycoords='axes fraction', size=17)
ax2.annotate(r'$\mathrm{Stars}$', xy=(0.431, 1.01), xycoords='axes fraction', color='black', size=17)

ax4.annotate(r'$\mathrm{(2,1)}$', xy=(0.82, 0.78), xycoords='axes fraction', size=17)
ax4.annotate(r'$\mathrm{Bulge}$', xy=(0.431, 1.01), xycoords='axes fraction', color='red', size=17)

ax6.annotate(r'$\mathrm{(2,3)}$', xy=(0.82, 0.78), xycoords='axes fraction', size=17)
ax6.annotate(r'$\mathrm{Disc}$', xy=(0.444, 1.01), xycoords='axes fraction', color='blue', size=17)

# Save the figure #
plt.savefig('/lustre/scratch/astro/di43/Python/output/58/SMF_Comp_58-' + date + '.png', bbox_inches='tight')
