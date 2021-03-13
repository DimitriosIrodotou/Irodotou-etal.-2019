# Load data arrays #
Sfr = np.load(SavePath + 'Sfr.npy')
Type = np.load(SavePath + 'Type.npy')
DiskMass = np.load(SavePath + 'DiskMass.npy')
BulgeMass = np.load(SavePath + 'BulgeMass.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')
SfrHWT15 = np.load(SavePath + 'HWT15/' + 'Sfr.npy')
TypeHWT15 = np.load(SavePath + 'HWT15/' + 'Type.npy')
DiskMassHWT15 = np.load(SavePath + 'HWT15/' + 'DiskMass.npy')
BulgeMassHWT15 = np.load(SavePath + 'HWT15/' + 'BulgeMass.npy')
StellarMassHWT15 = np.load(SavePath + 'HWT15/' + 'StellarMass.npy')

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(10, 7.5))
figure.subplots_adjust(hspace=0, wspace=0)

# Figure parameters #
ax1.set_ylim(0.0, 1.0)
ax2.set_ylim(0.0, 1.0)
ax3.set_ylim(0.0, 1.0)
ax4.set_ylim(0.0, 1.0)

ax1.set_xscale('log')
ax2.set_xscale('log')
ax3.set_xscale('log')
ax4.set_xscale('log')

plt.ylabel(r'$\mathrm{f(B/T>0.5)}$')
plt.xlabel(r'$\mathrm{M_{\bigstar} / M_{\odot}}$')

ax1.set_ylabel(r'$\mathrm{f(B/T>0.5)}$')
ax2.set_ylabel(r'$\mathrm{f(B/T>0.5)}$')
ax3.set_ylabel(r'$\mathrm{f(B/T>0.5)}$')
ax4.set_ylabel(r'$\mathrm{f(B/T>0.5)}$')
ax1.set_xlabel(r'$\mathrm{M_{\bigstar} / M_{\odot}}$')
ax2.set_xlabel(r'$\mathrm{M_{\bigstar} / M_{\odot}}$')
ax3.set_xlabel(r'$\mathrm{M_{\bigstar} / M_{\odot}}$')
ax4.set_xlabel(r'$\mathrm{M_{\bigstar} / M_{\odot}}$')

ax1.tick_params(direction='in', which='both', top='on', right='on')
ax2.tick_params(direction='in', which='both', top='on', right='on')
ax3.tick_params(direction='in', which='both', top='on', right='on')
ax4.tick_params(direction='in', which='both', top='on', right='on')

# Change labels's position #
ax2.yaxis.tick_right()
ax4.yaxis.tick_right()

# Change ticks's position #
ax2.yaxis.set_label_position("right")
ax4.yaxis.set_label_position("right")

# Remove unwanted ticks/labels from subplots #
ax1.set_xticklabels([])
ax2.set_xticklabels([])

######################################################################################################################################################

# Convert masses to solar units #
offset = 1e10 / hubble

index = np.where((np.divide(Sfr, StellarMass) > 1e-2) & (Type == 0))
BM = BulgeMass[index] * offset
SM = StellarMass[index] * offset

index = np.where((np.divide(SfrHWT15, StellarMassHWT15) > 1e-2) & (TypeHWT15 == 0))
BMHWT15 = BulgeMassHWT15[index] * offset
SMHWT15 = StellarMassHWT15[index] * offset

# Divisions between bulge classes #
logRatio = 0.5

# Bins for histogram and plotting #
bins = np.logspace(8.9, 11.8 + 0.1, 20)

# Put galaxies into bins #
indBin = np.digitize(SM, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((BM[indThisBin] / SM[indThisBin]) > logRatio)[0]) / float(allBin)

# Plot L-Galaxies data #
ax1.plot(x, yBulge, color='blue', lw=2, label="$\mathrm{Irodotou+18}$")

# Put HWT15 galaxies into bins #
indBin = np.digitize(SMHWT15, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((BMHWT15[indThisBin] / SMHWT15[indThisBin]) > logRatio)[0]) / float(allBin)

# Plot L-Galaxies data #
ax1.plot(x, yBulge, color='blue', lw=2, linestyle='dotted', label="$\mathrm{Henriques+15}$")

# Save the figure #
ax1.legend(frameon=False, loc=2)
ax1.annotate(r'$\mathrm{SF-C}$', xy=(0.01, 0.5), xycoords='axes fraction', color='red', size=15)

######################################################################################################################################################

index = np.where((np.divide(Sfr, StellarMass) < 1e-2) & (Type == 0))
BM = BulgeMass[index] * offset
SM = StellarMass[index] * offset

index = np.where((np.divide(SfrHWT15, StellarMassHWT15) < 1e-2) & (TypeHWT15 == 0))
BMHWT15 = BulgeMassHWT15[index] * offset
SMHWT15 = StellarMassHWT15[index] * offset

# Put galaxies into bins #
indBin = np.digitize(SM, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((BM[indThisBin] / SM[indThisBin]) > logRatio)[0]) / float(allBin)

# Plot L-Galaxies data #
ax2.plot(x, yBulge, color='blue', lw=2, label="$\mathrm{Irodotou+18}$")

# Put HWT15 galaxies into bins #
indBin = np.digitize(SMHWT15, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((BMHWT15[indThisBin] / SMHWT15[indThisBin]) > logRatio)[0]) / float(allBin)

# Plot L-Galaxies data #
ax2.plot(x, yBulge, color='blue', lw=2, linestyle='dotted', label="$\mathrm{Henriques+15}$")
ax2.annotate(r'$\mathrm{Q-C}$', xy=(0.01, 0.5), xycoords='axes fraction', color='red', size=15)

######################################################################################################################################################

index = np.where((np.divide(Sfr, StellarMass) > 1e-2) & (Type != 0))
BM = BulgeMass[index] * offset
SM = StellarMass[index] * offset

index = np.where((np.divide(SfrHWT15, StellarMassHWT15) > 1e-2) & (TypeHWT15 != 0))
BMHWT15 = BulgeMassHWT15[index] * offset
SMHWT15 = StellarMassHWT15[index] * offset

# Put galaxies into bins #
indBin = np.digitize(SM, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((BM[indThisBin] / SM[indThisBin]) > logRatio)[0]) / float(allBin)

# Plot L-Galaxies data #
ax3.plot(x, yBulge, color='blue', lw=2, label="$\mathrm{Irodotou+18}$")

# Put HWT15 galaxies into bins #
indBin = np.digitize(SMHWT15, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((BMHWT15[indThisBin] / SMHWT15[indThisBin]) > logRatio)[0]) / float(allBin)

# Plot L-Galaxies data #
ax3.plot(x, yBulge, color='blue', lw=2, linestyle='dotted', label="$\mathrm{Henriques+15}$")
ax3.annotate(r'$\mathrm{SF-S}$', xy=(0.01, 0.5), xycoords='axes fraction', color='red', size=15)

######################################################################################################################################################

index = np.where((np.divide(Sfr, StellarMass) < 1e-2) & (Type != 0))
BM = BulgeMass[index] * offset
SM = StellarMass[index] * offset

index = np.where((np.divide(SfrHWT15, StellarMassHWT15) < 1e-2) & (TypeHWT15 != 0))
BMHWT15 = BulgeMassHWT15[index] * offset
SMHWT15 = StellarMassHWT15[index] * offset

# Put galaxies into bins #
indBin = np.digitize(SM, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((BM[indThisBin] / SM[indThisBin]) > logRatio)[0]) / float(allBin)

# Plot L-Galaxies data #
ax4.plot(x, yBulge, color='blue', lw=2, label="$\mathrm{Irodotou+18}$")

# Put HWT15 galaxies into bins #
indBin = np.digitize(SMHWT15, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((BMHWT15[indThisBin] / SMHWT15[indThisBin]) > logRatio)[0]) / float(allBin)

# Plot L-Galaxies data #
ax4.plot(x, yBulge, color='blue', lw=2, linestyle='dotted', label="$\mathrm{Henriques+15}$")

ax4.annotate(r'$\mathrm{Q-S}$', xy=(0.01, 0.5), xycoords='axes fraction', color='red', size=15)

######################################################################################################################################################

# Save the figure #
plt.savefig('Fig9-' + date + '-' + snap + '.png')