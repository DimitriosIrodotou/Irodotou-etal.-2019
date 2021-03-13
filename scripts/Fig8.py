# Load data arrays #
BM = np.load(SavePath + 'BulgeMass.npy')
SM = np.load(SavePath + 'StellarMass.npy')
BulgeMassHWT15 = np.load(SavePath + 'HWT15/' + 'BulgeMass.npy')
StellarMassHWT15 = np.load(SavePath + 'HWT15/' + 'StellarMass.npy')

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure = plt.figure(0, figsize=(10, 7.5))

# Figure parameters #
plt.ylim(0.0, 1.0)
plt.xlim(3e8, 1e12)
plt.xscale('log')
plt.ylabel(r'$\mathrm{f(B/T>0.5)}$')
plt.xlabel(r'$\mathrm{M_{\bigstar} / M_{\odot}}$')

plt.tick_params(direction='in', which='both', top='on', right='on')

######################################################################################################################################################

# Convert masses to solar units #
offset = 1e10 / hubble

Bulge = BM * offset
StellarMass = SM * offset
StellarMassObs = StellarMass * np.exp((np.random.randn(len(SM)) * 0.08 * (1 + redshift)))

BulgeHWT15 = BulgeMassHWT15 * offset
StellarMassHWT15 = StellarMassHWT15 * offset
StellarMassObsHWT15 = StellarMassHWT15 * np.exp(np.random.randn(len(StellarMassHWT15)) * 0.08 * (1 + redshift))

# Divisions between bulge classes #
logRatio1 = 0.5

# Bins for histogram and plotting #
bins = np.logspace(8.9, 11.8 + 0.1, 20)

# Put galaxies into bins #
indBin = np.digitize(StellarMassObs, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((Bulge[indThisBin] / StellarMass[indThisBin]) > logRatio1)[0]) / float(allBin)

# Plot L-Galaxies data #
plt.plot(x, yBulge, color='blue', lw=2, label="$\mathrm{Irodotou+18}$")

# Put HWT15 galaxies into bins #
indBin = np.digitize(StellarMassObsHWT15, bins)
nBin = len(bins) - 1
x = np.empty(nBin)
yBulge = np.empty(nBin)

# Loop over bins, counting fractions in each class #
for iBin in range(nBin):
    x[iBin] = 0.5 * (bins[iBin] + bins[iBin + 1])
    indThisBin = np.where(indBin == iBin + 1)[0]
    allBin = len(indThisBin)

    # Bulges
    yBulge[iBin] = len(np.where((BulgeHWT15[indThisBin] / StellarMassHWT15[indThisBin]) > logRatio1)[0]) / float(allBin)

# Plot L-Galaxies data #
plt.plot(x, yBulge, color='blue', lw=2, linestyle='dotted', label="$\mathrm{Henriques+15}$")

# Create the legends #
plt.legend(frameon=False, loc=2)

######################################################################################################################################################

# Save the figure #
plt.savefig('Fig8-' + date + '-' + snap + '.png')