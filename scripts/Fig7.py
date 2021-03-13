# Load data arrays #
BM = np.load(SavePath + 'BulgeMass.npy')
SM = np.load(SavePath + 'StellarMass.npy')

######################################################################################################################################################

# Generate initial figure #
plt.close()
fig = plt.figure(0, figsize=(10, 7.5))

# Figure parameters #
plt.ylim(0.0, 0.6)
plt.ylabel(r'$\mathrm{f(B/T)_{\bigstar}}$')
plt.xlabel(r'$\mathrm{(B/T)_{\bigstar}}$')

plt.tick_params(direction='in', which='both', top='on', right='on')

######################################################################################################################################################

# Convert masses to solar units #
MassUnits = 1e10 / hubble
BM = BM * MassUnits
SM = SM * MassUnits

# Find galaxies with masses > 10^10 Msun #
index = np.where(SM > 1e10)

# Calculate the B/T ratio #
Ratio = np.divide(BM[index], SM[index])

# Plots BBT19 bar's midpoints #
BBT19 = np.genfromtxt('./Obs_Data/BBT19.csv', delimiter=',', names=['BT', 'f'])
# plt.scatter(BBT19['BT'], BBT19['f'], color='red', s=size, marker='_', zorder=2, label="$\mathrm{Bluck+19}$")

# Weight each bin by the total number of values and make a histogram #
Weights = np.divide(np.ones_like(Ratio), float(len(Ratio)))
plt.hist(Ratio, weights=Weights, edgecolor='black', bins=20)

# Create the legends #
plt.legend(frameon=False, loc=2)

######################################################################################################################################################

# Save the figure #
plt.savefig('Fig7-' + date + '-' + snap + '.png')