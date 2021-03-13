# Load data arrays #
DiskMass = np.load(SavePath + 'DiskMass.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')
SHMR = np.load(SavePath + 'StellarHalfMassRadius.npy')

DiskMassHWT15 = np.load(SavePath + 'NoDI/' + 'DiskMass.npy')
StellarMassHWT15 = np.load(SavePath + 'NoDI/' + 'StellarMass.npy')
SHMRHWT15 = np.load(SavePath + 'NoDI/' + 'StellarHalfMassRadius.npy')

DiskMassNoGS = np.load(SavePath + 'NoDI/' + 'DiskMass.npy')
StellarMassNoGS = np.load(SavePath + 'NoDI/' + 'StellarMass.npy')
SHMRNoGS = np.load(SavePath + 'NoDI/' + 'StellarHalfMassRadius.npy')

DiskMassNoDI = np.load(SavePath + 'NoDI/' + 'DiskMass.npy')
StellarMassNoDI = np.load(SavePath + 'NoDI/' + 'StellarMass.npy')
SHMRNoDI = np.load(SavePath + 'NoDI/' + 'StellarHalfMassRadius.npy')

# Trim the data #
index = np.where(DiskMass > 0.7 * StellarMass)
indexHWT15 = np.where(DiskMassHWT15 > 0.7 * StellarMassHWT15)
indexNoDI = np.where(DiskMassNoDI > 0.7 * StellarMassNoDI)
indexNoGS = np.where(DiskMassNoGS > 0.7 * StellarMassNoGS)

S_H_M_R = SHMR[index] * LengthUnits
Stellar_Mass = StellarMass[index] * MassUnits

Stellar_Mass_HWT15 = DiskMassHWT15[indexHWT15] * MassUnits
S_H_M_R_HWT15 = SHMRHWT15[indexHWT15] * LengthUnits

Stellar_Mass_NoDI = DiskMassNoDI[indexNoDI] * MassUnits
S_H_M_R_NoDI = SHMRNoDI[indexNoDI] * LengthUnits

Stellar_Mass_NoGS = DiskMassNoGS[indexNoGS] * MassUnits
S_H_M_R_NoGS = SHMRNoGS[indexNoGS] * LengthUnits

dlog10 = 0.3

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure = plt.figure(0, figsize=(10, 7.5))

# Figure parameters #
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e9, 1e12)
plt.ylim(1e-1, 1e2)
plt.ylabel(r'$\mathrm{R_{HM} / kpc}$')
plt.xlabel(r'$\mathrm{M_{\bigstar} / M_{\odot}}$')

plt.tick_params(direction='in', which='both', top='on', right='on')

######################################################################################################################################################

# Calculate median and 1-sigma #
log10X = np.log10(Stellar_Mass_HWT15)
log10XMax = np.log10(max(Stellar_Mass_HWT15))
log10XMin = np.log10(min(Stellar_Mass_HWT15))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Stellar_Mass_HWT15[index])
    if len(index) > 0:
        median[i] = np.median(S_H_M_R_HWT15[index])
        slow[i] = np.percentile(S_H_M_R_HWT15[index], 15.87)
        shigh[i] = np.percentile(S_H_M_R_HWT15[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = plt.plot(X, median, color='red', lw=lw, linestyle='dotted', label=r'$\mathrm{HWT15: Median}$')

# Calculate median and 1-sigma #
log10X = np.log10(Stellar_Mass_NoDI)
log10XMax = np.log10(max(Stellar_Mass_NoDI))
log10XMin = np.log10(min(Stellar_Mass_NoDI))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Stellar_Mass_NoDI[index])
    if len(index) > 0:
        median[i] = np.median(S_H_M_R_NoDI[index])
        slow[i] = np.percentile(S_H_M_R_NoDI[index], 15.87)
        shigh[i] = np.percentile(S_H_M_R_NoDI[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = plt.plot(X, median, color='blue', lw=lw, linestyle='dashdot', label=r'$\mathrm{This\; work: old\; DI}$')

# Calculate median and 1-sigma #
log10X = np.log10(Stellar_Mass_NoGS)
log10XMax = np.log10(max(Stellar_Mass_NoGS))
log10XMin = np.log10(min(Stellar_Mass_NoGS))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Stellar_Mass_NoGS[index])
    if len(index) > 0:
        median[i] = np.median(S_H_M_R_NoGS[index])
        slow[i] = np.percentile(S_H_M_R_NoGS[index], 15.87)
        shigh[i] = np.percentile(S_H_M_R_NoGS[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = plt.plot(X, median, color='green', lw=lw, linestyle='dashdot', label=r'$\mathrm{This\; work: f=1.0}$')

# Calculate median and 1-sigma #
log10X = np.log10(Stellar_Mass)
log10XMax = np.log10(max(Stellar_Mass))
log10XMin = np.log10(min(Stellar_Mass))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(np.absolute(Stellar_Mass)[index])
    if len(index) > 0:
        median[i] = np.nanmedian(S_H_M_R[index])
        slow[i] = np.nanpercentile(S_H_M_R[index], 15.87)
        shigh[i] = np.nanpercentile(S_H_M_R[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = plt.plot(X, median, color='red', lw=lw, label=r'$\mathrm{This\; work: Median}$')

# Create the legends #
plt.legend(frameon=False, ncol=1, loc=2)

######################################################################################################################################################

# Save the figure #
plt.savefig('DM_Vs_SHMR2_58-' + date + '.png', bbox_inches='tight')