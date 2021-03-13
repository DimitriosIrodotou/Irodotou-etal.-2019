# Load data arrays #
DiskMass = np.load(SavePath + 'DiskMass.npy')
DiskSpin = np.load(SavePath + 'DiskSpin.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')

DiskMassHWT15 = np.load(SavePath + 'HWT15/' + 'DiskMass.npy')
DiskSpinHWT15 = np.load(SavePath + 'HWT15/' + 'DiskSpin.npy')
StellarMassHWT15 = np.load(SavePath + 'HWT15/' + 'StellarMass.npy')

DiskMassNoGS = np.load(SavePath + 'NoGS/' + 'DiskMass.npy')
DiskSpinNoGS = np.load(SavePath + 'NoGS/' + 'DiskSpin.npy')
StellarMassNoGS = np.load(SavePath + 'NoGS/' + 'StellarMass.npy')

DiskMassNoDI = np.load(SavePath + 'NoDI/' + 'DiskMass.npy')
DiskSpinNoDI = np.load(SavePath + 'NoDI/' + 'DiskSpin.npy')
StellarMassNoDI = np.load(SavePath + 'NoDI/' + 'StellarMass.npy')

# Trim the data #
index = np.where(DiskMass > 0.7 * StellarMass)
indexHWT15 = np.where(DiskMassHWT15 > 0.7 * StellarMassHWT15)
indexNoDI = np.where(DiskMassNoDI > 0.7 * StellarMassNoDI)
indexNoGS = np.where(DiskMassNoGS > 0.7 * StellarMassNoGS)

Disk_Mass = DiskMass[index] * MassUnits
Disk_Spin = np.linalg.norm(DiskSpin[index], axis=1) * SpinUnits

Disk_Mass_HWT15 = DiskMassHWT15[indexHWT15] * MassUnits
Disk_Spin_HWT15 = np.linalg.norm(DiskSpinHWT15[indexHWT15], axis=1) * SpinUnits

Disk_Mass_NoDI = DiskMassNoDI[indexNoDI] * MassUnits
Disk_Spin_NoDI = np.linalg.norm(DiskSpinNoDI[indexNoDI], axis=1) * SpinUnits

Disk_Mass_NoGS = DiskMassNoGS[indexNoGS] * MassUnits
Disk_Spin_NoGS = np.linalg.norm(DiskSpinNoGS[indexNoGS], axis=1) * SpinUnits

dlog10 = 0.1

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure = plt.figure(0, figsize=(10, 7.5))

# Figure parameters #
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e1, 1e5)
plt.xlim(1e9, 1e12)
plt.xlabel(r'$\mathrm{M_{d,\bigstar} / M_{\odot}}$')
plt.ylabel(r'$\mathrm{|j_{d,\bigstar}| / (km \cdot s^{-1} \cdot kpc)}$')
plt.tick_params(direction='in', which='both', top='on', right='on')

######################################################################################################################################################

# Calculate median and 1-sigma #
log10X = np.log10(Disk_Mass_HWT15)
log10XMax = np.log10(max(Disk_Mass_HWT15))
log10XMin = np.log10(min(Disk_Mass_HWT15))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Disk_Mass_HWT15[index])
    if len(index) > 0:
        median[i] = np.median(Disk_Spin_HWT15[index])
        slow[i] = np.percentile(Disk_Spin_HWT15[index], 15.87)
        shigh[i] = np.percentile(Disk_Spin_HWT15[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = plt.plot(X, median, color='black', lw=lw, linestyle='dotted', label=r'$\mathrm{HWT15: Median}$')

# Calculate median and 1-sigma #
log10X = np.log10(Disk_Mass_NoDI)
log10XMax = np.log10(max(Disk_Mass_NoDI))
log10XMin = np.log10(min(Disk_Mass_NoDI))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Disk_Mass_NoDI[index])
    if len(index) > 0:
        median[i] = np.median(Disk_Spin_NoDI[index])
        slow[i] = np.percentile(Disk_Spin_NoDI[index], 15.87)
        shigh[i] = np.percentile(Disk_Spin_NoDI[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = plt.plot(X, median, color='blue', lw=lw, linestyle='dashdot', label=r'$\mathrm{This\; work: old\; DI}$')

# Calculate median and 1-sigma #
log10X = np.log10(Disk_Mass_NoGS)
log10XMax = np.log10(max(Disk_Mass_NoGS))
log10XMin = np.log10(min(Disk_Mass_NoGS))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(Disk_Mass_NoGS[index])
    if len(index) > 0:
        median[i] = np.median(Disk_Spin_NoGS[index])
        slow[i] = np.percentile(Disk_Spin_NoGS[index], 15.87)
        shigh[i] = np.percentile(Disk_Spin_NoGS[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = plt.plot(X, median, color='red', lw=lw, linestyle='dashdot', label=r'$\mathrm{This\; work: f=1.0}$')

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
median, = plt.plot(X, median, color='black', lw=lw, label=r'$\mathrm{This\; work: Median}$')

# Create the legends #
plt.legend(frameon=False, ncol=1, loc=2)

######################################################################################################################################################

# Save the figure #
plt.savefig('DM_Vs_DS2_58-' + date + '.png', bbox_inches='tight')