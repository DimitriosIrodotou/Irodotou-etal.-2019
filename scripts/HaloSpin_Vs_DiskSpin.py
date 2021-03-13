# Load data arrays #
Type = np.load(SavePath + 'Type.npy')
HaloSpin = np.load(SavePath + 'HaloSpin.npy')
DiskSpin = np.load(SavePath + 'DiskSpin.npy')
ColdGasSpin = np.load(SavePath + 'ColdGasSpin.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')

index = np.where((Type == 0) & (StellarMass * MassUnits > 1e9))

Halo_Spin = np.linalg.norm(HaloSpin[index], axis=1) * SpinUnits
Disk_Spin = np.linalg.norm(DiskSpin[index], axis=1) * SpinUnits
Cold_Gas_Spin = np.linalg.norm(ColdGasSpin[index], axis=1) * SpinUnits

dlog10 = 0.2

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20, 7.5))
figure.subplots_adjust(hspace=0, wspace=0.3)

# Figure parameters #
ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_xscale('log')
ax2.set_yscale('log')

ax1.set_xlim(2.2e2, 1e6)
ax1.set_ylim(2.2e2, 1e6)
ax2.set_xlim(2.2e2, 1e6)
ax2.set_ylim(2.2e2, 1e6)

ax1.set_xlabel(r'$\mathrm{|\vec{J}_{halo}| / (km \cdot s^{-1} \cdot kpc)}$')
ax2.set_xlabel(r'$\mathrm{|\vec{J}_{halo}| / (km \cdot s^{-1} \cdot kpc)}$')
ax1.set_ylabel(r'$\mathrm{|\vec{J}_{d,gas}| / (km \cdot s^{-1} \cdot kpc)}$')
ax2.set_ylabel(r'$\mathrm{|\vec{J}_{d,\bigstar}| / (km \cdot s^{-1} \cdot kpc)}$')

ax1.tick_params(direction='in', which='both', top='on', right='on')
ax2.tick_params(direction='in', which='both', top='on', right='on')

######################################################################################################################################################

# Calculate median and 1-sigma #
log10X = np.log10(Halo_Spin)
log10XMax = np.log10(max(Halo_Spin))
log10XMin = np.log10(min(Halo_Spin))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(np.absolute(Halo_Spin)[index])
    if len(index) > 0:
        median[i] = np.nanmedian(Cold_Gas_Spin[index])
        slow[i] = np.nanpercentile(Cold_Gas_Spin[index], 15.87)
        shigh[i] = np.nanpercentile(Cold_Gas_Spin[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
medianGas, = ax1.plot(X, median, color='red', lw=lw, label=r'$\mathrm{Median}$')

# ax1.scatter(Halo_Spin, Cold_Gas_Spin, s=1, c='black')

# Calculate median and 1-sigma #
log10X = np.log10(Halo_Spin)
log10XMax = np.log10(max(Halo_Spin))
log10XMin = np.log10(min(Halo_Spin))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(np.absolute(Halo_Spin)[index])
    if len(index) > 0:
        median[i] = np.nanmedian(Disk_Spin[index])
        slow[i] = np.nanpercentile(Disk_Spin[index], 15.87)
        shigh[i] = np.nanpercentile(Disk_Spin[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
medianGas, = ax2.plot(X, median, color='red', lw=lw, label=r'$\mathrm{Median}$')

# ax2.scatter(Halo_Spin, Cold_Gas_Spin, s=1, c='black')

# Plot a 1:1 line
x = np.logspace(1, 6, 10)
ax1.plot(x, x, color='blue', label=r'$\mathrm{1:1\, line}$')
ax2.plot(x, x, color='blue', label=r'$\mathrm{1:1\, line}$')

ax1.legend(frameon=False, loc=2)
ax2.legend(frameon=False, loc=2)

# Save the figure #
plt.savefig('HS_Vs_DS-58' + date + '.png', bbox_inches='tight')