# Load data arrays #
Mvir = np.load(SavePath + 'Mvir.npy')
SM = np.load(SavePath + 'StellarMass.npy')

MvirH = np.load(SavePath + 'HWT15/' + 'Mvir.npy')
SMH = np.load(SavePath + 'HWT15/' + 'StellarMass.npy')

SM = SM * MassUnits
Mvir = Mvir * MassUnits

SMH = SMH * MassUnits
MvirH = MvirH * MassUnits

dlog10 = 0.1

######################################################################################################################################################

# Generate initial figure #
plt.close()
fig = plt.figure(0, figsize=(10, 7.5))

# Figure parameters #
plt.xscale('log')
plt.xlim(1e11, 1e16)
plt.ylabel(r'$\mathrm{M_{\bigstar} / M_{vir}}$')
plt.xlabel(r'$\mathrm{M_{vir} / M_{\odot}}$')
plt.tick_params(direction='in', which='both', top='on', right='on')

Ratio = SM / Mvir
RatioH = SMH / MvirH

######################################################################################################################################################

# Calculate median and 1-sigma #
log10X = np.log10(Mvir)
log10XMax = np.log10(max(Mvir))
log10XMin = np.log10(min(Mvir))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(np.absolute(Mvir)[index])
    if len(index) > 0:
        median[i] = np.nanmedian(Ratio[index])
        slow[i] = np.nanpercentile(Ratio[index], 15.87)
        shigh[i] = np.nanpercentile(Ratio[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = plt.plot(X, median, color='red', lw=lw, label=r'$\mathrm{This\;work}$')
plt.fill_between(X, shigh, slow, color='red', alpha='0.5')
fill, = plt.fill(np.NaN, np.NaN, color='red', alpha=0.5)

# Calculate median and 1-sigma #
log10X = np.log10(MvirH)
log10XMax = np.log10(max(MvirH))
log10XMin = np.log10(min(MvirH))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(np.absolute(MvirH)[index])
    if len(index) > 0:
        median[i] = np.nanmedian(RatioH[index])
        slow[i] = np.nanpercentile(RatioH[index], 15.87)
        shigh[i] = np.nanpercentile(RatioH[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
median, = plt.plot(X, median, color='green', lw=lw, label=r'$\mathrm{HWT15}$')
plt.fill_between(X, shigh, slow, color='green', alpha='0.5')
fill, = plt.fill(np.NaN, np.NaN, color='green', alpha=0.5)

# Create the legends #
plt.legend(frameon=False, loc=2)

######################################################################################################################################################

# Save the figure #
plt.savefig('SM_Vs_MVir-' + date + '-' + snap + '.png')