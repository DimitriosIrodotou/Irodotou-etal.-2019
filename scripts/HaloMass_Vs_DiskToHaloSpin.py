# Load data arrays #
Type = np.load(SavePath + 'Type.npy')
Mvir = np.load(SavePath + 'Mvir.npy')
HaloSpin = np.load(SavePath + 'HaloSpin.npy')
DiskSpin = np.load(SavePath + 'DiskSpin.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')

SavePath = '/Volumes/BAM-BLACK/output/output_ITH_Off/58/'
TypeHWT15 = np.load(SavePath + 'Type.npy')
MvirHWT15 = np.load(SavePath + 'Mvir.npy')
HaloSpinHWT15 = np.load(SavePath + 'HaloSpin.npy')
DiskSpinHWT15 = np.load(SavePath + 'DiskSpin.npy')
StellarMassHWT15 = np.load(SavePath + 'StellarMass.npy')
index = np.where((Type == 0) & (StellarMass * MassUnits > 1e9) & (Mvir * MassUnits > 1e11))
indexHWT15 = np.where((TypeHWT15 == 0) & (StellarMassHWT15 * MassUnits > 1e9) & (MvirHWT15 * MassUnits > 1e11))

M_vir = Mvir[index] * MassUnits
Halo_Spin = np.linalg.norm(HaloSpin[index], axis=1) * SpinUnits
Disk_Spin = np.linalg.norm(DiskSpin[index], axis=1) * SpinUnits

Mvir_HWT15 = MvirHWT15[indexHWT15] * MassUnits
HaloSpin_HWT15 = np.linalg.norm(HaloSpinHWT15[indexHWT15], axis=1) * SpinUnits
DiskSpin_HWT15 = np.linalg.norm(DiskSpinHWT15[indexHWT15], axis=1) * SpinUnits

RatioStars = np.divide(Disk_Spin, Halo_Spin)
RatioStarsHWT15 = np.divide(DiskSpin_HWT15, HaloSpin_HWT15)

dlog10 = 0.2

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure = plt.figure(figsize=(10, 7.5))
grid = plt.GridSpec(6, 6, hspace=0.0, wspace=0.0)
main_plot = figure.add_subplot(grid[:-1, 1:-1])
y_HistStars = figure.add_subplot(grid[:-1, 0])
x_hist = figure.add_subplot(grid[-1, 1:-1])

# Figure parameters #
x_hist.set_xscale('log')
x_hist.set_yscale('log')
main_plot.set_xscale('log')
main_plot.set_yscale('log')
y_HistStars.set_xscale('log')
y_HistStars.set_yscale('log')

x_hist.set_ylim(1e0, 5e6)
x_hist.set_xlim(1e11, 4e15)
main_plot.set_ylim(2e-3, 2e2)
main_plot.set_xlim(1e11, 4e15)
y_HistStars.set_xlim(1e0, 5e6)
y_HistStars.set_ylim(2e-3, 2e2)

main_plot.set_xticklabels([])
main_plot.set_yticklabels([])

x_hist.yaxis.set_ticks_position("right")
y_HistStars.xaxis.set_ticks_position("top")

x_hist.set_xlabel(r'$\mathrm{M_{vir} / M_{\odot}}$')
y_HistStars.set_ylabel(r'$\mathrm{|\vec{J}_{d,\bigstar}|/|\vec{J}_{halo}|}$')

x_hist.tick_params(direction='in', which='both', top='on', left='on', )
main_plot.tick_params(direction='in', which='both', top='on', right='on')
y_HistStars.tick_params(direction='in', which='both', bottom='on', right='on')

######################################################################################################################################################

# Calculate median and 1-sigma #
log10X = np.log10(M_vir)
log10XMax = np.log10(max(M_vir))
log10XMin = np.log10(min(M_vir))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(np.absolute(M_vir)[index])
    if len(index) > 0:
        median[i] = np.nanmedian(RatioStars[index])
        slow[i] = np.nanpercentile(RatioStars[index], 15.87)
        shigh[i] = np.nanpercentile(RatioStars[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
medianStars, = main_plot.plot(X, median, color='black', lw=lw)
main_plot.fill_between(X, shigh, slow, color='black', alpha='0.2', zorder=2)
fill, = main_plot.fill(np.NaN, np.NaN, color='black', alpha=0.5)

legend1 = main_plot.legend([medianStars, fill],
                           [r'$\mathrm{This\; work: Median}$', r'$\mathrm{This\; work: 16^{th}-84^{th}\,\%ile}$'],
                           frameon=False, loc=2)

# Calculate median and 1-sigma #
log10X = np.log10(Mvir_HWT15)
log10XMax = np.log10(max(Mvir_HWT15))
log10XMin = np.log10(min(Mvir_HWT15))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(np.absolute(Mvir_HWT15)[index])
    if len(index) > 0:
        median[i] = np.nanmedian(RatioStarsHWT15[index])
        slow[i] = np.nanpercentile(RatioStarsHWT15[index], 15.87)
        shigh[i] = np.nanpercentile(RatioStarsHWT15[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
medianStarsHWT15, = main_plot.plot(X, median, color='black', lw=lw, linestyle="dotted")
main_plot.fill_between(X, shigh, slow, color='black', hatch="\\", alpha='0.1', zorder=2)
fillHWT15, = main_plot.fill(np.NaN, np.NaN, color='black', hatch="\\", alpha=0.3)

# Create the legends #
legend2 = main_plot.legend([medianStarsHWT15, fillHWT15],
                           [r'$\mathrm{HWT15: Median}$', r'$\mathrm{HWT15: 16^{th}-84^{th}\,\%ile}$'], frameon=False,
                           loc=4)

main_plot.add_artist(legend1)
main_plot.add_artist(legend2)

# Create the legends #
plt.legend(frameon=False, loc=2)

######################################################################################################################################################

x_hist.hist(M_vir, bins=np.logspace(9, 16, 50), histtype='step', orientation='vertical', color='black')
x_hist.hist(Mvir_HWT15, bins=np.logspace(9, 16, 50), histtype='step', orientation='vertical', color='black',
            linestyle='dotted')
x_hist.invert_yaxis()

y_HistStars.hist(RatioStars, bins=np.logspace(-3, 3, 50), histtype='step', orientation='horizontal', color='black')
y_HistStars.hist(RatioStarsHWT15, bins=np.logspace(-3, 3, 50), histtype='step', orientation='horizontal', color='black',
                 linestyle='dotted')
y_HistStars.invert_xaxis()

######################################################################################################################################################

# Save the figure #
plt.savefig('Mvir_Vs_DTHS-58' + date + '.pdf', bbox_inches='tight')
