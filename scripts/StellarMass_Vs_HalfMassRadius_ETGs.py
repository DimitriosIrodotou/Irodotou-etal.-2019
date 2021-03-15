# Load data arrays #
CBMassMajor = np.load(SavePath + 'CBMassMajor.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')
SHMR = np.load(SavePath + 'StellarHalfMassRadius.npy')

SavePath = '/Volumes/BAM-BLACK/output/output_ITH_Off/58/'
CBMassMajorHWT15 = np.load(SavePath + 'CBMassMajor.npy')
StellarMassHWT15 = np.load(SavePath + 'StellarMass.npy')
SHMRHWT15 = np.load(SavePath + 'StellarHalfMassRadius.npy')

# Trim the data #
index = np.where(CBMassMajor > 0.7 * StellarMass)
indexHWT15 = np.where(CBMassMajorHWT15 > 0.7 * StellarMassHWT15)

S_H_M_R = SHMR[index] * LengthUnits
Stellar_Mass = StellarMass[index] * MassUnits
S_H_M_RHWT15 = SHMRHWT15[indexHWT15] * LengthUnits
Stellar_MassHWT15 = StellarMassHWT15[indexHWT15] * MassUnits

dlog10 = 0.3

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure, ax = plt.subplots()
figure, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20, 7.5))
figure.subplots_adjust(hspace=0, wspace=0.37)

# Figure parameters #
ax1.set_xscale('log')
ax2.set_xscale('log')
ax1.set_yscale('log')
ax2.set_yscale('log')

ax1.set_xlim(1e9, 1e12)
ax2.set_xlim(1e9, 1e12)
ax1.set_ylim(1e-1, 1e2)
ax2.set_ylim(1e-1, 1e2)

ax1.set_ylabel(r'$\mathrm{R_{HM} / kpc}$', size=25)
ax2.set_ylabel(r'$\mathrm{R_{HM} / kpc}$', size=25)
ax1.set_xlabel(r'$\mathrm{M_{\bigstar} / M_{\odot}}$', size=25)
ax2.set_xlabel(r'$\mathrm{M_{\bigstar} / M_{\odot}}$', size=25)

ax1.set_xticklabels([])
ax1.tick_params(direction='in', which='both', top='on', right='on', labelsize=25)
ax2.tick_params(direction='in', which='both', top='on', right='on', labelsize=25)

######################################################################################################################################################

# Read observational data from CCW10, FSS17, CMA13, KCR17SP, ZY17, CALIFAETGs, SMW03ETG and G09 #
CCW10 = np.genfromtxt('./Obs_Data/CCW10.csv', delimiter=',', names=['Mstar', 'R50'])
FSS17 = np.genfromtxt('./Obs_Data/FSS17.csv', delimiter=',', names=['Mstar', 'R50'])
CMA13 = np.genfromtxt('./Obs_Data/CMA13.csv', delimiter=',', names=['Mstar', 'R50'])
# KCR17SP = np.genfromtxt('./Obs_Data/KCR17SP.csv', delimiter=',', names=['Mstar', 'R50'])
ZY17 = np.genfromtxt('./Obs_Data/ZY17Mass_Red.csv', delimiter=',', names=['Mstar', 'R50'])
# CALIFAETGs = np.genfromtxt('./Obs_Data/CALIFAETGs.csv', delimiter=',', names=['Mstar', 'R50'])
SMW03ETGst = np.genfromtxt('./Obs_Data/SMW03_ETGst.csv', delimiter=',', names=['Mstar', 'R50'])
SMW03ETGsm = np.genfromtxt('./Obs_Data/SMW03_ETGsm.csv', delimiter=',', names=['Mstar', 'R50'])
SMW03ETGsb = np.genfromtxt('./Obs_Data/SMW03_ETGsb.csv', delimiter=',', names=['Mstar', 'R50'])
ytop = np.array(SMW03ETGst['R50']) - np.array(SMW03ETGsm['R50'])
ybot = np.array(SMW03ETGsm['R50']) - np.array(SMW03ETGsb['R50'])
G09 = np.genfromtxt('./Obs_Data/G09.txt', delimiter=';', skip_header=1,
                    dtype=[('Type', 'U15'), ('Type2', 'U15'), ('h', 'f8'), ('re', 'f8'), ('BtoT', 'f8'),
                           ('MtoLb', 'f8'), ('MtoLd', 'f8'), ('Mb', 'f8'), ('Md', 'f8')])

# Plot observational data from SMW03ETG,CCW10, FSS17, ZY17, CALIFAETGs, CMA13 and KCR17SP#
scatter1 = ax1.errorbar(np.power(10, SMW03ETGsm['Mstar']), SMW03ETGsm['R50'], yerr=(ybot, ytop), color='lime',
                        marker='o', linestyle="None", elinewidth=1, capsize=4, capthick=1, zorder=3)
scatter2 = ax1.scatter(np.power(10, CCW10['Mstar']), np.power(10, CCW10['R50']), edgecolor='black', color='red',
                       s=2 * 50, marker='*', zorder=2)
scatter3 = ax1.scatter(np.power(10, ZY17['Mstar']) / hubble ** 2, ZY17['R50'], edgecolor='black', color='cyan', s=50,
                       marker='d', zorder=3)
scatter4 = ax1.scatter(np.power(10, FSS17['Mstar']), np.power(10, FSS17['R50']), edgecolor='black', color='blue',
                       s=50, marker='p', zorder=3)
# scatter5 = ax1.scatter(np.power(10, CALIFAETGs['Mstar']), CALIFAETGs['R50'], c='g', marker='s', s=50/2, zorder=2)
scatter6 = ax1.scatter(np.power(10, CMA13['Mstar']), CMA13['R50'], edgecolor='black', color='darkorange', s=50,
                       marker='+', zorder=3)
# scatter7 = ax1.scatter(np.power(10, KCR17SP['Mstar']), np.power(10, KCR17SP['R50']), color='magenta', s=50,
#                        marker='^', zorder=3)

index = np.where(G09['Type'] == 'elliptical  ')
G09_Stellar_Mass = G09['Mb'][index]
scatter8 = ax1.scatter(G09_Stellar_Mass, G09['re'][index], edgecolor='black', color='Purple', s=50, marker='s',
                       zorder=2)

# Plot L-Galaxies data - 2D histogram #
hexbin = ax1.hexbin(Stellar_Mass, S_H_M_R, xscale='log', yscale='log', bins='log', cmap='Greys', mincnt=2)

# Adjust the color bar #
cbaxes = figure.add_axes([0.452, 0.11, 0.01, 0.77])
cb = plt.colorbar(hexbin, cax=cbaxes)
cbaxes.tick_params(direction='out', which='both', right='on', labelsize=25)
cb.set_label('$\mathrm{Counts\; per\; hexbin}$', size=25)

# Create the legends #
legend1 = ax1.legend([scatter1, scatter2, scatter3, scatter4, scatter6, scatter8],
                     [r'$\mathrm{Shen+03}$', r'$\mathrm{Chen+10}$', r'$\mathrm{Zhang+19}$', r'$\mathrm{Forbes+17}$',
                      r'$\mathrm{Cappellari+13b}$', r'$\mathrm{Gadotti\, 09: Ellipticals}$'], ncol=2, scatterpoints=3,
                     frameon=False, loc=2)

colors = ['black', 'grey', 'lightgrey']
squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
legend2 = ax1.legend([squares], [r'$\mathrm{This\;work:M_{cb(ma)} / M_{\bigstar} > 0.7}}$'], scatterpoints=len(colors),
                     scatteryoffsets=[.5], handlelength=len(colors), markerscale=2, frameon=False, loc=4)

ax1.add_artist(legend1)
ax1.add_artist(legend2)

######################################################################################################################################################

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
median, = ax2.plot(X, median, color='black', lw=lw, zorder=5)
ax2.fill_between(X, shigh, slow, color='black', alpha='0.5')
fill, = plt.fill(np.NaN, np.NaN, color='black', alpha=0.5)

# Create the legends #
legend4 = ax2.legend([median, fill],
                     [r'$\mathrm{This\; work: Median}$', r'$\mathrm{This\; work: 16^{th}-84^{th}\,\%ile}$'],
                     frameon=False, loc=2)

# Calculate median and 1-sigma #
log10X = np.log10(Stellar_MassHWT15)
log10XMax = np.log10(max(Stellar_MassHWT15))
log10XMin = np.log10(min(Stellar_MassHWT15))
nbin = int((log10XMax - log10XMin) / dlog10)
X = np.empty(nbin)
median = np.empty(nbin)
slow = np.empty(nbin)
shigh = np.empty(nbin)
log10XLow = log10XMin
for i in range(nbin):
    index = np.where((log10X >= log10XLow) & (log10X < log10XLow + dlog10))[0]
    X[i] = np.mean(np.absolute(Stellar_MassHWT15)[index])
    if len(index) > 0:
        median[i] = np.nanmedian(S_H_M_RHWT15[index])
        slow[i] = np.nanpercentile(S_H_M_RHWT15[index], 15.87)
        shigh[i] = np.nanpercentile(S_H_M_RHWT15[index], 84.13)
    log10XLow += dlog10

# Plot median and 1-sigma lines #
medianHWT15, = ax2.plot(X, median, color='black', lw=lw, linestyle='dotted')

# Read observational data from BDL11E and LDR15E #
# BDL11Em = np.genfromtxt('./Obs_Data/BDL11Em.csv', delimiter=',', names=['Mstar', 'R50'])
# BDL11Eb = np.genfromtxt('./Obs_Data/BDL11Eb.csv', delimiter=',', names=['Mstar', 'R50'])
# BDL11Et = np.genfromtxt('./Obs_Data/BDL11Et.csv', delimiter=',', names=['Mstar', 'R50'])
LDR15Eb = np.genfromtxt('./Obs_Data/LDR15Eb.csv', delimiter=',', names=['Mstar', 'R50'])
LDR15Ea = np.genfromtxt('./Obs_Data/LDR15Ea.csv', delimiter=',', names=['Mstar', 'R50'])
G09b = np.genfromtxt('./Obs_Data/G09b.csv', delimiter=',', names=['Mstar', 'R50'])

# Plot observational data from BDL11E and LDR15E #
line1, = ax2.plot(10 ** G09b['Mstar'], G09b['R50'], color='blue', lw=3)
line2, = ax2.plot(LDR15Eb['Mstar'], LDR15Eb['R50'], color='lime', lw=3)
line3, = ax2.plot(LDR15Ea['Mstar'], LDR15Ea['R50'], color='red', lw=3)
# ax2.plot(np.power(10, BDL11Em['Mstar']), np.power(10, (BDL11Em['R50'] - 3)), color='black', lw=lw)
# ax2.plot(np.power(10, BDL11Eb['Mstar']), np.power(10, (BDL11Eb['R50'] - 3)), color='black', lw=lw, linestyle='dashed')
# ax2.plot(np.power(10, BDL11Et['Mstar']), np.power(10, (BDL11Et['R50'] - 3)), color='black', lw=lw, linestyle='dashed')

# Create the legends #
legend5 = ax2.legend([medianHWT15], [r'$\mathrm{HWT15: Median}$'], frameon=False, loc=1)


# Create an object to combine median and 1-sigma lines from G09 in the legend #
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.9 * height, 0.9 * height], color='black', linestyle='dotted')
        l2 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], color='black')
        l3 = plt.Line2D([x0, y0 + width], [0.1 * height, 0.1 * height], color='black', linestyle='dotted')
        return [l1, l2, l3]


legend6 = ax2.legend([line1, line2, line3],
                     [r'$\mathrm{Gadotti\, 09: Ellipticals}$', r'$\mathrm{Lange+15:Single\; p.l.}$',
                      r'$\mathrm{Lange+15:Double\; p.l.}$'], handler_map={object:AnyObjectHandler()}, frameon=False,
                     loc=4)

ax2.add_artist(legend6)
ax2.add_artist(legend4)
ax2.add_artist(legend5)

######################################################################################################################################################

# Save the figure #
plt.savefig('SM_Vs_SHMR_ETGs_55-' + date + '.pdf', bbox_inches='tight')
