# Load data arrays #
DiskMass = np.load(SavePath + 'DiskMass.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')
SHMR = np.load(SavePath + 'StellarHalfMassRadius.npy')

SavePath = '/Volumes/BAM-BLACK/output/output_ITH_Off/58/'
DiskMassHWT15 = np.load(SavePath + 'DiskMass.npy')
StellarMassHWT15 = np.load(SavePath + 'StellarMass.npy')
SHMRHWT15 = np.load(SavePath + 'StellarHalfMassRadius.npy')

# Trim the data #
index = np.where(DiskMass > 0.7 * StellarMass)
indexHWT15 = np.where(DiskMassHWT15 > 0.7 * StellarMassHWT15)

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

# Read observational data from SMW03, ZY17,SRD18 and KCR17SR #
SMW03LTGst = np.genfromtxt('./Obs_Data/SMW03_LTGst.csv', delimiter=',', names=['Mstar', 'R50'])
SMW03LTGsm = np.genfromtxt('./Obs_Data/SMW03_LTGsm.csv', delimiter=',', names=['Mstar', 'R50'])
SMW03LTGsb = np.genfromtxt('./Obs_Data/SMW03_LTGsb.csv', delimiter=',', names=['Mstar', 'R50'])
ytop = np.array(SMW03LTGst['R50']) - np.array(SMW03LTGsm['R50'])
ybot = np.array(SMW03LTGsm['R50']) - np.array(SMW03LTGsb['R50'])
ZY17 = np.genfromtxt('./Obs_Data/ZY17Mass_Blue.csv', delimiter=',', names=['Mstar', 'R50'])
# SRD18 = np.genfromtxt('./Obs_Data/SRD18.csv', delimiter=',', names=['Mstar', 'R50'])
KCR17SR = np.genfromtxt('./Obs_Data/KCR17SR.csv', delimiter=',', names=['Mstar', 'R50'])

# Plot observational data from SMW03, ZY17,SRD18 and KCR17SR #
scatter1 = ax1.errorbar(np.power(10, SMW03LTGsm['Mstar']), SMW03LTGsm['R50'], yerr=(ybot, ytop), color='blue',
                        marker='o', linestyle="None", elinewidth=2, capsize=8, capthick=2, zorder=3)
scatter2 = ax1.scatter(np.power(10, ZY17['Mstar']) / hubble ** 2, ZY17['R50'], edgecolor='black', color='lime', s=50,
                       marker='p', zorder=2)
# scatter3 = ax1.scatter(np.power(10, SRD18['Mstar']), SRD18['R50'], color='blue', s=size, marker='*', zorder=2)
scatter4 = ax1.scatter(np.power(10, KCR17SR['Mstar']), np.power(10, KCR17SR['R50']), edgecolor='black', color='red',
                       s=50, marker='^', zorder=2)

# Plot L-Galaxies data - 2D histogram #
hexbin = ax1.hexbin(Stellar_Mass, S_H_M_R, xscale='log', yscale='log', bins='log', cmap='Greys', mincnt=2)

# Adjust the color bar #
cbaxes = figure.add_axes([0.452, 0.11, 0.01, 0.77])
cb = plt.colorbar(hexbin, cax=cbaxes)
cbaxes.tick_params(direction='out', which='both', right='on', labelsize=25)
cb.set_label('$\mathrm{Counts\; per\; hexbin}$', size=25)

# Create the legends #
legend1 = ax1.legend([scatter1, scatter2, scatter4],
                     [r'$\mathrm{Shen+03}$', r'$\mathrm{Zhang+19}$', r'$\mathrm{Kalinova+17(SR)}$'], frameon=False,
                     loc=4, scatterpoints=3)

colors = ['black', 'grey', 'lightgrey']
squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
legend2 = ax1.legend([squares], [r'$\mathrm{This\;work:M_{d,\bigstar} / M_{\bigstar}> 0.7}$'],
                     scatterpoints=len(colors), scatteryoffsets=[.5], handlelength=len(colors), markerscale=2,
                     frameon=False, loc=2)

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
ax2.fill_between(X, shigh, slow, color='black', alpha=0.5)
fill, = plt.fill(np.NaN, np.NaN, color='black', alpha=0.5)

# Create the legends #
legend3 = ax2.legend([median, fill],
                     [r'$\mathrm{This\; work: Median}$', r'$\mathrm{This\; work:16^{th}-84^{th}\,\%ile}$'],
                     frameon=False, loc=2)

# Read observational data from BDL11L and LDR15L #
BDL11Lm = np.genfromtxt('./Obs_Data/BDL11Lm.csv', delimiter=',', names=['Mstar', 'R50'])
BDL11Lb = np.genfromtxt('./Obs_Data/BDL11Lb.csv', delimiter=',', names=['Mstar', 'R50'])
BDL11Lt = np.genfromtxt('./Obs_Data/BDL11Lt.csv', delimiter=',', names=['Mstar', 'R50'])
LDR15Lb = np.genfromtxt('./Obs_Data/LDR15Lb.csv', delimiter=',', names=['Mstar', 'R50'])
LDR15La = np.genfromtxt('./Obs_Data/LDR15La.csv', delimiter=',', names=['Mstar', 'R50'])

ax2.plot(np.power(10, BDL11Lt['Mstar']), np.power(10, (BDL11Lt['R50'] - 3)), color='blue', lw=lw, linestyle='dashed')
ax2.plot(np.power(10, BDL11Lb['Mstar']), np.power(10, (BDL11Lb['R50'] - 3)), color='blue', lw=lw, linestyle='dashed')
ax2.plot(np.power(10, BDL11Lm['Mstar']), np.power(10, (BDL11Lm['R50'] - 3)), color='blue', lw=lw)
line1, = ax2.plot(LDR15Lb['Mstar'], LDR15Lb['R50'], color='lime', lw=lw)
line2, = ax2.plot(LDR15La['Mstar'], LDR15La['R50'], color='red', lw=lw)

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


# Create an object to combine median and 1-sigma lines from Baldry+12 in the legend #
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [1.0 * height, 1.0 * height], color='blue', linestyle='dashed')
        l2 = plt.Line2D([x0, y0 + width], [0.5 * height, 0.5 * height], color='blue')
        l3 = plt.Line2D([x0, y0 + width], [0.0 * height, 0.0 * height], color='blue', linestyle='dashed')
        return [l1, l2, l3]


legend4 = ax2.legend([object, line1, line2],
                     [r'$\mathrm{Baldry+12:Fit\; &\; 1\operatorname{-}\sigma}$', r'$\mathrm{Lange+15:Single\; p.l.}$',
                      r'$\mathrm{Lange+15:Double\; p.l.}$'], handler_map={object:AnyObjectHandler()}, frameon=False,
                     loc=4)

legend5 = ax2.legend([medianHWT15], [r'$\mathrm{HWT15: Median}$'], frameon=False, loc=1)

ax2.add_artist(legend3)
ax2.add_artist(legend4)
ax2.add_artist(legend5)

######################################################################################################################################################

# Save the figure #
plt.savefig('SM_Vs_SHMR_LTGs_55-' + date + '.pdf', bbox_inches='tight')
