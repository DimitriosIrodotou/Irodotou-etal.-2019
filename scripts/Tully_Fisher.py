# Load data arrays #
sfr = np.load(SavePath + 'Sfr.npy')
Vmax = np.load(SavePath + 'Vmax.npy')
Type = np.load(SavePath + 'Type.npy')
ColdGas = np.load(SavePath + 'ColdGas.npy')
DiskMass = np.load(SavePath + 'DiskMass.npy')
StellarMass = np.load(SavePath + 'StellarMass.npy')

dlog10 = 0.02

######################################################################################################################################################

# Generate initial figure #
if snap == '58':
    plt.close()
    figure, ax = plt.subplots()
    figure, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(10, 15))
    figure.subplots_adjust(hspace=0.2)

    # Figure parameters #
    ax1.set_ylim(9.01, 12.1)
    ax2.set_ylim(9.01, 12.1)
    ax1.set_xlim(1.6, 3.0)
    ax2.set_xlim(1.6, 3.0)

    ax1.set_xlabel(r'$\mathrm{log_{10}(V_{c} / (km \cdot s^{-1}))}$', size=25)
    ax2.set_xlabel(r'$\mathrm{log_{10}(V_{c} / (km \cdot s^{-1}))}$', size=25)
    ax1.set_ylabel(r'$\mathrm{log_{10}((M_{\bigstar} + M_{d,gas}) / M_{\odot})}$', size=25)
    ax2.set_ylabel(r'$\mathrm{log_{10}((M_{\bigstar} + M_{d,gas}) / M_{\odot})}$', size=25)

    ax1.tick_params(direction='in', which='both', top='on', right='on', labelsize=25)
    ax2.tick_params(direction='in', which='both', top='on', right='on', labelsize=25)

    ##################################################################################################################################################

    # Trim the data #
    index = np.where((ColdGas > StellarMass) & (Type == 0))

    V_max = Vmax[index] * VelocityUnits
    Baryons = (ColdGas[index] + StellarMass[index]) * MassUnits

    # Read observational data from McGaugh12 #
    McGaugh12 = np.genfromtxt('./Obs_Data/McGaugh12.csv', delimiter=',', names=['x', 's', 'g'])

    # Plot observational data from McGaugh12 #
    ax1.scatter(np.log10(McGaugh12['x']), np.log10(np.power(10, McGaugh12['s']) + np.power(10, McGaugh12['g'])),
                edgecolor='black', color='red', s=50, marker='^', label=r'$\mathrm{McGaugh12}$', zorder=2)

    # Plot L-Galaxies data - 2D histogram #
    h = ax1.hexbin(np.log10(V_max), np.log10(Baryons), bins='log', cmap='Greys', mincnt=2)

    # Adjust the color bar #
    cbaxes = figure.add_axes([0.9, 0.53, 0.02, 0.35])
    cb = plt.colorbar(h, cax=cbaxes)
    cbaxes.tick_params(direction='out', which='both', right='on', labelsize=25)
    cb.set_label('$\mathrm{Counts\; per\; hexbin}$', size=25)

    # Create the legends #
    colors = ['black', 'grey', 'lightgrey']
    squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
    legend1 = ax1.legend([squares], [r'$\mathrm{This\; work:M_{d,gas} > M_{\bigstar}}$'], scatterpoints=len(colors),
                         scatteryoffsets=[.5], handlelength=len(colors), markerscale=2, frameon=False, loc=2)

    ax1.add_artist(legend1)
    # ax1.add_artist(legend2)
    ax1.legend(frameon=False, scatterpoints=3, loc=4)

    ##################################################################################################################################################

    # Trim the data #
    index = np.where((DiskMass > 0.7 * StellarMass) & (Type == 0))

    V_max = Vmax[index] * VelocityUnits
    Baryons = (ColdGas[index] + StellarMass[index]) * MassUnits

    # Read observational data from AZF08, TEA11 and AVS18 #
    AZF08 = np.genfromtxt('./Obs_Data/AZF08.csv', delimiter=',', names=['x', 'y'])
    TEA11 = np.genfromtxt('./Obs_Data/TEA11.csv', delimiter=',', names=['x', 'y'])
    AVS18 = np.genfromtxt('./Obs_Data/AVS18.csv', delimiter=',', names=['x', 'y'])

    # Plot observational data from AZF08, AVS18 #
    ax2.scatter(AZF08['y'], AZF08['x'], edgecolor='black', color='lime', s=50, zorder=2,
                label=r'$\mathrm{Avila-Reese+08}$')
    ax2.scatter(TEA11['x'], TEA11['y'], edgecolor='black', color='blue', s=50, marker='s', zorder=2,
                label=r'$\mathrm{Torres-Flores+11}$')

    # Plot L-Galaxies data - 2D histogram #
    hexbin = ax2.hexbin(np.log10(V_max), np.log10(Baryons), bins='log', cmap='Greys', mincnt=2)
    # Adjust the color bar #
    cbaxes = figure.add_axes([0.9, 0.11, 0.02, 0.35])
    cb = plt.colorbar(hexbin, cax=cbaxes)
    cbaxes.tick_params(direction='out', which='both', right='on', labelsize=25)
    cb.set_label('$\mathrm{Counts\; per\; hexbin}$', size=25)

    # Create the legends #
    colors = ['black', 'grey', 'lightgrey']
    squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
    legend3 = ax2.legend([squares], [r'$\mathrm{This\; work:M_{d,\bigstar}/M_{\bigstar}>0.7}$'],
                         scatterpoints=len(colors), scatteryoffsets=[.5], handlelength=len(colors), markerscale=2,
                         frameon=False, loc=2)

    ax2.add_artist(legend3)
    ax2.legend(frameon=False, scatterpoints=2, loc=4)

    # Save the figure #
    plt.savefig('TF_58-' + date + '.pdf', bbox_inches='tight')

######################################################################################################################################################
#
# if snap == '58':
#     # Generate initial figure #
#     plt.close()
#     figure, ax = plt.subplots()
#     figure, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharex=True, figsize=(30, 7.5))
#     figure.subplots_adjust(hspace=0, wspace=0.35)
#
#     # 2D histogram plot parameters #
#     ax1.set_ylim(9.01, 12.1)
#     ax2.set_ylim(9.01, 12.1)
#     ax3.set_ylim(9.01, 12.1)
#     ax1.set_xlim(1.6, 2.8)
#     ax2.set_xlim(1.6, 2.8)
#     ax3.set_xlim(1.6, 2.8)
#
#     ax1.set_xlabel(r'$\mathrm{log_{10}(V_{c} / (km \cdot s^{-1}))}$')
#     ax2.set_xlabel(r'$\mathrm{log_{10}(V_{c} / (km \cdot s^{-1}))}$')
#     ax3.set_xlabel(r'$\mathrm{log_{10}(V_{c} / (km \cdot s^{-1}))}$')
#     ax1.set_ylabel(r'$\mathrm{log_{10}(M_{\bigstar} / M_{\odot})}$')
#     ax2.set_ylabel(r'$\mathrm{log_{10}(M_{\bigstar} / M_{\odot})}$')
#     ax3.set_ylabel(r'$\mathrm{log_{10}(M_{\bigstar} / M_{\odot})}$')
#
#     ax1.tick_params(direction='in', which='both', top='on', right='on')
#     ax2.tick_params(direction='in', which='both', top='on', right='on')
#     ax3.tick_params(direction='in', which='both', top='on', right='on')
#
#     ##################################################################################################################################################
#
#     # Trim the data #
#     index = np.where((DiskMass > 0.7 * StellarMass) & (Type == 0))
#
#     V_max = Vmax[index] * VelocityUnits
#     Stellar_Mass = StellarMass[index] * MassUnits
#
#     # Read observational data from AZF08, TEA11 and AVS18 #
#     CBE05L07T = np.genfromtxt('./Obs_Data/CBE05L07T.csv', delimiter=',', names=['x', 'y'])
#     CBE05L07M = np.genfromtxt('./Obs_Data/CBE05L07M.csv', delimiter=',', names=['x', 'y'])
#     CBE05L07B = np.genfromtxt('./Obs_Data/CBE05L07B.csv', delimiter=',', names=['x', 'y'])
#
#     # Plot observational data from AZF08, AVS18 #
#     ax1.plot(CBE05L07T['x'], CBE05L07T['y'], color='green', zorder=2, label=r'$\mathrm{Conselice+05}$',
#     linestyle='dashed')
#     ax1.plot(CBE05L07M['x'], CBE05L07M['y'], color='green', zorder=2, label=r'$\mathrm{Conselice+05}$')
#     ax1.plot(CBE05L07B['x'], CBE05L07B['y'], color='green', zorder=2, label=r'$\mathrm{Conselice+05}$',
#     linestyle='dashed')
#
#     # Plot L-Galaxies data - 2D histogram #
#     h = ax1.hexbin(np.log10(V_max), np.log10(Stellar_Mass), bins='log', cmap='Greys', gridsize=gs, mincnt=mc)
#
#     # Adjust the color bar #
#     cbaxes = figure.add_axes([0.334, 0.11, 0.01, 0.77])
#     cb = plt.colorbar(h, cax=cbaxes)
#     cb.set_label('$\mathrm{Counts\; per\; hexbin}$')
#
#     # Create the legends #
#     colors = ['black', 'grey', 'lightgrey']
#     squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
#     legend1 = ax1.legend([squares], [r'$\mathrm{This\; work:M_{d,gas} > M_{\bigstar}}$'], scatterpoints=len(
#     colors), scatteryoffsets=[.5],
#                          handlelength=len(colors), markerscale=2, frameon=False, loc=2)
#
#     ax1.add_artist(legend1)
#     ax1.legend(frameon=False, scatterpoints=sp, loc=4)
#
#     ax1.annotate(r'$\mathrm{z \sim 0.0}$', xy=(0.036, 0.85), xycoords='axes fraction', color='black', size=18)
#
# ######################################################################################################################################################
#
# if snap == '38':
#     # Trim the data #
#     index = np.where((DiskMass > 0.7 * StellarMass) & (Type == 0))
#
#     V_max = Vmax[index] * VelocityUnits
#     Stellar_Mass = StellarMass[index] * MassUnits
#
#     # Read observational data from AZF08, TEA11 and AVS18 #
#     CBE05B07 = np.genfromtxt('./Obs_Data/CBE05B07.csv', delimiter=',', names=['x', 'y'])
#
#     # Plot observational data from AZF08, AVS18 #
#     ax2.scatter(CBE05B07['x'], CBE05B07['y'], color='green', zorder=2, label=r'$\mathrm{Conselice+05}$')
#
#     # Plot L-Galaxies data - 2D histogram #
#     hexbin = ax2.hexbin(np.log10(V_max), np.log10(Stellar_Mass), bins='log', cmap='Greys', gridsize=gs, mincnt=mc)
#
#     # Adjust the color bar #
#     cbaxes = figure.add_axes([0.617, 0.11, 0.01, 0.77])
#     cb = plt.colorbar(hexbin, cax=cbaxes)
#     cb.set_label('$\mathrm{Counts\; per\; hexbin}$')
#
#     # Create the legends #
#     colors = ['black', 'grey', 'lightgrey']
#     squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
#     legend3 = ax2.legend([squares], [r'$\mathrm{This\; work:M_{d,\bigstar}/M_{\bigstar}>0.7}$'], scatterpoints=len(
#     colors), scatteryoffsets=[.5],
#                          handlelength=len(colors), markerscale=2, frameon=False, loc=2)
#
#     ax2.add_artist(legend3)
#     ax2.legend(frameon=False, scatterpoints=sp, loc=4)
#
#     ax2.annotate(r'$\mathrm{z \sim 1.0}$', xy=(0.036, 0.85), xycoords='axes fraction', color='black', size=18)
#
# ######################################################################################################################################################
#
# if snap == '30':
#     # Trim the data #
#     index = np.where((sfr > 10) & (Type == 0))
#
#     V_max = Vmax[index] * VelocityUnits
#     Stellar_Mass = StellarMass[index] * MassUnits
#
#     # Read observational data from McGaugh12 #
#     CHG09 = np.loadtxt('./Obs_Data/CHG09.txt')
#
#     # Plot observational data from McGaugh12 #
#     ax3.scatter(np.log10(CHG09[:, 1]), np.log10(CHG09[:, 0] * 1e10), color='green', s=size, marker='^',
#     label=r'$\mathrm{Cresci+09}$', zorder=2)
#
#     # Plot L-Galaxies data - 2D histogram #
#     h = ax3.hexbin(np.log10(V_max), np.log10(Stellar_Mass), bins='log', cmap='Greys', gridsize=gs, mincnt=mc)
#
#     # Adjust the color bar #
#     cbaxes = figure.add_axes([0.9, 0.11, 0.01, 0.77])
#     cb = plt.colorbar(h, cax=cbaxes)
#     cb.set_label('$\mathrm{Counts\; per\; hexbin}$')
#
#     # Create the legends #
#     colors = ['black', 'grey', 'lightgrey']
#     squares = collections.RegularPolyCollection(numsides=6, sizes=(20,), facecolors=colors)
#     legend3 = ax3.legend([squares], [r'$\mathrm{This\; work:M_{d,gas} > M_{\bigstar}}$'], scatterpoints=len(
#     colors), scatteryoffsets=[.5],
#                          handlelength=len(colors), markerscale=2, frameon=False, loc=2)
#
#     ax3.add_artist(legend3)
#     ax3.legend(frameon=False, scatterpoints=sp, loc=4)
#
#     ax3.annotate(r'$\mathrm{z \sim 2.0}$', xy=(0.036, 0.85), xycoords='axes fraction', color='black', size=18)
#
#     # Save the figure #
#     plt.savefig('TF_30-' + date + '.png', bbox_inches='tight')
