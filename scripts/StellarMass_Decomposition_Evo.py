# Load data arrays #
CBMass = np.load(SavePath + '/CBMass.npy')
PBMass = np.load(SavePath + '/PBMass.npy')
DiskMass = np.load(SavePath + '/DiskMass.npy')
BulgeMass = np.load(SavePath + '/BulgeMass.npy')
CBMassMajor = np.load(SavePath + '/CBMassMajor.npy')
CBMassMinor = np.load(SavePath + '/CBMassMinor.npy')
StellarMass = np.load(SavePath + '/StellarMass.npy')

# Trim the data #
CB_Mass = CBMass * MassUnits
PB_Mass = PBMass * MassUnits
Disk_Mass = DiskMass * MassUnits
Bulge_Mass = BulgeMass * MassUnits
Stellar_Mass = StellarMass * MassUnits
CB_Mass_Major = CBMassMajor * MassUnits
CB_Mass_Minor = CBMassMinor * MassUnits

######################################################################################################################################################

# Generate initial figure #
if snap == '58':
    plt.close()
    figure, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(10, 7.5))
    figure.subplots_adjust(hspace=0, wspace=0)

    # Figure parameters #
    ax1.set_xlim(8.9, 12.1)
    ax2.set_xlim(8.9, 12.1)
    ax3.set_xlim(8.9, 12.1)
    ax4.set_xlim(8.9, 12.1)
    ax1.set_ylim(-0.05, 1.05)
    ax2.set_ylim(-0.05, 1.05)
    ax3.set_ylim(-0.05, 1.05)
    ax4.set_ylim(-0.05, 1.05)
    ax1.set_ylabel(r'$\mathrm{M_{comp.}/M_{\bigstar}}$')
    ax2.set_ylabel(r'$\mathrm{M_{comp.}/M_{\bigstar}}$')
    ax3.set_ylabel(r'$\mathrm{M_{comp.}/M_{\bigstar}}$')
    ax4.set_ylabel(r'$\mathrm{M_{comp.}/M_{\bigstar}}$')
    ax1.set_xlabel(r'$\mathrm{log_{10}(M_{\bigstar}/M_{\odot})}$')
    ax2.set_xlabel(r'$\mathrm{log_{10}(M_{\bigstar}/M_{\odot})}$')
    ax3.set_xlabel(r'$\mathrm{log_{10}(M_{\bigstar}/M_{\odot})}$')
    ax4.set_xlabel(r'$\mathrm{log_{10}(M_{\bigstar}/M_{\odot})}$')

    # Change labels's position #
    ax2.yaxis.tick_right()
    ax4.yaxis.tick_right()

    # Change ticks's position #
    ax2.yaxis.set_label_position("right")
    ax4.yaxis.set_label_position("right")

    ax1.tick_params(direction='in', which='both', top='on', right='on')
    ax2.tick_params(direction='in', which='both', top='on', left='on')
    ax3.tick_params(direction='in', which='both', top='on', right='on')
    ax4.tick_params(direction='in', which='both', top='on', left='on')

    # Stellar masses bins #
    x = np.arange(9, 12.0, 0.2)

    ##################################################################################################################################################

    # Split stellar masses into bins and create stacked bars #
    if snap == '58':
        for i in range(0, len(x)):
            if i < len(x) - 1:
                index = np.where((np.log10(Stellar_Mass) > x[i]) & (np.log10(Stellar_Mass) < x[i + 1]))
            else:
                index = np.where(np.log10(Stellar_Mass) > x[i])
            Disk_Ratio = np.divide(Disk_Mass[index], Stellar_Mass[index])
            ID_Ratio = np.divide(PB_Mass[index], Stellar_Mass[index])
            MiMD_Ratio = np.divide(CB_Mass_Minor[index], Stellar_Mass[index])
            MaMD_Ratio = np.divide(CB_Mass_Major[index], Stellar_Mass[index])
            ax1.bar(x[i], np.nanmean(Disk_Ratio), align='edge', width=0.2, color='blue', alpha=0.6, edgecolor='None')
            ax1.bar(x[i], np.nanmean(ID_Ratio), bottom=np.nanmean(Disk_Ratio), align='edge', width=0.2, color='green', alpha=0.6, edgecolor='None')
            ax1.bar(x[i], np.nanmean(MaMD_Ratio), bottom=np.nanmean(ID_Ratio) + np.nanmean(Disk_Ratio), align='edge', width=0.2, color='cyan',
                    alpha=0.6, edgecolor='None')
            ax1.bar(x[i], np.nanmean(MiMD_Ratio), bottom=np.nanmean(ID_Ratio) + np.nanmean(Disk_Ratio) + np.nanmean(MaMD_Ratio), align='edge',
                    width=0.2, color='magenta', alpha=0.6, edgecolor='None')

        ax1.annotate(r'$\mathrm{z \sim 0.0}$', xy=(0.03, 0.1), xycoords='axes fraction', color='white', size=20)

if snap == '38':
    # Split stellar masses into bins and create stacked bars #
    for i in range(0, len(x)):
        if i < len(x) - 1:
            index = np.where((np.log10(Stellar_Mass) > x[i]) & (np.log10(Stellar_Mass) < x[i + 1]))
        else:
            index = np.where(np.log10(Stellar_Mass) > x[i])
        Disk_Ratio = np.divide(Disk_Mass[index], Stellar_Mass[index])
        ID_Ratio = np.divide(PB_Mass[index], Stellar_Mass[index])
        MiMD_Ratio = np.divide(CB_Mass_Minor[index], Stellar_Mass[index])
        MaMD_Ratio = np.divide(CB_Mass_Major[index], Stellar_Mass[index])
        ax2.bar(x[i], np.nanmean(Disk_Ratio), align='edge', width=0.2, color='blue', alpha=0.6, edgecolor='None')
        ax2.bar(x[i], np.nanmean(ID_Ratio), bottom=np.nanmean(Disk_Ratio), align='edge', width=0.2, color='green', alpha=0.6, edgecolor='None')
        ax2.bar(x[i], np.nanmean(MaMD_Ratio), bottom=np.nanmean(ID_Ratio) + np.nanmean(Disk_Ratio), align='edge', width=0.2, color='cyan', alpha=0.6,
                edgecolor='None')
        ax2.bar(x[i], np.nanmean(MiMD_Ratio), bottom=np.nanmean(ID_Ratio) + np.nanmean(Disk_Ratio) + np.nanmean(MaMD_Ratio), align='edge', width=0.2,
                color='magenta', alpha=0.6, edgecolor='None')

    ax2.annotate(r'$\mathrm{z \sim 1.0}$', xy=(0.03, 0.1), xycoords='axes fraction', color='white', size=20)

if snap == '30':
    # Split stellar masses into bins and create stacked bars #
    for i in range(0, len(x)):
        if i < len(x) - 1:
            index = np.where((np.log10(Stellar_Mass) > x[i]) & (np.log10(Stellar_Mass) < x[i + 1]))
        else:
            index = np.where(np.log10(Stellar_Mass) > x[i])
        Disk_Ratio = np.divide(Disk_Mass[index], Stellar_Mass[index])
        ID_Ratio = np.divide(PB_Mass[index], Stellar_Mass[index])
        MiMD_Ratio = np.divide(CB_Mass_Minor[index], Stellar_Mass[index])
        MaMD_Ratio = np.divide(CB_Mass_Major[index], Stellar_Mass[index])
        ax3.bar(x[i], np.nanmean(Disk_Ratio), align='edge', width=0.2, color='blue', alpha=0.6, edgecolor='None')
        ax3.bar(x[i], np.nanmean(ID_Ratio), bottom=np.nanmean(Disk_Ratio), align='edge', width=0.2, color='green', alpha=0.6, edgecolor='None')
        ax3.bar(x[i], np.nanmean(MaMD_Ratio), bottom=np.nanmean(ID_Ratio) + np.nanmean(Disk_Ratio), align='edge', width=0.2, color='cyan', alpha=0.6,
                edgecolor='None')
        ax3.bar(x[i], np.nanmean(MiMD_Ratio), bottom=np.nanmean(ID_Ratio) + np.nanmean(Disk_Ratio) + np.nanmean(MaMD_Ratio), align='edge', width=0.2,
                color='magenta', alpha=0.6, edgecolor='None')

    ax3.annotate(r'$\mathrm{z \sim 2.0}$', xy=(0.03, 0.1), xycoords='axes fraction', color='white', size=20)

if snap == '25':
    # Split stellar masses into bins and create stacked bars #
    for i in range(0, len(x)):
        if i < len(x) - 1:
            index = np.where((np.log10(Stellar_Mass) > x[i]) & (np.log10(Stellar_Mass) < x[i + 1]))
        else:
            index = np.where(np.log10(Stellar_Mass) > x[i])
        Disk_Ratio = np.divide(Disk_Mass[index], Stellar_Mass[index])
        ID_Ratio = np.divide(PB_Mass[index], Stellar_Mass[index])
        MiMD_Ratio = np.divide(CB_Mass_Minor[index], Stellar_Mass[index])
        MaMD_Ratio = np.divide(CB_Mass_Major[index], Stellar_Mass[index])
        ax4.bar(x[i], np.nanmean(Disk_Ratio), align='edge', width=0.2, color='blue', alpha=0.6, edgecolor='None')
        ax4.bar(x[i], np.nanmean(ID_Ratio), bottom=np.nanmean(Disk_Ratio), align='edge', width=0.2, color='green', alpha=0.6, edgecolor='None')
        ax4.bar(x[i], np.nanmean(MaMD_Ratio), bottom=np.nanmean(ID_Ratio) + np.nanmean(Disk_Ratio), align='edge', width=0.2, color='cyan', alpha=0.6,
                edgecolor='None')
        ax4.bar(x[i], np.nanmean(MiMD_Ratio), bottom=np.nanmean(ID_Ratio) + np.nanmean(Disk_Ratio) + np.nanmean(MaMD_Ratio), align='edge', width=0.2,
                color='magenta', alpha=0.6, edgecolor='None')

    ax4.annotate(r'$\mathrm{z \sim 3.0}$', xy=(0.03, 0.1), xycoords='axes fraction', color='white', size=20)

    ##################################################################################################################################################

    # Save the figure #
    plt.savefig('CompMassRatio_Vs_SM_Evo-' + date + '.pdf', bbox_inches='tight')