# Load data arrays #
HaloSpin = np.load(SavePath + 'HaloSpin.npy')
ColdGasSpin = np.load(SavePath + 'ColdGasSpin.npy')

index = np.where((np.linalg.norm(ColdGasSpin, axis=1) > 0.0) & (np.linalg.norm(HaloSpin, axis=1) > 0.0))

Halo_Spin = np.linalg.norm(HaloSpin[index], axis=1) * SpinUnits
ColdGas_Spin = np.linalg.norm(ColdGasSpin[index], axis=1) * SpinUnits

dlog10 = 0.1

######################################################################################################################################################

# Generate initial figure #
plt.close()
figure, ax = plt.subplots(1, figsize=(10, 7.5))

# Figure parameters #
plt.xlim(0, 2)
plt.ylabel(r'$\mathrm{Probability\, density\, function}$')
plt.xlabel(r'$\mathrm{|\vec{J}_{d,gas}|/|\vec{J}_{halo}|}$')
plt.tick_params(direction='in', which='both', top='on', right='on')

######################################################################################################################################################

Ratio = np.divide(ColdGas_Spin, Halo_Spin)
weights = np.ones_like(Ratio) / float(len(ColdGas_Spin / Halo_Spin))
plt.hist(Ratio, range=(0, 2), bins=50, color="black", density='True', alpha=0.6, edgecolor='None')

MeanArrow = plt.arrow(np.mean(Ratio), 0.25, 0, -0.25, length_includes_head=True, head_width=0.05, head_length=0.1,
                      ec='red', fc='red')
MedianArrow = plt.arrow(np.median(Ratio), 0.25, 0, -0.25, length_includes_head=True, head_width=0.05, head_length=0.1,
                        ec='blue', fc='blue')


def make_legend_arrow(legend, orig_handle, xdescent, ydescent, width, height, fontsize):
    p = mpatches.FancyArrow(0, 0.5 * height, width, 0, length_includes_head=True, head_width=0.75 * height)
    return p


plt.legend([MeanArrow, MedianArrow], [r'$\mathrm{Mean}$', r'$\mathrm{Median}$'],
           handler_map={mpatches.FancyArrow:HandlerPatch(patch_func=make_legend_arrow)}, frameon=False, loc=1)

######################################################################################################################################################

# Create a zoomed-up inset plot #
axins = inset_axes(ax, width="30%", height="20%", loc=7)
axins.set_xlim(0.721, 0.749)
axins.set_ylim(0.0, 0.25)
axins.set_xticks([0.73, 0.74])
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.1")

weights = np.ones_like(Ratio) / float(len(ColdGas_Spin / Halo_Spin))
axins.hist(Ratio, range=(0, 2), bins=50, color="black", density='True', alpha=0.6, edgecolor='None')

MeanArrow = axins.arrow(np.mean(Ratio), 0.25, 0, -0.25, length_includes_head=True, head_width=0.005, head_length=0.1,
                        ec='red', fc='red')
MedianArrow = axins.arrow(np.median(Ratio), 0.25, 0, -0.25, length_includes_head=True, head_width=0.005,
                          head_length=0.1, ec='blue', fc='blue')

######################################################################################################################################################

# Save the figure #
plt.savefig('PDF-58' + date + '.pdf', bbox_inches='tight')
