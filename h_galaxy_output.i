# 1 "./code/h_galaxy_output.h"
# 1 "/lustre/scratch/astro/di43/Development_Branch//"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "./code/h_galaxy_output.h"
/**
 * Galaxy structure for output
 */
# 41 "./code/h_galaxy_output.h"
#pragma pack(1)
struct GALAXY_OUTPUT {
# 72 "./code/h_galaxy_output.h"
	int Type; // //Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.

	int HaloIndex; // // ?Unique ID of MPA halo containing this galaxy


	float HaloM_Mean200; // 1e10 Msun/h // M200 cf mean last time this halo was a type 0
	float HaloM_Crit200; // 1e10 Msun/h // M200 cf critical last time this halo was a type 0
	float HaloM_TopHat; // 1e10 Msun/h // Virial mass last time this halo was a type 0
	float HaloPos[3]; // Mpc/h // Comoving position of halo.
	float HaloVel[3]; // km/s // Mean velocity of halo.
	float HaloVelDisp; // km/s // Velocity dispersion of halo.
	float HaloVmax; // km/s // Maximum circular velocity of halo.
	float HaloSpin[3]; // km/s Mpc/h // specific spin of the halo.

	int   SnapNum; // //The snapshot number where this galaxy was identified.
	float LookBackTimeToSnap; //yr //The time from a given snapshot to z=0
	float CentralMvir; // 10^10/h Msun // virial mass of background (FOF) halo containing this galaxy
	float CentralRvir; // Mpc/h // Proper[?] R200 cf critical of background (FOF) halo containing this galaxy
	float DistanceToCentralGal[3]; // Mpc/h // Proper[?] components of the distance between this galaxy and the galaxy at the centre of the FoF group.
	float Pos[3]; // 1/h Mpc // Comoving galaxy/subhalo position
	float Vel[3]; // km/s // Galaxy/subhalo peculiar velocity
	int   Len; // // Number of particles in the associated subhalo
	/* properties of subhalo at the last time this galaxy was a central galaxy */
	float Mvir; // 10^10/h Msun // M200 cf critical of the halo last time galaxy was type 0
	float Rvir; // Mpc/h // R200 cf critical of the halo last time galaxy was type 0
	float Vvir; // km/s // Virial velocity of the halo last time galaxy was type 0
	float Vmax; // km/s // Maximum rotational velocity of the subhalo, or the last value for type 2's galaxies.
	float ColdGasSpin[3]; // Mpc/h km/s // The specific angular momentum of the cold gas disk
	float DiskSpin[3]; // Mpc/h km/s // The specific angular momentum of the stellar disk
	float BulgeSpin[3]; // Mpc/h km/s // The specific angular momentum of the bulge
	float InfallVmax; // km/s // Maximum rotational velocity of the host halo of this galaxy at infall (ie last time a type 0)
	float InfallVmaxPeak; // km/s // ? Peak Vmax along past history
	int   InfallSnap; // // Most recent (largest) snapnum at which this galaxy's type changed from 0 to 1 or 2
	float InfallHotGas; // 10^10 Msun/h // Mass in hot gas at the time of infall (same as hotGas for type 0 galaxies).
	float
	      HotRadius; // Mpc/h // Proper[?] radius out to which hot gas extends: rvir for type 0; 0 for type 2; maximum radius out to which hot gas is not stripped for type 1.
	/*dynamical friction merger time*/

	float OriMergTime; // yr // Estimated dynamical friction time when the merger clock is set.
	float MergTime; //yr // Estimated remaining merging time.
# 119 "./code/h_galaxy_output.h"
	float NMajorMergers;
	float NMinorMergers;
	/* baryonic reservoirs */
	float ColdGas; // 10^10/h Msun // Mass in cold gas.





	float StellarMass; // 10^10/h Msun // Total mass in stars in the disk and the bulge combined
	float DiskMass; // 10^10/h Msun // Mass of stars in the disk
	float BulgeMass; // 10^10/h Msun // Mass of stars in the bulge
	float PBMass; // 10^10/h Msun // Mass of stars in the bulge due to DI
	float CBMassMajor; // 10^10/h Msun // Mass of stars in the bulge due to major mergers
	float CBMassMinor; // 10^10/h Msun // Mass of stars in the bulge due to minor mergers
	float CBMass; // 10^10/h Msun // Mass of stars in the bulge due to mergers.






	float HotGas; // 10^10/h Msun // Mass in hot gas
	float ReheatedGas; // 10^10/h Msun // Mass in reheated gas
	float EjectedMass; // 10^10/h Msun // Mass in ejected gas



	float BlackHoleMass; // 10^10/h Msun // Mass of central black hole
	//float BlackHoleGas; // 10^10/h Msun // Mass in BH accretion disk
	/* ICL magnitude and mass*/
	float ICM; //10^10/h Msun //Total mass in metals in intra-cluster stars, for type 0,1
# 184 "./code/h_galaxy_output.h"
	float MetalsColdGas; // 10^10/h Msun // Mass in metals in cold gas.



	float MetalsStellarMass; // 10^10/h Msun //	Mass in metals in the bulge+disk
	float MetalsDiskMass; // 10^10/h Msun // Mass in metals in the disk
	float MetalsBulgeMass; // 10^10/h Msun // Mass in metals in the bulge






	float MetalsHotGas; // 10^10/h Msun // Mass in metals in the hot gas
	//float MetalsReheatedGas; // 10^10/h Msun // Mass in metals in the Reheated gas
	float MetalsEjectedMass; // 10^10/h Msun // Mass in metals in the ejected gas



	float MetalsICM; // 10^10/h Msun // Mass in metals in intra-cluster stars, for type 0,1




	/* misc */
	float PrimordialAccretionRate; // Msun/yr // Accretion rate of primordial gas.
	float CoolingRadius; // Mpc/h // The radius within which the cooling time scale is shorter than the dynamical timescale
	//float CoolingGas; // Mpc/h // Mass of cooling gas
	float CoolingRate; // Msun/yr // Cooling rate of the hot gas
	float CoolingRate_beforeAGN; // Msun/yr // What the cooling rate of the hot gas would have been if there was no AGN feedback.
	float QuasarAccretionRate; // Msun/yr // Rate at which cold gas is accreted into the central black hole in the quasar mode.
	float RadioAccretionRate; // Msun/yr // Rate at which hot gas is accreted into the central black hole in the radio mode.
	float Sfr; // Msun/yr // Star formation rate



	float SfrBulge; // Msun/yr // Star formation rate in bulge.
	float XrayLum; // log10(erg/sec) // (log_10 of) X-Ray luminosity
	float BulgeSize; // Mpc/h // Half mass radius of bulge
	float DiskRadius; // Mpc/h // Size of the stellar disk, 3x the scale length.
	float ColdGasRadius; // Mpc/h // Size of the gas disk, 3x the scale length.
	float StellarHalfMassRadius; // Mpc/h // stellar Half mass radius
	float StellarHalfLightRadius; // Mpc/h // stellar Half light radius
	float StellarNinetyLightRadius; // Mpc/h // stellar Half light radius
	float CosInclination; // deg // Inclination of the galaxy. Derived from the angle between the stellar spins of the galaxy and the z-axis
	int
	      DisruptOn; // // 0: galaxy merged onto merger center 1: galaxy was disrupted before merging onto its descendant, matter went into ICM of merger center;

	int
			MergeOn; // // 0: standard delucia-like merger behaviour for type 1 galaxy; 1: galaxy mass > halo mass, separate dynamical friction time calculated ....



	/* magnitudes in various bands */


	float MagDust[40]; // // dust corrected, rest-frame absolute mags
	float Mag[40]; // // rest-frame absolute mags
	float MagBulge[40]; // // rest-frame absolute mags for the bulge
# 273 "./code/h_galaxy_output.h"
	float MassWeightAge; //10^9yr //The age of this galaxy weighted by mass of its components.
# 339 "./code/h_galaxy_output.h"
};

// next only of interest to DB output, which generally requires complete tree
# 386 "./code/h_galaxy_output.h"
#pragma pack()
