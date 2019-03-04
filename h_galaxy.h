/** @brief Structure with all the data associated with galaxies (this is not the same as the output!)                 */
struct GALAXY            // Galaxy data
{
	int HeapIndex;
	int GalTreeIndex;
	int NextGalaxy;
	
	#ifdef GALAXYTREE
	int FirstProgGal;
	#endif // GALAXYTREE
	
	int       Type;
	int       HaloNr;
	long long MostBoundID;
	int       SnapNum;
	int       CentralGal;        // Galaxy at the centre of this sub-halo. Own ID for types 0 & 1.
	int       FOFCentralGal;     // Galaxy at centre of FOF group; ie the type 0.
	int       MergerCentralGal;  // Galaxy this galaxy will merge onto.  Own ID for types 0 and 1, unless 1's have little dark matter and already merging to type 0.
	float     CentralMvir;
	// Properties of sub-halo at the last time this galaxy was a central galaxy
	float     Pos[3];
	float     MergCentralPos[3];
	float     Vel[3];
	float     Pos_notupdated[3];
	float     Vel_notupdated[3];
	
	#ifdef HALOPROPERTIES
	float HaloM_Mean200;
	float HaloM_Crit200;
	float HaloM_TopHat;
	float HaloPos[3];
	float HaloVel[3];
	float HaloVelDisp;
	float HaloVmax;
	#endif // HALOPROPERTIES
	
	double HaloSpin[3];
	double ColdGasSpin[3];
	double DiskSpin[3];
	double BulgeSpin[3];
	int    Len;
	float  Mvir;
	float  Rvir;
	float  Vvir;
	float  Vmax;
	float  InfallVmax;
	float  InfallVmaxPeak; // km/s - Max previous Vmax at infall
	int    InfallSnap;
	float  InfallHotGas;
	float  InfallHotGasRadius;
	float  HotRadius;
	double ColdGas;
	
	#ifdef H2_AND_RINGS
	// Double to make higher precisions for multirings
	 double ColdGasRings[RNUM];
	 double H2fraction;
	 double H2fractionRings[RNUM];
	#endif // H2_AND_RINGS
	
	double DiskMass;
	double BulgeMass;
	double PBMass;
	double CBMassMajor;
	double CBMassMinor;
	double CBMass;
	#ifdef H2_AND_RINGS
	double DiskMassRings[RNUM];
	#ifdef RINGS_IN_BULGES
	double BulgeMassRings[RNUM];
	#endif // RINGS_IN_BULGES
	#endif // H2_AND_RINGS
	
	double HotGas;
	double EjectedMass;
	
	#ifdef EXCESS_MASS
	float ExcessMass;
	#endif // EXCESS_MASS
	
	float BlackHoleMass;
	float BlackHoleGas;
	
	#ifdef DETAILED_METALS_AND_MASS_RETURN
	struct metals MetalsColdGas;
	#ifdef H2_AND_RINGS
	struct metals MetalsColdGasRings[RNUM];
	#endif // H2_AND_RINGS
	struct metals MetalsDiskMass;
	struct metals MetalsBulgeMass;
	#ifdef H2_AND_RINGS
	struct metals MetalsDiskMassRings[RNUM];
	#ifdef RINGS_IN_BULGES
	struct metals MetalsBulgeMassRings[RNUM];
	#endif // RINGS_IN_BULGES
	#endif // H2_AND_RINGS
	struct metals MetalsHotGas;
	struct metals MetalsEjectedMass;
	#ifdef EXCESS_MASS
	struct metals MetalsExcessMass;
	#endif // EXCESS_MASS
	#ifdef METALS_SELF
	struct metals MetalsHotGasSelf;
	#endif // METALS_SELF
	#else // DETAILED_METALS_AND_MASS_RETURN
	double MetalsColdGas;
	#ifdef H2_AND_RINGS
	double MetalsColdGasRings[RNUM];
	#endif // H2_AND_RINGS
	double MetalsDiskMass;
	double MetalsBulgeMass;
	#ifdef H2_AND_RINGS
	double MetalsDiskMassRings[RNUM];
	#ifdef RINGS_IN_BULGES
	double MetalsBulgeMassRings[RNUM];
	#endif // RINGS_IN_BULGES
	#endif // H2_AND_RINGS
	double MetalsHotGas;
	double MetalsEjectedMass;
	#ifdef EXCESS_MASS
	float MetalsExcessMass;
	#endif // EXCESS_MASS
	#ifdef METALS_SELF
	float MetalsHotGasSelf;
	#endif // METALS_SELF
	#endif //DETAILED_METALS_AND_MASS_RETURN
	
	#ifdef TRACK_MASSGROWTH_CHANNELS
	float MassFromInSitu;
	float MassFromMergers;
	float MassFromBursts;
	#endif // TRACK_MASSGROWTH_CHANNELS
	
	#ifdef TRACK_BURST
	float BurstMass;
	#endif //TRACK_BURST
	
	// Misc
	float  PrimordialAccretionRate;
	float  CoolingRate;
	float  CoolingRate_beforeAGN;
	float  CoolingRadius;
	double CoolingGas;
	float  QuasarAccretionRate;
	float  RadioAccretionRate;
	float  AGNheatingFromCentral;
	
	#ifdef H2_AND_RINGS
	double Sfr;
	double SfrRings[RNUM];
	#else // H2_AND_RINGS
	float Sfr;
	#endif // H2_AND_RINGS
	
	float SfrBulge;
	float StarMerge;
	float XrayLum;
	float BulgeSize;
	float DiskRadius;
	float ColdGasRadius;
	float StellarHalfMassRadius;
	
	#ifdef GALAXYTREE
	int   DisruptOn;
	#endif // GALAXYTREE
	
	float CosInclination; // Angle between galaxy spin and the z-axis
	
	#ifndef HT09_DISRUPTION
	float OriMergTime;
	float MergTime;
	float OriMvir;
	float OriRvir;
	#else // HT09_DISRUPTION
	float OriMergRadius;
	float MergRadius;
	float OriMergmass;
	#endif // HT09_DISRUPTION
	
	float MergeSat;
	float DistanceToCentralGal[3];
	int   MergeOn;
	
	#ifdef TRACK_SPLASHBACKS
	int flagSplashBack;
	float TimeSinceSplashBack;
	#endif // TRACK_SPLASHBACKS
	float NMajorMergers;
	float NMinorMergers;
	float ICM;
	#ifdef DETAILED_METALS_AND_MASS_RETURN
	struct metals MetalsICM;
	#else // DETAILED_METALS_AND_MASS_RETURN
	float MetalsICM;
	#endif // DETAILED_METALS_AND_MASS_RETURN
	
	// Luminosities in various bands.
	#ifdef COMPUTE_SPECPHOT_PROPERTIES
	#ifndef POST_PROCESS_MAGS
	#ifdef OUTPUT_REST_MAGS
	float Lum[NMAG][NOUT];
	float YLum[NMAG][NOUT];
	float LumBulge[NMAG][NOUT];
	float YLumBulge[NMAG][NOUT];
	float LumDust[NMAG][NOUT];
	#ifdef ICL
	float ICLLum[NMAG][NOUT];
	#endif // ICL
	#endif //OUTPUT_REST_MAGS
	#ifdef COMPUTE_OBS_MAGS
	float ObsLum[NMAG][NOUT];
	float ObsYLum[NMAG][NOUT];
	float ObsLumBulge[NMAG][NOUT];
	float ObsYLumBulge[NMAG][NOUT];
	float ObsLumDust[NMAG][NOUT];
	#ifdef ICL
	float ObsICL[NMAG][NOUT];
	#endif // ICL
	#ifdef OUTPUT_MOMAF_INPUTS
	float dObsLum[NMAG][NOUT];
	float dObsYLum[NMAG][NOUT];
	float dObsLumBulge[NMAG][NOUT];
	float dObsYLumBulge[NMAG][NOUT];
	float dObsLumDust[NMAG][NOUT];
	#ifdef ICL
	float dObsICL[NMAG][NOUT];
	#endif // ICL
	#endif // OUTPUT_MOMAF_INPUTS
	#endif // COMPUTE_OBS_MAGS
	#endif // POST_PROCESS_MAGS
	#endif // COMPUTE_SPECPHOT_PROPERTIES
	
	float MassWeightAge[NOUT];
	#ifdef STAR_FORMATION_HISTORY
	int sfh_ibin; // Index of highest bin currently in use.
	float sfh_age; // Time in years of last call to sph_update_bins.
	float sfh_dt[SFH_NBIN]; // Size of time interval in units of years.
	float sfh_t[SFH_NBIN]; // Time at low-redshift edge of bin in same units.
	int sfh_Nbins[SFH_NBIN]; // Number of bins on the time interval.
	double sfh_DiskMass[SFH_NBIN]; // Stellar mass in disk, in bin in standard units.
	double sfh_BulgeMass[SFH_NBIN]; // Stellar mass in bulge, in bin in standard units.
	#ifdef H2_AND_RINGS
	double sfh_DiskMassRings[RNUM][SFH_NBIN]; // Stellar mass in disk RINGS, in bin in standard units.
	#ifdef RINGS_IN_BULGES
	double sfh_BulgeMassRings[RNUM][SFH_NBIN]; // Stellar mass in disk RINGS, in bin in standard units.
	#endif // RINGS_IN_BULGES
	#endif // H2_AND_RINGS
	float sfh_ICM[SFH_NBIN]; // Stellar mass in ICM, in bin in standard units.
	#ifdef DETAILED_METALS_AND_MASS_RETURN
	struct metals sfh_MetalsDiskMass[SFH_NBIN]; // Metals locked up in stars in disk.
	struct metals sfh_MetalsBulgeMass[SFH_NBIN]; // Metals locked up in stars in bulge.
	struct metals sfh_MetalsICM[SFH_NBIN]; // Metals locked up in stars in ICM.
	#else // DETAILED_METALS_AND_MASS_RETURN
	double sfh_MetalsDiskMass[SFH_NBIN]; // Metals locked up in stars in disk.
	double sfh_MetalsBulgeMass[SFH_NBIN]; // Metals locked up in stars in bulge.
	double sfh_MetalsICM[SFH_NBIN]; // Metals locked up in stars in ICM.
	#endif // DETAILED_METALS_AND_MASS_RETURN
	#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
	double sfh_MassFromInSitu[SFH_NBIN]; // Stellar mass formed in situ, in standard units.
	double sfh_MassFromMergers[SFH_NBIN]; // Stellar mass accreted from mergers, in standard units.
	double sfh_MassFromBursts[SFH_NBIN]; // Stellar mass formed in bursts, in standard units.
	#endif // TRACK_SFH_MASSGROWTH_CHANNELS
	#ifdef TRACK_BURST
	double sfh_BurstMass[SFH_NBIN]; // Stellar mass formed in bursts, in standard units.
	#endif //TRACK_BURST
	#endif //STAR_FORMATION_HISTORY
	
	#ifdef INDIVIDUAL_ELEMENTS
	//All: [H][He][Cb][N][O][Ne][Mg][Si][S][Ca][Fe] or //Only [H][He][O][Mg][Fe]
	#ifdef STAR_FORMATION_HISTORY
	float sfh_DiskMass_elements[SFH_NBIN][NUM_ELEMENTS];
	float sfh_BulgeMass_elements[SFH_NBIN][NUM_ELEMENTS];
	float sfh_ICM_elements[SFH_NBIN][NUM_ELEMENTS];
	#endif // STAR_FORMATION_HISTORY
	float DiskMass_elements[NUM_ELEMENTS];
	float BulgeMass_elements[NUM_ELEMENTS];
	float ColdGas_elements[NUM_ELEMENTS];
	float HotGas_elements[NUM_ELEMENTS];
	float ICM_elements[NUM_ELEMENTS];
	float EjectedMass_elements[NUM_ELEMENTS];
	#ifdef EXCESS_MASS
	float ExcessMass_elements[NUM_ELEMENTS];
	#endif // EXCESS_MASS
	#ifdef H2_AND_RINGS
	float DiskMassRings_elements[RNUM][NUM_ELEMENTS];
	#ifdef RINGS_IN_BULGES
	float BulgeMassRings_elements[RNUM][NUM_ELEMENTS];
	#endif // RINGS_IN_BULGES
	float ColdGasRings_elements[RNUM][NUM_ELEMENTS];
	#endif // H2_AND_RINGS
	#endif //INDIVIDUAL_ELEMENTS
} *Gal, *HaloGal;