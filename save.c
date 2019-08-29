#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"


/**@file save.c
 * @brief Copies the relevant properties in Galaxy structure into Galaxy_Output structure and saves them into the output
 * files (SA_z**_**) - redshift/filenr.
 *
 * There are two distinct procedures to write the output depending on whether GALAXY_TREE option is turned on or off.
 * If it is on the full galaxy tree is written using save_galaxy_tree. If it is off, the output is only written for the
 * chosen output snap numbers using save_galaxies.
 *
 * If USE_MEMORY_TO_MINIMIZE_IO ON, these routines copy the galaxy data from the working structures into pointers until
 * that has been done for all the tree in a given file.
 *
 * After all the galaxy trees are written finalize_galaxy_file is called in main.c to include an header in the output
 * files. If GALAXY_TREE=1 three numbers are written: 1 (int); size_of_struct(Galaxy_Output) (int);UnitTime_in_years
 * and TotGalCount(int). If GALAXY_TREE=0 then the header is: Ntrees (int); total number of galaxies for the snapshot
 * corresponding to the current file -> TotGalaxies[n] (int); and the number of galaxies on each tree on the current
 * snapshot -> TreeNgals[n] (int*Ntrees).
 *
 * If USE_MEMORY_TO_MINIMIZE_IO ON, finalize_galaxy_file also copies the galaxy data stored in pointers into the output
 * files, so that it is all done in one go for a given file. This is done using either write_galaxy_data_snap (for SNAP
 * output), write_all_galaxy_data (for GALAXYTREE option) or write_galaxy_for_momaf (for MOMAF option).
 *
 * If UPDATETYPE2 is defined, the positions of type 2 galaxies (satellites without a dark matter halo) will be updated
 * before output to contain the subsequent positions of the most bound dark matter particle at disruption time (using
 * get_coordinates).                                                                                                  */



















/** @brief Create output files - snapshot option																	  */
void create_galaxy_files (int filenr) {
	int n, i;

#ifdef HDF5_OUTPUT
	open_hdf5_file(filenr);
#else //HDF5_OUTPUT
	// create output files - snapshot option
	char buf[1000];
#endif //HDF5_OUTPUT

	for (n = 0; n < NOUT; n ++) {
		for (i = 0; i < Ntrees; i ++) {
			TreeNgals[n][i] = 0;
		}
#ifdef HDF5_OUTPUT
		create_hdf5_table(n);
#else // HDF5_OUTPUT
		sprintf (buf, "%s/%s_z%1.2f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[n]], filenr);
		if (! (FdGalDumps[n] = fopen (buf, "w+"))) {
			char sbuf[1000];
			sprintf (sbuf, "can't open file `%s'\n", buf);
			terminate(sbuf);
		}
		fseek (FdGalDumps[n], (2 + Ntrees) * sizeof (int), SEEK_SET);
		TotGalaxies[n] = 0;
#endif //HDF5_OUTPUT
	}
}
// end create_galaxy_files



















/** @brief 																											  */
void close_galaxy_files (void) {
	int n;

#ifdef HDF5_OUTPUT
	for ( n = 0; n < NOUT; n++ ) {
		hdf5_append_data(n, galaxy_output_hdf5[n], b[n]);
	} // Output the final galaxies
	hdf5_close();
#else //HDF5_OUTPUT
	for (n = 0; n < NOUT; n ++) {
		fseek (FdGalDumps[n], 0, SEEK_SET);
		myfwrite (&Ntrees, sizeof (int), 1, FdGalDumps[n]);    //Number of trees
		myfwrite (&TotGalaxies[n], sizeof (int), 1, FdGalDumps[n]);    // total number of galaxies
		myfwrite (TreeNgals[n], sizeof (int), Ntrees, FdGalDumps[n]);    // Number of galaxies in each tree
		fclose (FdGalDumps[n]);
	}
#endif //HDF5_OUTPUT
}
// end close_galaxy_files



















/** @brief Saves the Galaxy_Output structure for all the galaxies in the current tree into the current output file (one
 * for each input dark matter file) for the chosen snapshots.
 *
 * If UPDATETYPETWO=1 then the positions and velocities of type 2 galaxies are updated from the most bound dark matter
 * particle. After that the GALAXY_OUTPUT structure is created and written. input: int file number (current file where
 * the output is being written), int tree number (tree being currently treated).
 *
 * If USE_MEMORY_TO_MINIMIZE_IO ON, this write statements in this routine copy the galaxy data from the working
 * structures into pointers until that has been done for all the tree in a given file.								  */
void save_galaxy_append (int tree, int i, int n) {
	struct GALAXY_OUTPUT galaxy_output;

#ifndef NORMALIZEDDB
	prepare_galaxy_for_output (n, &HaloGal[i], &galaxy_output);

#ifdef HDF5_OUTPUT
	if ( b[n] < NRECORDS_APP ) {
		galaxy_output_hdf5[n][b[n]] = galaxy_output;
		b[n]++;
	}
	else {
		// Append the data to the HDF5 table if b[n]==NRECORDS_APP
		hdf5_append_data(n, galaxy_output_hdf5[n], NRECORDS_APP);
		b[n] = 0;
	}
#else //HDF5_OUTPUT
	myfwrite (&galaxy_output, sizeof (struct GALAXY_OUTPUT), 1, FdGalDumps[n]);
#endif //HDF5_OUTPUT
#endif //NORMALIZEDDB

#ifndef HDF5_OUTPUT
	TotGalaxies[n] ++;      // Total number of galaxies
	TreeNgals[n][tree] ++;  // Number of galaxies in each tree
#endif //HDF5_OUTPUT
}
// end save_galaxy_append



















/** @brief Copies all the relevant properties from the Galaxy structure into the Galaxy output structure, some units
 * are corrected.																									  */
#ifdef NORMALIZEDDB
void prepare_galaxy_for_output(int n, struct GALAXY *g, struct GALAXY_OUTPUT *o, struct SFH_BIN *sfh_bin)
#else // NORMALIZEDDB

void prepare_galaxy_for_output (int n, struct GALAXY*g, struct GALAXY_OUTPUT*o)
#endif // NORMALIZEDDB
{
	int j;
#ifndef GALAXYTREE
#ifdef OUTPUT_RINGS
	int jj;
#endif
#ifdef OUTPUT_ELEMENTS
	int kk;
#endif
#endif

#ifndef NO_PROPS_OUTPUTS
	o->Type        = g->Type;
	o->SnapNum     = g->SnapNum;
	o->CentralMvir = get_virial_mass (Halo[g->HaloNr].FirstHaloInFOFgroup);
	o->CentralRvir = get_virial_radius (Halo[g->HaloNr].FirstHaloInFOFgroup);
	o->Mvir        = g->Mvir;
	o->Rvir        = g->Rvir;
	o->Vvir        = g->Vvir;

	for (j = 0; j < 3; j ++) {
		o->Pos[j]                  = g->Pos[j];
		o->DistanceToCentralGal[j] = wrap(Halo[Halo[g->HaloNr].FirstHaloInFOFgroup].Pos[j] - g->Pos[j], BoxSize);
	}

	o->ColdGas     = g->ColdGas;
	o->StellarMass = g->BulgeMass + g->DiskMass;

	o->DiskMass      = g->DiskMass;
	o->BulgeMass     = g->BulgeMass;
	o->PBMass        = g->PBMass;
	o->CBMassMajor   = g->CBMassMajor;
	o->CBMassMinor   = g->CBMassMinor;
	o->CBMass        = g->CBMass;
	o->HotGas        = g->HotGas;
	o->BlackHoleMass = g->BlackHoleMass;
#ifdef OUTPUT_RINGS
	o->H2fraction                 = g->H2fraction;
#endif // OUTPUT_RINGS
#endif // NO_PROPS_OUTPUTS

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
#ifdef OUTPUT_REST_MAGS
	// Luminosities are converted into Mags in various bands.
	for ( j = 0; j < NMAG; j++ ) {
		o->Mag[j] = lum_to_mag(g->Lum[j][n]);
	}
#endif // OUTPUT_REST_MAGS
#endif // POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

#ifndef LIGHT_OUTPUT
#ifndef NO_PROPS_OUTPUTS
#ifdef GALAXYTREE
	o->HaloID   = HaloIDs[g->HaloNr].HaloID;
	o->Redshift = ZZ[g->SnapNum];

	int ii = (int) floor(o->Pos[0] * ScaleFactor);
	int jj = (int) floor(o->Pos[1] * ScaleFactor);
	int kk = (int) floor(o->Pos[2] * ScaleFactor);

	o->PeanoKey = peano_hilbert_key(ii, jj, kk, Hashbits);

	o->SubID = calc_big_db_subid_index(g->SnapNum, Halo[g->HaloNr].FileNr, Halo[g->HaloNr].SubhaloIndex);

#ifdef NIFTY
	o->FOFCentralID=HaloIDs[g->HaloNr].FirstHaloInFOFgroup;
#endif // NIFTY

	int tmpfirst = Halo[g->HaloNr].FirstHaloInFOFgroup;
	int lenmax   = 0;
	int next     = tmpfirst;
	while (next != -1) {
		if ( Halo[next].Len > lenmax ) {
			lenmax   = Halo[next].Len;
			tmpfirst = next;
		}
		next = Halo[next].NextHaloInFOFgroup;
	}

	o->MMSubID = calc_big_db_subid_index(g->SnapNum, Halo[tmpfirst].FileNr, Halo[tmpfirst].SubhaloIndex);
#endif // GALAXYTREE

	o->LookBackTimeToSnap = NumToTime (g->SnapNum) * UnitTime_in_years / Hubble_h;
	o->InfallVmax         = g->InfallVmax;
	o->InfallVmaxPeak     = g->InfallVmaxPeak;
	o->InfallSnap         = g->InfallSnap;
	o->InfallHotGas       = g->InfallHotGas;
	o->HotRadius          = g->HotRadius;
#ifdef HALOPROPERTIES
	o->HaloM_Mean200      = g->HaloM_Mean200;
	o->HaloM_Crit200      = g->HaloM_Crit200;
	o->HaloM_TopHat       = g->HaloM_TopHat;
	o->HaloVelDisp        = g->HaloVelDisp;
	o->HaloVmax           = g->HaloVmax;
#endif // HALOPROPERTIES

	o->Len         = g->Len;
	o->Vmax        = g->Vmax;
	o->EjectedMass = CORRECTDBFLOAT(g->EjectedMass);

	for (j = 0; j < 3; j ++) {
		o->Vel[j] = g->Vel[j];
#ifdef HALOSPIN
		o->HaloSpin[j] = g->HaloSpin[j];
#endif // HALOSPIN
		o->ColdGasSpin[j] = g->ColdGasSpin[j];
		o->DiskSpin[j]    = g->DiskSpin[j];
		o->BulgeSpin[j]   = g->BulgeSpin[j];

#ifdef HALOPROPERTIES
		o->HaloPos[j] = g->HaloPos[j];
		o->HaloVel[j] = g->HaloVel[j];
		o->HaloSpin[j] = g->HaloSpin[j];
#endif // HALOPROPERTIES
	}

	o->XrayLum       = g->XrayLum;
	o->ColdGasRadius = g->ColdGasRadius;
	o->DiskRadius    = g->DiskRadius;

#ifndef RINGS_IN_BULGES
	o->BulgeSize = g->BulgeSize;
#else // RINGS_IN_BULGES
	if ( g->BulgeSize < 1.0e-6 ) {
		o->BulgeSize = RingRadius[0];
	}
	o->BulgeSize = 0.5 * RingRadius[0] * g->BulgeMassRings[0];
	for ( jj = 1; jj < RNUM; jj ++ ) {
		o->BulgeSize += (0.5 * (RingRadius[jj - 1] + RingRadius[jj]) * g->BulgeMassRings[jj]);
	}
	o->BulgeSize = 3. * o->BulgeSize / g->BulgeMass / 2.;
#endif //OUTPUT_RINGS

	o->StellarHalfMassRadius = g->StellarHalfMassRadius;
	o->CoolingRadius         = g->CoolingRadius;
	o->ICM                   = g->ICM;
	o->MetalsICM             = g->MetalsICM;
	o->QuasarAccretionRate   = g->QuasarAccretionRate * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
	o->RadioAccretionRate    = g->RadioAccretionRate * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
	o->CosInclination        = g->CosInclination;

#ifndef HT09_DISRUPTION
	if (g->Type == 2 || (g->Type == 1 && g->MergeOn == 1)) {
		o->OriMergTime = g->OriMergTime * UnitTime_in_years / Hubble_h;
		o->MergTime    = g->MergTime * UnitTime_in_years / Hubble_h;
	} else {
		o->OriMergTime = 0.0;
		o->MergTime    = 0.0;
	}
#else // HT09_DISRUPTION
	if(g->Type == 2 || (g->Type == 1 && g->MergeOn == 1)) {
	   o->OriMergRadius=g->OriMergRadius;
	   o->MergRadius = g->MergRadius;
	 }
	 else {
	   o->OriMergRadius=0.0;
	   o->MergRadius = 0.0;
	 }
#endif // HT09_DISRUPTION

#ifdef TRACK_SPLASHBACKS
	o->flagSplashBack=g->flagSplashBack;
	o->TimeSinceSplashBack=g->TimeSinceSplashBack* UnitTime_in_Megayears / 1000. / Hubble_h;
#endif // TRACK_SPLASHBACKS

#ifndef GALAXYTREE
	o->HaloIndex = g->HaloNr;
#endif // GALAXYTREE
#ifdef MBPID
	o->MostBoundID = g->MostBoundID;
#endif // MBPID

#ifdef GALAXYTREE
	o->DisruptOn = g->DisruptOn;
#endif // GALAXYTREE
#ifdef MERGE01
	o->MergeOn = g->MergeOn;
#endif // MERGE01


	//METALS
#ifndef   DETAILED_METALS_AND_MASS_RETURN
	o->MetalsColdGas     = CORRECTDBFLOAT(g->MetalsColdGas);
	o->MetalsStellarMass = CORRECTDBFLOAT(g->MetalsDiskMass) + CORRECTDBFLOAT(g->MetalsBulgeMass);
	o->MetalsDiskMass    = CORRECTDBFLOAT(g->MetalsDiskMass);
	o->MetalsBulgeMass   = CORRECTDBFLOAT(g->MetalsBulgeMass);
	o->MetalsHotGas      = CORRECTDBFLOAT(g->MetalsHotGas);
	o->MetalsEjectedMass = CORRECTDBFLOAT(g->MetalsEjectedMass);
#ifdef METALS_SELF
	o->MetalsHotGasSelf = CORRECTDBFLOAT(g->MetalsHotGasSelf);
#endif // METALS_SELF
#else // DETAILED_METALS_AND_MASS_RETURN
	o->MetalsColdGas = g->MetalsColdGas;
	o->MetalsStellarMass =  metals_add(g->MetalsDiskMass,g->MetalsBulgeMass,1.0);
	o->MetalsDiskMass = g->MetalsDiskMass;
	o->MetalsBulgeMass = g->MetalsBulgeMass;
	o->MetalsHotGas = g->MetalsHotGas;
	o->MetalsEjectedMass = g->MetalsEjectedMass;
#ifdef METALS_SELF
	o->MetalsHotGasSelf = g->MetalsHotGasSelf;
#endif // METALS_SELF
#endif // DETAILED_METALS_AND_MASS_RETURN

#ifdef TRACK_MASSGROWTH_CHANNELS
	o->MassFromInSitu=g->MassFromInSitu;
	o->MassFromMergers=g->MassFromMergers;
	o->MassFromBursts=g->MassFromBursts;
#endif // TRACK_MASSGROWTH_CHANNELS

#ifdef TRACK_BURST
	o->BurstMass=g->BurstMass;
#endif // TRACK_BURST

#ifdef OUTPUT_RINGS
	for(jj=0; jj<RNUM; jj++)
	{
		o->H2fractionRings[jj] = g -> H2fractionRings[jj];
		o->ColdGasRings[jj] = g->ColdGasRings[jj];
		o->MetalsColdGasRings[jj] = g->MetalsColdGasRings[jj];
		o->DiskMassRings[jj] = g->DiskMassRings[jj];
		o->MetalsDiskMassRings[jj] = g->MetalsDiskMassRings[jj];
#ifdef RINGS_IN_BULGES
		o->BulgeMassRings[jj] = g->BulgeMassRings[jj];
		o->MetalsBulgeMassRings[jj] = g->MetalsBulgeMassRings[jj];
#endif // RINGS_IN_BULGES
	}
#endif // OUTPUT_RINGS


	//STAR FORMATION HISTORIES / RATES

#ifdef OUTPUT_SFH
	o->sfh_ibin = g->sfh_ibin;
	for ( j = 0; j <= o->sfh_ibin; j ++ ) {
#ifndef NORMALIZEDDB
		// ROB: Lookback time to middle of SFH bin, in years. Now use LookBackTimeToSnap + sfh_time instead.
		//o->sfh_time[j]=(g->sfh_t[j]+g->sfh_dt[j]/2.)*UnitTime_in_years/Hubble_h;
		//Time from middle of this sfh bin to snapshot - converted from code units to years
		//o->sfh_time[j]=(g->sfh_t[j]+g->sfh_dt[j]/2.-NumToTime(g->SnapNum))*UnitTime_in_years/Hubble_h;
		//o->sfh_dt[j]=g->sfh_dt[j]*UnitTime_in_years/Hubble_h;
		o->sfh_DiskMass[j]  = g->sfh_DiskMass[j];
		o->sfh_BulgeMass[j] = g->sfh_BulgeMass[j];
#ifdef OUTPUT_RINGS
		for(jj=0; jj<RNUM; jj++)
		  {
			o->sfh_DiskMassRings[jj][j]=g->sfh_DiskMassRings[jj][j];
#ifdef RINGS_IN_BULGES
			o->sfh_BulgeMassRings[jj][j]=g->sfh_BulgeMassRings[jj][j];
#endif // RINGS_IN_BULGES
		  }
#endif // OUTPUT_RINGS
		o->sfh_ICM[j]             = g->sfh_ICM[j];
		o->sfh_MetalsDiskMass[j]  = g->sfh_MetalsDiskMass[j];
		o->sfh_MetalsBulgeMass[j] = g->sfh_MetalsBulgeMass[j];
		o->sfh_MetalsICM[j]       = g->sfh_MetalsICM[j];
		//#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef OUTPUT_ELEMENTS
		for(kk=0;kk<NUM_ELEMENTS;kk++)
		  {
			o->sfh_DiskMass_elements[j][kk]=g->sfh_DiskMass_elements[j][kk];
			o->sfh_BulgeMass_elements[j][kk]=g->sfh_BulgeMass_elements[j][kk];
			o->sfh_ICM_elements[j][kk]=g->sfh_ICM_elements[j][kk];
		  }
#endif // OUTPUT_ELEMENTS
		//#endif // DETAILED_METALS_AND_MASS_RETURN
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
		o->sfh_MassFromInSitu[j]=g->sfh_MassFromInSitu[j];
		o->sfh_MassFromMergers[j]=g->sfh_MassFromMergers[j];
		o->sfh_MassFromBursts[j]=g->sfh_MassFromBursts[j];
#endif // TRACK_SFH_MASSGROWTH_CHANNELS

#ifdef TRACK_BURST
		o->sfh_BurstMass[j]=g->sfh_BurstMass[j];
#endif // TRACK_BURST
#else // NORMALIZEDDB
		sfh_bin[j].sfh_DiskMass  = g->sfh_DiskMass[j];
		sfh_bin[j].sfh_BulgeMass = g->sfh_BulgeMass[j];
#ifdef OUTPUT_RINGS
		for(jj=0; jj<RNUM; jj++)
		  {
			sfh_bin[j].sfh_DiskMassRings[jj] = g->sfh_DiskMassRings[jj][j];
#ifdef RINGS_IN_BULGES
			sfh_bin[j].sfh_BulgeMassRings[jj] = g->sfh_BulgeMassRings[jj][j];
#endif // RINGS_IN_BULGES
		  }
#endif // OUTPUT_RINGS
		sfh_bin[j].sfh_ICM             = g->sfh_ICM[j];
		sfh_bin[j].sfh_MetalsDiskMass  = g->sfh_MetalsDiskMass[j];
		sfh_bin[j].sfh_MetalsBulgeMass = g->sfh_MetalsBulgeMass[j];
		sfh_bin[j].sfh_MetalsICM       = g->sfh_MetalsICM[j];
		sfh_bin[j].sfh_ibin            = j;
		sfh_bin[j].snapnum             = g->SnapNum;
		sfh_bin[j].GalID               = g->GalTreeIndex; // TODO must be reset
		o->GalID                       = g->GalTreeIndex;
#endif // NORMALIZEDDB
	}

	//Set all non-used array elements to zero: important if we want to read files in database that all values are valid
	// SQLServer floats
	for ( j = o->sfh_ibin + 1;
		  j < SFH_NBIN;
		  j ++ ) {
#ifndef NORMALIZEDDB
		//o->sfh_time[j]=0.;
		//o->sfh_dt[j]=0.;
		o->sfh_DiskMass[j]  = 0.;
		o->sfh_BulgeMass[j] = 0.;
#ifdef OUTPUT_RINGS
		for(jj=0; jj<RNUM; jj++)
		  {
			o->sfh_DiskMassRings[jj][j]=0.;
#ifdef RINGS_IN_BULGES
			o->sfh_BulgeMassRings[jj][j]=0.;
#endif // RINGS_IN_BULGES
		  }
#endif // OUTPUT_RINGS
		o->sfh_ICM[j]            = 0.;
		o->sfh_MetalsDiskMass[j] = metals_init();

		o->sfh_MetalsBulgeMass[j] = metals_init();

		o->sfh_MetalsICM[j] = metals_init();

#ifdef OUTPUT_ELEMENTS
		for(kk=0;kk<NUM_ELEMENTS;kk++)
		  {
			o->sfh_DiskMass_elements[j][kk]=0.;
			o->sfh_BulgeMass_elements[j][kk]=0.;
			o->sfh_ICM_elements[j][kk]=0.;
		  }
#endif // OUTPUT_ELEMENTS
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
		o->sfh_MassFromInSitu[j]=0.;
		o->sfh_MassFromMergers[j]=0.;
		o->sfh_MassFromBursts[j]=0.;
#endif // TRACK_SFH_MASSGROWTH_CHANNELS
#ifdef TRACK_BURST
		o->sfh_BurstMass[j]=0.;
#endif // TRACK_BURST
#else // NORMALIZEDDB
		sfh_bin[j].sfh_DiskMass=0;
		sfh_bin[j].sfh_BulgeMass=0;
#ifdef OUTPUT_RINGS
		for(jj=0; jj<RNUM; jj++)
		  {
			sfh_bin[j].sfh_DiskMassRings[jj]=0.;
#ifdef RINGS_IN_BULGES
			sfh_bin[j].sfh_BulgeMassRings[jj]=0.;
#endif // RINGS_IN_BULGES
		  }
#endif // OUTPUT_RINGS
		sfh_bin[j].sfh_ICM=0;
		sfh_bin[j].sfh_ibin = 0;
		sfh_bin[j].snapnum = g->SnapNum;
		// TODO other elements not important, are not being written anyway. Or are they used elsewhere?
#endif // NORMALIZEDDB
	}
#endif //OUTPUT_SFH

#ifdef OUTPUT_ELEMENTS
	for(kk=0;kk<NUM_ELEMENTS;kk++)
	  {
		o->DiskMass_elements[kk] = g->DiskMass_elements[kk];
		o->BulgeMass_elements[kk] = g->BulgeMass_elements[kk];
		o->ColdGas_elements[kk] = g->ColdGas_elements[kk];
		o->HotGas_elements[kk] = g->HotGas_elements[kk];
		o->EjectedMass_elements[kk] = g->EjectedMass_elements[kk];
		o->ICM_elements[kk] = g->ICM_elements[kk];
#ifdef OUTPUT_RINGS
		for(jj=0; jj<RNUM; jj++)
	  {
		o->DiskMassRings_elements[jj][kk] = g -> DiskMassRings_elements[jj][kk];
#ifdef RINGS_IN_BULGES
		o->BulgeMassRings_elements[jj][kk] = g -> BulgeMassRings_elements[jj][kk];
#endif // RINGS_IN_BULGES
		o->ColdGasRings_elements[jj][kk] = g->ColdGasRings_elements[jj][kk];
	  }
#endif // OUTPUT_RINGS
	  }
#endif //OUTPUT_ELEMENTS

	o->PrimordialAccretionRate = CORRECTDBFLOAT(
			                             g->PrimordialAccretionRate * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR /
			                             SOLAR_MASS);
	o->CoolingRate             = CORRECTDBFLOAT(
			                             g->CoolingRate * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
	o->CoolingRate_beforeAGN   = CORRECTDBFLOAT(
			                             g->CoolingRate_beforeAGN * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR /
			                             SOLAR_MASS);

	//o->MajorMergerRate = CORRECTDBFLOAT(g->MajorMergerRate / UnitTime_in_s * SEC_PER_YEAR);
	//o->MinorMergerRate = CORRECTDBFLOAT(g->MinorMergerRate / UnitTime_in_s * SEC_PER_YEAR);
	o->NMajorMergers = g->NMajorMergers;
	o->NMinorMergers = g->NMinorMergers;

	//NOTE: in Msun/yr
	o->Sfr      = CORRECTDBFLOAT(g->Sfr * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
	o->SfrBulge = CORRECTDBFLOAT(g->SfrBulge * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
#ifdef OUTPUT_RINGS
	for(jj=0; jj<RNUM; jj++)
	  o->SfrRings[jj] = g->SfrRings[jj] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
#endif // OUTPUT_RINGS
#endif //NO_PROPS_OUTPUTS


	//MAGNITUDES
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef POST_PROCESS_MAGS
	
	/*   int N_vespa_files=29, N_vespa_AgeBins=16;
		 double vespa_age[17]={0.000125893, 0.02000, 0.03000, 0.04800, 0.07400, 0.11500,
							0.17700, 0.27500, 0.42500, 0.65800, 1.02000, 1.57000,
							2.44000, 3.78000, 5.84000, 9.04000, 14.0000};
		 int vespa_sfh_IDs[29]={0,   1, 12, 16,  2, 22, 23, 27, 29, 36,
								42, 44, 50, 53, 54, 55, 56, 57, 61, 62,
								63, 64, 65, 71, 81, 82, 90, 94, 99};
		 double vespa_sfh[N_vespa_AgeBins], vespa_metal[N_vespa_AgeBins], dumb[6];
		 int ii, jj, kk;
		 char buf[1000], sbuf[1000];
		 FILE *fa;
	
		 for (ii=1;ii<N_vespa_files;ii++)
		  {
	
			 sprintf(buf, "./devel/vespa_sfh/disc_good_%d.txt", vespa_sfh_IDs[ii]);
			 if(!(fa = fopen(buf, "r")))
			 {
				 char sbuf[1000];
				 sprintf(sbuf, "can't open file `%s'\n", buf);
				 terminate(sbuf);
			 }
	
			 for(jj=0;jj<6;jj++)
				 fgets(buf, 300, fa);
	
			 //read vespa SFH and metallicities
			 for (jj=0;jj<N_vespa_AgeBins;jj++)
			 {
				 vespa_sfh[jj]=0.;
				 vespa_metal[jj]=0.;
				 fscanf(fa,"%lg %lg %lg %lg %lg %lg %lg %lg\n", &dumb[0], &dumb[1], &vespa_sfh[jj], &dumb[2], &vespa_metal[jj],
														  &dumb[3], &dumb[4], &dumb[5]);
			 }
			 fclose(fa);
	
			 for (jj=0;jj<SFH_NBIN;jj++)
			  {
				 if(jj<N_vespa_AgeBins)
				 {
					 o->sfh_time[jj]=(vespa_age[jj]+vespa_age[jj+1])/2.*1.e9;
					 o->sfh_dt[jj]=(vespa_age[jj+1]-vespa_age[jj])/2.*1.e9;
				   o->sfh_DiskMass[jj]=vespa_sfh[jj];
				   o->sfh_BulgeMass[jj]=0.;
				   o->sfh_MetalsDiskMass[jj]=vespa_metal[jj]*o->sfh_DiskMass[jj];
				   o->sfh_MetalsBulgeMass[jj]=metals_init();
				 }
				 else
				 {
					 o->sfh_time[jj]=0.;
					 o->sfh_dt[jj]=0.;
					 o->sfh_DiskMass[jj]=0.;
					 o->sfh_BulgeMass[jj]=0.;
					 o->sfh_MetalsDiskMass[jj]=metals_init();
					 o->sfh_MetalsBulgeMass[jj]=metals_init();
				 }
			  }
	
			  post_process_spec_mags(o);
	
			  sprintf(buf, "./devel/vespa_sfh/output_spectradisc_good_%d.txt", vespa_sfh_IDs[ii]);
			  if(!(fa = fopen(buf, "w")))
			  {
				  char sbuf[1000];
				  sprintf(sbuf, "can't open file `%s'\n", buf);
				  terminate(sbuf);
			  }
	
			  for(jj=0;jj<NMAG;jj++)
				  fprintf(fa,"%e\n",o->Mag[jj]);
			  fclose(fa);
	
		   exit(0);
		 }
	*/
		// Convert recorded star formation histories into mags
#ifdef NORMALIZEDDB
	  post_process_spec_mags(o, &(sfh_bin[0]));
#else // NORMALIZEDDB
	  post_process_spec_mags(o);
#endif // NORMALIZEDDB
#else // POST_PROCESS_MAGS

#ifdef OUTPUT_REST_MAGS
	// Luminosities are converted into Mags in various bands
	for(j = 0; j < NMAG; j++)
	  {
		//o->Mag[j] = lum_to_mag(g->Lum[j][n]); -> DONE ON TOP FOR LIGHT_OUTPUT AS WELL
		o->MagBulge[j] = lum_to_mag(g->LumBulge[j][n]);
		o->MagDust[j] = lum_to_mag(g->LumDust[j][n]);
#ifdef ICL
		o->MagICL[j] = lum_to_mag(g->ICLLum[j][n]);
#endif // ICL
	  }


#endif //OUTPUT_REST_MAGS
#ifdef OUTPUT_OBS_MAGS
#ifdef COMPUTE_OBS_MAGS
	// Luminosities in various bands
	for(j = 0; j < NMAG; j++)
	  {
		o->ObsMag[j] = lum_to_mag(g->ObsLum[j][n]);
		o->ObsMagBulge[j] = lum_to_mag(g->ObsLumBulge[j][n]);
		o->ObsMagDust[j] = lum_to_mag(g->ObsLumDust[j][n]);
#ifdef ICL
		o->ObsMagICL[j] = lum_to_mag(g->ObsICL[j][n]);
#endif // ICL

#ifdef OUTPUT_MOMAF_INPUTS
		o->dObsMag[j] = lum_to_mag(g->dObsLum[j][n]);
		o->dObsMagBulge[j] = lum_to_mag(g->dObsLumBulge[j][n]);
		o->dObsMagDust[j] = lum_to_mag(g->dObsLumDust[j][n]);
#ifdef ICL
		o->dObsMagICL[j] = lum_to_mag(g->dObsICL[j][n]);
#endif // ICL
#endif // OUTPUT_MOMAF_INPUTS
	  }
#endif // COMPUTE_OBS_MAGS
#endif // OUTPUT_OBS_MAGS
#endif // POST_PROCESS_MAGS
#endif // COMPUTE_SPECPHOT_PROPERTIES

#ifndef NO_PROPS_OUTPUTS
	if ((g->DiskMass + g->BulgeMass) > 0.0) {
		o->MassWeightAge            = g->MassWeightAge[n] / (g->DiskMass + g->BulgeMass);
		o->MassWeightAge            = o->MassWeightAge / 1000. * UnitTime_in_Megayears / Hubble_h; // Age in Gyr
		o->StellarHalfLightRadius   = stellar_half_light_radius (o);
		o->StellarNinetyLightRadius = stellar_ninety_light_radius (o);
	} else {
		o->MassWeightAge            = 0.;
		o->StellarHalfLightRadius   = 0.;
		o->StellarNinetyLightRadius = 0.;
	}
#endif // NO_PROPS_OUTPUTS

#ifdef FIX_OUTPUT_UNITS
	fix_units_for_ouput(o);
#endif // FIX_OUTPUT_UNITS
#endif // LIGHT_OUTPUT
}
// end prepare_galaxy_for_output



















/**@brief Removes h from units of galaxy properties. If desired (makefile option FIX_OUTPUT_UNITS is set), the output
 * properties of the galaxies can be scaled to physical units excluding any factors of h == Hubble/100 km/s/Mpc. 	  */
#ifdef FIX_OUTPUT_UNITS
#ifndef NO_PROPS_OUTPUTS
#ifdef LIGHT_OUTPUT
void fix_units_for_ouput(struct GALAXY_OUTPUT *o)
{
  o->Pos[0] /= Hubble_h;
  o->Pos[1] /= Hubble_h;
  o->Pos[2] /= Hubble_h;
  o->Mvir /= Hubble_h;
  o->Rvir /= Hubble_h;
  o->ColdGas /= Hubble_h;
  o->DiskMass /= Hubble_h;
  o->BulgeMass /= Hubble_h;
  o->PBMass /= Hubble_h;
  o->CBMassMajor /= Hubble_h;
  o->CBMassMinor /= Hubble_h;
  o->CBMass /= Hubble_h;
  o->HotGas /= Hubble_h;
  o->BlackHoleMass /= Hubble_h;

}
// end fix_units_for_ouput
#else // LIGHT_OUTPUT

/**@brief Removes h from units of galaxy properties. If desired (makefile option FIX_OUTPUT_UNITS is set), the output
properties of the galaxies can be scaled to physical units excluding any factors of h == Hubble/100 km/s/Mpc.         */
void fix_units_for_ouput(struct GALAXY_OUTPUT *o)
{
  int j;

  o->Pos[0] /= Hubble_h;
  o->Pos[1] /= Hubble_h;
  o->Pos[2] /= Hubble_h;
  o->CentralMvir /= Hubble_h;
  o->Mvir /= Hubble_h;
  o->Rvir /= Hubble_h;
  
#ifdef HALOPROPERTIES
  o->HaloM_Mean200 /= Hubble_h;
  o->HaloM_Crit200 /= Hubble_h;
  o->HaloM_TopHat /= Hubble_h;
  o->HaloPos[0] /= Hubble_h;
  o->HaloPos[1] /= Hubble_h;
  o->HaloPos[2] /= Hubble_h;
  o->HaloSpin[0] /= Hubble_h;
  o->HaloSpin[1] /= Hubble_h;
  o->HaloSpin[2] /= Hubble_h;
#endif // HALOPROPERTIES

  o->ColdGasSpin[0] /= Hubble_h;
  o->ColdGasSpin[1] /= Hubble_h;
  o->ColdGasSpin[2] /= Hubble_h;
  o->DiskSpin[0] /= Hubble_h;
  o->DiskSpin[1] /= Hubble_h;
  o->DiskSpin[2] /= Hubble_h;
  o->BulgeSpin[0] /= Hubble_h;
  o->BulgeSpin[1] /= Hubble_h;
  o->BulgeSpin[2] /= Hubble_h;
  o->HotRadius /= Hubble_h;
  o->ColdGas /= Hubble_h;
  o->DiskMass /= Hubble_h;
  o->BulgeMass /= Hubble_h;
  o->PBMass /= Hubble_h;
  o->CBMassMajor /= Hubble_h;
  o->CBMassMinor /= Hubble_h;
  o->CBMass /= Hubble_h;
  o->HotGas /= Hubble_h;
  o->EjectedMass /= Hubble_h;
  o->BlackHoleMass /= Hubble_h;
  o->ICM /= Hubble_h;
  o->BulgeSize /= Hubble_h;
  o->DiskRadius /= Hubble_h;
  o->ColdGasRadius /= Hubble_h;
  o->CoolingRadius /= sqrt(Hubble_h);

#ifdef TRACK_MASSGROWTH_CHANNELS
  o->MassFromInSitu/= Hubble_h;
  o->MassFromMergers/= Hubble_h;
  o->MassFromBursts/= Hubble_h;
#endif // TRACK_MASSGROWTH_CHANNELS

#ifdef TRACK_BURST
  o->BurstMass /= Hubble_h;
#endif // TRACK_BURST

  o->MetalsColdGas=metals_add(metals_init(),o->MetalsColdGas,1./Hubble_h);
  o->MetalsDiskMass=metals_add(metals_init(),o->MetalsDiskMass,1./Hubble_h);
  o->MetalsBulgeMass=metals_add(metals_init(),o->MetalsBulgeMass,1./Hubble_h);
  o->MetalsHotGas=metals_add(metals_init(),o->MetalsHotGas,1./Hubble_h);
  o->MetalsEjectedMass=metals_add(metals_init(),o->MetalsEjectedMass,1./Hubble_h);
  o->MetalsICM=metals_add(metals_init(),o->MetalsICM,1./Hubble_h);

#ifdef OUTPUT_SFH
  for(j=0;j<=o->sfh_ibin;j++) {
	o->sfh_DiskMass[j] /= Hubble_h;
	o->sfh_BulgeMass[j] /= Hubble_h;
	o->sfh_ICM[j] /= Hubble_h;
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
	o->sfh_MassFromInSitu[j]/= Hubble_h;
	o->sfh_MassFromMergers[j]/= Hubble_h;
	o->sfh_MassFromBursts[j]/= Hubble_h;
#endif // TRACK_SFH_MASSGROWTH_CHANNELS
#ifdef TRACK_BURST
	o->sfh_BurstMass[j] /= Hubble_h;
#endif // TRACK_BURST
	o->sfh_MetalsDiskMass[j]=metals_add(metals_init(),o->sfh_MetalsDiskMass[j],1./Hubble_h);
	o->sfh_MetalsBulgeMass[j]=metals_add(metals_init(),o->sfh_MetalsBulgeMass[j],1./Hubble_h);
	o->sfh_MetalsICM[j]=metals_add(metals_init(),o->sfh_MetalsICM[j],1./Hubble_h);
  }
#endif //OUTPUT_SFH
}
#endif // LIGHT_OUTPUT
#endif // NO_PROPS_OUTPUTS
#endif // FIX_OUTPUT_UNITS
// end fix_units_for_ouput



















/** @brief																											  */
long long calc_big_db_offset (int filenr, int treenr) {
	long long big;
#ifdef MRII
	// TODO assumes < 10^9 galaxies per tree, is that always correct?
	big = (((filenr * (long long) 1000000) + treenr) * (long long) 1000000000);
#else // MRII
	big = (((filenr * (long long) 1000000) + treenr) * (long long) 1000000);
#endif // MRII
	return big;
}
// end calc_big_db_offset



















/** @brief																											  */
long long calc_big_db_subid_index (int snapnum, int filenr, int subhaloindex) {
	long long big;
#ifdef MRII
	big = snapnum * (long long) 10000000000000 + filenr * (long long) 1000000000 + subhaloindex;
#else // MRII
	big = snapnum * (long long) 1000000000000 + filenr * (long long) 100000000 + subhaloindex;
#endif // MRII
	return big;
}
// end calc_big_db_subid_index



















/** @brief																											  */
#ifdef OUTPUT_SFH
void write_sfh_bins()
{
	FILE* SFH_Bins_File;
	struct SFH_Time *sfh_times;
	char buf[1000];

  int snap,j, nbins,ibin;
  nbins = 0;
  for(snap = 0; snap < MAXSNAPS; snap++)
  {
	  for(j=0;j < SFH_ibin[snap][0];j++)
	  {
		  nbins++;
	  }
  }

  sfh_times = (struct SFH_Time *) mymalloc("sfh_times", sizeof(struct SFH_Time) * nbins);
  ibin = 0;
  for(snap = 0; snap < MAXSNAPS; snap++)
  {
	  for(j=0;j < SFH_ibin[snap][0];j++)
	  {
		  sfh_times[ibin].snapnum = snap;
		  sfh_times[ibin].bin = j;
		  sfh_times[ibin].lookbacktime = (SFH_t[snap][0][j]+SFH_dt[snap][0][j]/2.-NumToTime(snap))*UnitTime_in_years/Hubble_h;
		  sfh_times[ibin].dt=SFH_dt[snap][0][j]*UnitTime_in_years/Hubble_h;
		  sfh_times[ibin].nbins = SFH_Nbins[snap][0][j];
		  ibin++;
	  }
  }
	sprintf(buf, "%s/SFH_Bins", FinalOutputDir);
	if(!(SFH_Bins_File = fopen(buf, "w")))
	{
	  char sbuf[1000];
	  sprintf(sbuf, "can't open file `%s'\n", buf);
	  terminate(sbuf);
	} else
	  printf("writing sfh bins to %s\n", buf);

	// write # bins
	myfwrite(&nbins, sizeof(int), 1, SFH_Bins_File);	// write 1
	// skip header
	myfwrite(sfh_times, sizeof(struct SFH_Time), nbins, SFH_Bins_File);	// size of an output structure (Galaxy_Output)
	fflush(SFH_Bins_File);
	fclose(SFH_Bins_File);
}
#endif // OUTPUT_SFH
// end OUTPUT_SFH