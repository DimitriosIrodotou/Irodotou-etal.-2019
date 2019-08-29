#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

/**@brief Calculates the merging time, the central galaxy (for type 1's), adds galaxies together, calculates SF from
 * bursts and grows black holes.
 *
 * <B>estimate_merging_time</B> sets up a merger clock. Originally this was done only for type 2 galaxies. The
 * positions of galaxies in the model are normally given by the position of their dark matter halo. However, when
 * galaxies become satellites, their dark matter haloes are stripped by the central object to the point where there is
 * none left. At this point, the dark matter of the satellite becomes part of the main halo, but the galaxy's position
 * should continue to evolve due to the dynamical friction force caused by the dark matter aroundit.
 *
 * The code does keep track of the position of the most bounded particle when the satellite's halo was disrupted, but
 * this is not used to track the galaxy position. Instead a clock is set, giving the time left until the satellite
 * mergers with the central galaxy. Before, this was only done for type 2's (satellites that lost a halo). Guo2010
 * included the MERGE01 option, that sets the clock also for type 1's (satellites with a halo), since for the highest
 * resolution millennium 2, their haloes can be very small and orbit forever around the central companion.
 *
 * This time is computed using the Chandrasekhar's formula for dynamical friction, as in Binney & Tremaine 1987:
 *       \f$F_{\rm{df}}=
 *               -\frac{4\pi {\rm{G}}^2 m^2_{\rm{sat}} \ln(\Lambda) \rho B(x)}
 *                     {v^2_{\rm{rel}}}.\f$
 *
 *       Which gives (B&T Eq. 7.26):
 *
 *       \f$t_{df}\approx1.17\frac{V_{\rm{vir}}r^2_{\rm{sat}}}
 *                                {{\rm{G}}\,m_{\rm{sat}}\ln(\Lambda)},\f$
 *
 * that is afterwards multiplied by 2 (after Delucia2007 to fit the data). When the merging time reaches zero, the
 * satellite is assumed to merge with the central galaxy.
 *
 * <B>deal_with_galaxy_merger</B> deals with the process, according to the mass fraction of the merger. Major if
 * \f$M_{\rm{sat}}/M_{\rm{central}}>0.3\f$ and minor otherwise. It calls
 *
 * - add_galaxies_together - Add the cold and stellar phase of the merged galaxy to the central one. Also form a bulge
 * at the central galaxy with the stars from the satellite in a minor merger if BulgeFormationInMinorMergersOn=1 (Major
 * mergers are dealt later).
 * - Then calls grow_black_hole - Grows black hole through accretion from cold gas during mergers (due to the
 * instabilities triggered), as in Kauffmann & Haehnelt (2000). This is commonly referred as the quasar mode, main
 * responsible for the black hole growth. After Croton2006 this mode is active even in minor mergers:
 *       \f$\Delta m_{\rm BH,Q}=M_{\rm{BH,min}}
 *          \frac{f_{\rm BH}(m_{\rm sat}/m_{\rm central})\,m_{\rm cold}}
 *           {1+(280\,\mathrm{km\,s}^{-1}/V_{\rm vir})^2}.\f$
 
 * - Finally the burst of star formation due to the merger is treated.
 *
 * - If StarBurstModel = 0 (since Croton2006), the Somerville 2001 model of bursts is used
 * collisional_starburst_recipe(). The burst can happen for both major and minor mergers, with a fraction of the added
 * cold gas from the satellite and central being consumed:
 *
 * \f$\dot{m}_{\star}^{\rm{burst}}
 *             = 0.56 \left(\frac{m_{\rm{sat}}}{m_{\rm{central}}}\right)^{0.7}
 *               m_{\rm{gas}}\f$.
 *           SN Feedback from starformation is computed and the sizes of bulge
 *           and disk followed.
 *
 * - When a major merger occurs, the disk of both merging galaxies is completely destroyed to form a bulge. In either
 * type of mergers, the bulge size is updated using Eq. 33 in Guo2010:
 *       \f$C\frac{GM^2_{\rm{new,bulge}}}{R_{\rm{new,bulge}}}=
 *          C\frac{GM^2_1}{R_1}+C\frac{GM^2_2}{R_2}+\alpha_{\rm{inter}}
 *          \frac{GM_1M_2}{R_1+R_2}\f$                                                                                */



















/** @brief Calculates the merging time whenever a galaxy becomes a satellite. Binney & Tremaine 1987 - 7.26 merging
 * time for satellites due to dynamical friction. After Delucia2007 *2, shown to agree with Kolchin2008 simulations in
 * Delucia2010. This is set when a galaxy becomes a type 2 or being a type 1 \f$M_{\rm{star}}>M_{\rm{vir}}\f$. In
 * DeLucia2007 they could only merge into a type 0, now (after guo2010) they can merge into a type 1.                   */
double estimate_merging_time (int halonr, int mother_halonr, int p)
{
	int    central_halonr;
	double coulomb, mergtime, SatelliteMass, SatelliteRadius, MotherHaloRvir;

	// Recipe updated for more accurate merging time (see BT eq 7.26), now satellite radius at previous timestep is
	// included.
	central_halonr = Halo[Halo[halonr].Descendant].FirstProgenitor;
	if (Gal[p].Type == 1)
	{
		central_halonr = mother_halonr;
	}
	if (central_halonr == halonr)
	{
		terminate("can't be...!\n");
	}

	coulomb = log (Halo[mother_halonr].Len / ((double) Halo[halonr].Len) + 1);

	//  Should include stellar+cold gas in SatelliteMass!
	SatelliteMass   = get_virial_mass (halonr) + (Gal[p].DiskMass + Gal[p].BulgeMass);
	SatelliteRadius = separation_halo (central_halonr, halonr) / (1 + ZZ[Halo[halonr].SnapNum]);

	int j;
	for (j = 0; j < 3; j ++)
	{
		Gal[p].DistanceToCentralGal[j] = wrap(Halo[central_halonr].Pos[j] - Halo[halonr].Pos[j], BoxSize);
	}

	MotherHaloRvir = get_virial_radius (mother_halonr);
	if (SatelliteRadius > MotherHaloRvir)
	{
		SatelliteRadius = MotherHaloRvir;
	}

	if (SatelliteMass > 0.0)
	{
		mergtime = 1.17 * SatelliteRadius * SatelliteRadius * get_virial_velocity (mother_halonr) /
		           (coulomb * G * SatelliteMass); // Binney & Tremaine Eq.7.26
		// Change introduced by Delucia2007 to fit observations
		mergtime = MergerTimeMultiplier * mergtime;
	}
	else
	{
		mergtime = - 99999.9;
	}
	return mergtime;
}
// end estimate_merging_time



















/** @brief Deals with the physics triggered by mergers, according to the mass fraction of the merger
 * \f$(M_{\rm{sat}}/M_{\rm{central}}><0.3)\f$. Add the cold and stellar phase of the satellite galaxy to the central
 * one, form a bulge at the central galaxy with the stars from the satellite in a minor merger if
 * BulgeFormationInMinorMergersOn=1. Grows black hole through accretion from cold gas "quasar mode". If StarBurstModel
 * = 0, the Somerville 2001 model of bursts is used, SN Feedback from starformation is computed and the sizes of bulge
 * and disk followed. When a major merger occurs, the disk of both merging galaxies is completely destroyed to form
 * a bulge. New stars form of to the bulge                                                                              */
void deal_with_galaxy_merger (int p, double time, double deltaT, int nstep)
{
	int    FOF_centralgal, merger_centralgal;
	double mi, ma, mass_ratio, MgasCentral, MstarCentral, MbulgeCentral, MgasSat, MstarSat, MbulgeSat, frac;
	double RgasCentral, RStellarDiskCentral, RgasSat, RStellarDiskSat, HMRCentral, HMRSat;

	FOF_centralgal    = Gal[p].FOFCentralGal;
	merger_centralgal = Gal[p].MergerCentralGal;

	mass_checks (p, "model_mergers.c", __LINE__);
	mass_checks (merger_centralgal, "model_mergers.c", __LINE__);
	mass_checks (FOF_centralgal, "model_mergers.c", __LINE__);

#ifdef GALAXYTREE
	int q;

	q = Gal[merger_centralgal].FirstProgGal;
	if ( q >= 0 ) {
		while (GalTree[q].NextProgGal >= 0)
			q = GalTree[q].NextProgGal;
		GalTree[q].NextProgGal = Gal[p].FirstProgGal;
		if ( GalTree[q].NextProgGal >= NGalTree ) {
			printf("q=%d p=%d GalTree[q].NextProgGal=%d NGalTree=%d\n", q, p, GalTree[q].NextProgGal, NGalTree);
			terminate("problem");
		}
	}
	if ( q < 0 ) {
		terminate("q<0\n");
	}
	q = GalTree[q].NextProgGal;

	if ( HaloGal[GalTree[q].HaloGalIndex].GalTreeIndex != q ) {
		terminate("inconsistency");
	}

	HaloGal[GalTree[q].HaloGalIndex].MergeOn = 2;

	if ( Gal[p].Type == 1 ) {
		HaloGal[GalTree[q].HaloGalIndex].MergeOn = 3;
	}
#endif // GALAXYTREE

	// Flag galaxy as finished.
	Gal[p].Type = 3;

	//  Calculate mass ratio of merging galaxies.
	mi = Gal[p].DiskMass + Gal[p].BulgeMass + Gal[p].ColdGas;
	ma = Gal[merger_centralgal].DiskMass + Gal[merger_centralgal].BulgeMass + Gal[merger_centralgal].ColdGas;

	if (max(mi, ma) > 0.)
	{
		mass_ratio = min(mi, ma) / max(mi, ma);
	}
	else
	{
		mass_ratio = 1.0;
	}

	// Record the before merger gas and stellar component of merger central and satellite galaxies.
	MstarCentral  = (Gal[merger_centralgal].DiskMass + Gal[merger_centralgal].BulgeMass);
	MbulgeCentral = Gal[merger_centralgal].BulgeMass;
	MgasCentral   = Gal[merger_centralgal].ColdGas;

	MstarSat  = (Gal[p].DiskMass + Gal[p].BulgeMass);
	MbulgeSat = Gal[p].BulgeMass;
	MgasSat   = Gal[p].ColdGas;

	mass_checks (p, "model_mergers.c", __LINE__);
	mass_checks (merger_centralgal, "model_mergers.c", __LINE__);

	// Compute the sizes before galaxies are merged to use in the calculation of the resulting bulge size (the
	// calculation is only needed here because with H2_AND_RINGS the sizes of stellar and gaseous disks are computed
	// on the fly)
	RgasCentral         = get_gas_disk_radius (merger_centralgal) / 3.;
	RStellarDiskCentral = get_stellar_disk_radius (merger_centralgal) / 3.;
	HMRCentral          = Gal[merger_centralgal].StellarHalfMassRadius;

	RgasSat         = get_gas_disk_radius (p) / 3.;
	RStellarDiskSat = get_stellar_disk_radius (p) / 3.;
	HMRSat          = Gal[p].StellarHalfMassRadius;

	// Add the cold and stellar phase of the merged galaxy to the central one. Also form a bulge if
	// BulgeFormationInMinorMergersOn is set on (transfer stars from satellite disk to central bulge). In a major
	// merger (dealt at the make_bulge_from_burst) the disk of the central (now made up of central and satellite will
	// be moved to the bulge). Any new stars formed will go to the bulge.

	add_galaxies_together (merger_centralgal, p, deltaT);

	mass_checks (p, "model_mergers.c", __LINE__);
	mass_checks (merger_centralgal, "model_mergers.c", __LINE__);

	// Grow black hole through accretion from cold disk during mergers, as in Kauffmann & Haehnelt (2000) + minor
	// mergers - Quasar Mode
	if (AGNRadioModeModel < 4)
	{
		grow_black_hole (merger_centralgal, mass_ratio, deltaT);
	}

	mass_checks (p, "model_mergers.c", __LINE__);
	mass_checks (merger_centralgal, "model_mergers.c", __LINE__);

	// Starburst as in Somerville 2001, with feedback computed inside. All star formation happens in the disk, but in a
	// major merger this will then be destroyed with everything moved to the bulge.
	if (StarBurstModel == 0)
	{
		frac = collisional_starburst_recipe (mass_ratio, merger_centralgal, FOF_centralgal, time, deltaT);
#ifdef DI_BULGESIZE
		size_of_merger_remnants(mass_ratio, merger_centralgal,  MstarCentral, MbulgeCentral,
								MgasCentral, MstarSat,  MgasSat, frac,HMRCentral, HMRSat);
#else // DI_BULGESIZE
		bulgesize_from_merger (mass_ratio, merger_centralgal, p, MstarCentral, MbulgeCentral, MgasCentral,
		                       MstarSat, MbulgeSat, MgasSat, frac, RgasCentral, RStellarDiskCentral, RgasSat,
		                       RStellarDiskSat);
#endif // DI_BULGESIZE

		mass_checks (p, "model_mergers.c", __LINE__);
		mass_checks (merger_centralgal, "model_mergers.c", __LINE__);
		mass_checks (FOF_centralgal, "model_mergers.c", __LINE__);

		if (mass_ratio > ThreshMajorMerger)
		{
			make_bulge_from_burst (merger_centralgal);
		}
	}

#ifdef RINGS_IN_BULGES
	//Bulge Mass was added into the same place as the disk, it will now be redistributed
	//according to a Jaffe profile and after the new bulge size has been calculated

	double rb = Gal[merger_centralgal].BulgeSize, TotMassInsideRings = 0., fractionRings[RNUM];
	int    jj;

	if ( rb > 0. )
		TotMassInsideRings = (RingRadius[RNUM - 1] / rb) / (1 + RingRadius[RNUM - 1] / rb);

	if ( TotMassInsideRings > 0. ) {
		fractionRings[0] = (RingRadius[0] / rb) / (1 + RingRadius[0] / rb) / TotMassInsideRings;
		for ( jj = 1; jj < RNUM; jj++ )
			fractionRings[jj] = ((RingRadius[jj] / rb) / (1 + RingRadius[jj] / rb) -
								 (RingRadius[jj - 1] / rb) / (1 + RingRadius[jj - 1] / rb)) / TotMassInsideRings;
	}
	else
		for ( jj = 0; jj < RNUM; jj++ )
			fractionRings[jj] = 0.;

	for ( jj = 0; jj < RNUM; jj++ ) {
		Gal[merger_centralgal].BulgeMassRings[jj] = fractionRings[jj] * Gal[merger_centralgal].BulgeMass;
		Gal[merger_centralgal].MetalsBulgeMassRings[jj] = metals_add(metals_init(),
																	 Gal[merger_centralgal].MetalsBulgeMass,
																	 fractionRings[jj]);
#ifdef INDIVIDUAL_ELEMENTS
		int kk;
		for(kk=0;kk<NUM_ELEMENTS;kk++)
			Gal[merger_centralgal].BulgeMassRings_elements[jj][kk] = fractionRings[jj]*Gal[merger_centralgal].BulgeMass_elements[kk];
#endif // INDIVIDUAL_ELEMENTS
#ifdef STAR_FORMATION_HISTORY
		int ii;
		for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
			Gal[merger_centralgal].sfh_BulgeMassRings[jj][ii] = fractionRings[jj]*Gal[merger_centralgal].sfh_BulgeMass[ii];
#endif // STAR_FORMATION_HISTORY
	}
#endif //RINGS_IN_BULGES

	mass_checks (p, "model_mergers.c", __LINE__);
	mass_checks (merger_centralgal, "model_mergers.c", __LINE__);
	mass_checks (FOF_centralgal, "model_mergers.c", __LINE__);

	// If we are in the presence of a minor merger, check disk stability (the disk is completely destroyed in major
	// mergers).
	if (DiskInstabilityModel == 0)
	{
		if (mass_ratio < ThreshMajorMerger &&
		    (Gal[merger_centralgal].DiskMass + Gal[merger_centralgal].BulgeMass) > 0.0)
		{
			check_disk_instability (merger_centralgal, deltaT / STEPS, time);
		}
	}
	mass_checks (merger_centralgal, "model_mergers.c", __LINE__);

	if ((Gal[merger_centralgal].BulgeMass > 1.e-6 && Gal[merger_centralgal].BulgeSize == 0.0) ||
	    (Gal[merger_centralgal].BulgeMass == 0.0 && Gal[merger_centralgal].BulgeSize > 1.e-6))
	{
		char sbuf[1000];
		sprintf (sbuf, "2 central: stellarmass %f, bulgemass %f, bulgesize %f, stellardisksize %f \n",
		         (Gal[merger_centralgal].DiskMass + Gal[merger_centralgal].BulgeMass),
		         Gal[merger_centralgal].BulgeMass, Gal[merger_centralgal].BulgeSize, Gal[merger_centralgal].DiskRadius);
		terminate(sbuf);
	}

	if (DiskRadiusModel == 0)
	{
		Gal[merger_centralgal].ColdGasRadius = get_gas_disk_radius (merger_centralgal);
		Gal[merger_centralgal].DiskRadius    = get_stellar_disk_radius (merger_centralgal);

		// The following may not be necessary - depends what happens in collisonal_starburst_recipe.
		if (merger_centralgal != FOF_centralgal)
		{
			Gal[FOF_centralgal].ColdGasRadius = get_gas_disk_radius (FOF_centralgal);
			Gal[FOF_centralgal].DiskRadius    = get_stellar_disk_radius (FOF_centralgal);
		}
	}

	mass_checks (p, "model_mergers.c", __LINE__);
	mass_checks (merger_centralgal, "model_mergers.c", __LINE__);

}
// end deal_with_galaxy_merger



















/** @brief Grows black hole through accretion from cold gas during mergers, as in Kauffmann & Haehnelt (2000). In
 * addition, black holes can grow during minor mergers. BlackHoleGrowth == 0 gives instantaneous accretion onto the
 * black hole; BlackHoleGrowth == 1 instead feeds an accretion disk: accretion occurs in main.c                       */
void grow_black_hole (int merger_centralgal, double mass_ratio, double deltaT)
{
	double BHaccrete, fraction;
#ifdef H2_AND_RINGS
	int    jj;
	double fractionRings[RNUM];
#endif

	if (Gal[merger_centralgal].ColdGas > 0.0)
	{
		BHaccrete = BlackHoleGrowthRate * mass_ratio * Gal[merger_centralgal].ColdGas
		            / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[merger_centralgal].Vvir)));
		// Redshift dependent accretion, not published
		// BHaccrete = BlackHoleGrowthRate * (1.0 + ZZ[Halo[halonr].SnapNum]) * mass_ratio

		// Cannot accrete more gas than is available
		if (BHaccrete > Gal[merger_centralgal].ColdGas)
		{
			BHaccrete = Gal[merger_centralgal].ColdGas;
		}
		fraction  = BHaccrete / Gal[merger_centralgal].ColdGas;

#ifdef H2_AND_RINGS
		for ( jj = 0; jj < RNUM; jj++ )
			fractionRings[jj] = fraction;
		// Doesn't matter the destination, since there are no rings in BlackHoleMass or BlackHoleGas but maybe there
		// will be in the future.
		/*        double mass_to_transfer;
		 mass_to_transfer = BHaccrete;
		 for ( jj = 0; jj < RNUM; jj++ )
		 fractionRings[jj] = 0.;
		 //transfer from inner circles first
		 for ( jj = 0; jj < RNUM; jj++ ) {
		 RStellarDiskSat
		 if ( mass_to_transfer < Gal[merger_centralgal].ColdGasRings[jj] ) {
		 fractionRings[jj] = mass_to_transfer / Gal[merger_centralgal].ColdGasRings[jj];
		 break;
		 }
		 else {
		 fractionRings[jj] = 1.;
		 mass_to_transfer -= Gal[merger_centralgal].ColdGasRings[jj];
		 }
		 }*/
		if ( BlackHoleGrowth == 0 ) {
			transfer_material_with_rings(merger_centralgal, "BlackHoleMass", merger_centralgal, "ColdGas",
										 fractionRings, "model_mergers.c", __LINE__);
		}
		else if ( BlackHoleGrowth == 1 ) {
			transfer_material_with_rings(merger_centralgal, "BlackHoleGas", merger_centralgal, "ColdGas", fractionRings,
										 "model_mergers.c", __LINE__);
		}
#else // H2_AND_RINGS
		if (BlackHoleGrowth == 0)
		{
			transfer_material (merger_centralgal, "BlackHoleMass", merger_centralgal, "ColdGas", fraction,
			                   "model_mergers.c", __LINE__);
			Gal[merger_centralgal].QuasarAccretionRate += BHaccrete / deltaT;
		}
		else if (BlackHoleGrowth == 1)
		{
			transfer_material (merger_centralgal, "BlackHoleGas", merger_centralgal, "ColdGas", fraction,
			                   "model_mergers.c", __LINE__);
		}
#endif //H2_AND_RINGS

		Gal[merger_centralgal].QuasarAccretionRate += BHaccrete / deltaT;
	}
}
// end grow_black_hole




















/** @brief All the components of the satellite galaxy are added to the correspondent component of the central galaxy.
 * Cold gas spin is updated and a bulge is formed at the central galaxy, with the stars of the satellite if
 * BulgeFormationInMinorMergersOn=1. In case of a major merger, everything that was put in the disk of the central
 * galaxy will be moved into the bulge                                                                                  */
void add_galaxies_together (int t, int p, double deltaT)
{

	int    outputbin, j, ii;
	float  tspin[3], tmass, pmass, ptrans;
	double mass_ratio, tmass_total, pmass_total;
	double tbulge, tdisk, pbulge, pdisk; // t central, p satellite

	mass_checks (p, "model_mergers.c", __LINE__);
	mass_checks (t, "model_mergers.c", __LINE__);

	tmass           = Gal[t].ColdGas;
	pmass           = Gal[p].ColdGas;
	tbulge          = Gal[t].BulgeMass;
	tdisk           = Gal[t].DiskMass;
	pbulge          = Gal[p].BulgeMass;
	pdisk           = Gal[p].DiskMass;

	Gal[t].MergeSat += (Gal[p].DiskMass + Gal[p].BulgeMass);
	Gal[p].MergeSat = 0.;

	tmass_total          = Gal[t].DiskMass + Gal[t].BulgeMass + Gal[t].ColdGas;
	pmass_total          = Gal[p].DiskMass + Gal[p].BulgeMass + Gal[p].ColdGas;

	if (max(tmass_total, pmass_total) > 0.)
	{
		mass_ratio = min(tmass_total, pmass_total) / max(tmass_total, pmass_total);
	}
	else
	{
		mass_ratio = 1.0;
	}

	// Count the number of each type of merger.
	Gal[t].NMajorMergers += Gal[p].NMajorMergers;
	Gal[t].NMinorMergers += Gal[p].NMinorMergers;

	if (mass_ratio > ThreshMajorMerger)
	{
		Gal[t].NMajorMergers += 1.;  //
	}
	else
	{
		Gal[t].NMinorMergers += 1.; //
	}

#ifdef H2_AND_RINGS
	double ringtot, fractionRings[RNUM], rd;
	int    jj;

	// Distribute the satellite gas evenly through the rings and then assume the same profile as for infall when it
	// goes to the disk of the central galaxy
	rd      = get_initial_disk_radius(Gal[t].HaloNr, t) / 3.;
	ringtot = 1 - (1 + RingRadius[RNUM - 1] / rd) / exp(RingRadius[RNUM - 1] / rd);

	// FractionRings * RNUM gives the fraction in a ring as a function of the mass in a ring without * RNUM gives the
	// fraction in a ring as a function of the total mass

	// Normalized fraction of mass in each ring.
	fractionRings[0] = (1 - (1 + RingRadius[0] / rd) / exp(RingRadius[0] / rd)) / ringtot;
	for ( j = 1; j < RNUM; j++ ) {
		fractionRings[j] = ((1 + RingRadius[j - 1] / rd) / exp(RingRadius[j - 1] / rd) -
							(1 + RingRadius[j] / rd) / exp(RingRadius[j] / rd)) / ringtot;
	}

	for ( j = 0; j < RNUM; j++ ) {
		if ( fractionRings[j] < 0. ) {
			fractionRings[j] = 0.;
		}
	}

	for ( j              = 0; j < RNUM; j++ ) {
		Gal[p].ColdGasRings[j]       = Gal[p].ColdGas * fractionRings[j];
		Gal[p].MetalsColdGasRings[j] = metals_add(metals_init(), Gal[p].MetalsColdGas, fractionRings[j]);

#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef INDIVIDUAL_ELEMENTS //All: [H][He][Cb][N][O][Ne][Mg][Si][S][Ca][Fe] or //Only [H][He][O][Mg][Fe]
		int kk;
		for ( kk = 0; kk < NUM_ELEMENTS; kk++ ) {
			Gal[p].ColdGasRings_elements[j][kk] = Gal[p].ColdGas_elements[kk] * fractionRings[j];
		}
#endif // INDIVIDUAL_ELEMENTS
#endif // DETAILED_METALS_AND_MASS_RETURN
		fractionRings[j] = 1.;
	}
	transfer_material_with_rings(t, "ColdGas", p, "ColdGas", fractionRings, "model_mergers.c", __LINE__);
#else // H2_AND_RINGS
	transfer_material (t, "ColdGas", p, "ColdGas", 1., "model_mergers.c", __LINE__);
#endif // H2_AND_RINGS
	transfer_material (t, "HotGas", p, "HotGas", 1., "model_mergers.c", __LINE__);
	transfer_material (t, "EjectedMass", p, "EjectedMass", 1., "model_mergers.c", __LINE__);

	mass_checks (p, "model_mergers.c", __LINE__);
	mass_checks (t, "model_mergers.c", __LINE__);

#ifdef TRACK_MASSGROWTH_CHANNELS
	Gal[t].MassFromMergers += Gal[p].DiskMass + Gal[p].BulgeMass;
	Gal[p].MassFromInSitu  = 0.;
	Gal[p].MassFromMergers = 0.;
	Gal[p].MassFromBursts  = 0.;

#ifdef STAR_FORMATION_HISTORY
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
	for ( ii = 0; ii <= Gal[t].sfh_ibin; ii++ ) {
		Gal[t].sfh_MassFromMergers[ii] += Gal[t].sfh_DiskMass[ii] + Gal[t].sfh_BulgeMass[ii];
	}

	for ( ii = 0; ii <= Gal[p].sfh_ibin; ii++ ) {
		Gal[p].sfh_MassFromInSitu[ii] = 0.;
	}
	for ( ii = 0; ii <= Gal[p].sfh_ibin; ii++ ) {
		Gal[p].sfh_MassFromMergers[ii] = 0.;
	}
	for ( ii             = 0; ii <= Gal[p].sfh_ibin; ii++ ) {
		Gal[p].sfh_MassFromBursts[ii] = 0.;
	}
#endif // TRACK_SFH_MASSGROWTH_CHANNELS
#endif // STAR_FORMATION_HISTORY
#endif // TRACK_MASSGROWTH_CHANNELS

#ifdef TRACK_BURST
	// The whole burst component gets transferred.
	transfer_material(t, "BurstMass", p, "BurstMass", 1., "model_mergers.c", __LINE__);
#endif // TRACK_BURST
#ifdef H2_AND_RINGS
	for ( j              = 0; j < RNUM; j++ ) {
		fractionRings[j] = 1.;
	}
#endif // H2_AND_RINGS

	//Bulge occupies the same place as the disk it forms from, after the merger is finished the material will be
	// distributed in a Jaffe profile
#ifdef H2_AND_RINGS
	if ( BulgeFormationInMinorMergersOn ) {
		transfer_material_with_rings(t, "BulgeMass", p, "DiskMass", fractionRings, "model_mergers.c", __LINE__);
	}
	else {
		transfer_material_with_rings(t, "DiskMass", p, "DiskMass", fractionRings, "model_mergers.c", __LINE__);
	}

	if ( Gal[p].BulgeMass > 0. )
#ifdef RINGS_IN_BULGES
		transfer_material_with_rings(t, "BulgeMass", p, "BulgeMass", fractionRings, "model_mergers.c", __LINE__);
#else // RINGS_IN_BULGES
	transfer_material(t, "BulgeMass", p, "BulgeMass", 1., "model_mergers.c", __LINE__);
#endif // RINGS_IN_BULGES

#else // H2_AND_RINGS

#ifdef DI_MERGERS
	// Minor mergers
	if ( mass_ratio < ThreshMajorMerger ) {
		// on bulge-dominated galaxies
		if ( tbulge > 0.5 * tmass_total ) {
			transfer_material(t, "BulgeMass", p, "DiskMass", 1., "model_mergers.c", __LINE__);
			if ( pbulge > 0. ) {
				transfer_material(t, "BulgeMass", p, "BulgeMass", 1., "model_mergers.c", __LINE__);
			}
		}// on disk-dominated galaxies
		else if ( tdisk > 0.5 * tmass_total ) {
			transfer_material(t, "DiskMass", p, "DiskMass", 1., "model_mergers.c", __LINE__);
			if ( pbulge > 0. ) {
				transfer_material(t, "DiskMass", p, "BulgeMass", 1., "model_mergers.c", __LINE__);
			}
		}// on gas-dominated galaxies
		else {
			if ( tdisk > tbulge ) {
				transfer_material(t, "DiskMass", p, "DiskMass", 1., "model_mergers.c", __LINE__);
				if ( pbulge > 0. ) {
					transfer_material(t, "DiskMass", p, "BulgeMass", 1., "model_mergers.c", __LINE__);
				}
			}
			if ( tbulge > tdisk ) {
				transfer_material(t, "BulgeMass", p, "DiskMass", 1., "model_mergers.c", __LINE__);
				if ( pbulge > 0. ) {
					transfer_material(t, "BulgeMass", p, "BulgeMass", 1., "model_mergers.c", __LINE__);
				}
			}
		}
		// Count the number of minor mergers.
		Gal[t].NMinorMergers += 1;
		Gal[t].NMinorMergers += Gal[p].NMinorMergers;
	}
	// Major mergers
	if ( mass_ratio > ThreshMajorMerger ) {
		transfer_material(t, "BulgeMass", p, "DiskMass", 1., "model_mergers.c", __LINE__);
		if ( pbulge > 0.0 ) {
			transfer_material(t, "BulgeMass", p, "BulgeMass", 1., "model_mergers.c", __LINE__);
		}
		// Count the number of major mergers.
		Gal[t].NMajorMergers += 1;
		Gal[t].NMajorMergers += Gal[p].NMajorMergers;
	}

	// If the bulge mass of the central galaxy is still 0.0 set it to something tiny in order to avoid the
	// termination of you code (line 333).
	if ( Gal[t].BulgeMass == 0.0 ) {
		Gal[t].BulgeMass = 1e-10;
	}
#else // DI_MERGERS
	if (BulgeFormationInMinorMergersOn)
	{
		transfer_material (t, "BulgeMass", p, "DiskMass", 1., "model_mergers.c", __LINE__);
	}
	else
	{
		transfer_material (t, "DiskMass", p, "DiskMass", 1., "model_mergers.c", __LINE__);
	}
	if (Gal[p].BulgeMass > 0.)
	{
		transfer_material (t, "BulgeMass", p, "BulgeMass", 1., "model_mergers.c", __LINE__);
	}

#endif // DI_MERGERS
#endif // H2_AND_RINGS

	transfer_material (t, "ICM", p, "ICM", 1., "model_mergers.c", __LINE__);

	Gal[t].BlackHoleMass += Gal[p].BlackHoleMass;
	Gal[p].BlackHoleMass = 0.;
	Gal[t].BlackHoleGas += Gal[p].BlackHoleGas;
	Gal[p].BlackHoleGas  = 0.;
	Gal[t].StarMerge += Gal[p].StarMerge;
	Gal[p].StarMerge     = 0.;

	mass_checks (p, "model_mergers.c", __LINE__);
	mass_checks (t, "model_mergers.c", __LINE__);

	// Update the gas spin - gasdiskradius updated in the end of deal_with_galaxy_mergers
	for (ii = 0; ii < 3; ii ++)
	{
		tspin[ii] = Gal[t].ColdGasSpin[ii] * tmass + Gal[t].HaloSpin[ii] * pmass;
	}
	if (Gal[t].ColdGas != 0)
	{
		for (ii = 0; ii < 3; ii ++)
		{
			Gal[t].ColdGasSpin[ii] = tspin[ii] / (Gal[t].ColdGas);
		}
	}

	// Bulge spin
#ifdef DI_BULGESPIN

	if ( BulgeFormationInMinorMergersOn )
		ptrans = pdisk + pbulge;
	else
		ptrans = pbulge;
	for ( j    = 0; j < 3; j ++ ) {
		tspin[j] = Gal[t].BulgeSpin[j] * tbulge + Gal[t].HaloSpin[j] * ptrans;
	}
	if ( Gal[t].BulgeMass != 0 ) {
		for ( j = 0; j < 3; j ++ ) {
			Gal[t].BulgeSpin[j] = tspin[j] / Gal[t].BulgeMass;
		}
	}
#endif // DI_BULGESPIN
	Gal[t].Sfr += Gal[p].Sfr;
	if (BulgeFormationInMinorMergersOn)
	{
		Gal[t].SfrBulge += Gal[p].Sfr;
	}

	// Add the black hole accretion rates. This makes little sense but is not used if the superior BlackHoleGrowth==1
	// switch is on.
	Gal[t].QuasarAccretionRate += Gal[p].
			                                    QuasarAccretionRate;
	Gal[t].RadioAccretionRate += Gal[p].
			                                   RadioAccretionRate;

	for (
			outputbin = 0;
			outputbin < NOUT;
			outputbin ++)
	{
		Gal[t].MassWeightAge[outputbin] += Gal[p].MassWeightAge[outputbin];
	}
}
// end add_galaxies_together



















/** @brief In a major merger, both disks are destroyed and all the mass transferred to the bulge. The galaxies have
 * already been merged, so all we need to do here is transfer stars from disk to bulge of the central galaxy.           */
void make_bulge_from_burst (int p)
{
#ifdef DI_BULGESPIN
	int   j;
	float pbulge, pdisk, pspin[3];

	// Record masses before mass transfer
	pbulge = Gal[p].BulgeMass;
	pdisk  = Gal[p].DiskMass;
#endif // DI_BULGESPIN

#ifdef H2_AND_RINGS
	int    jj;
	double fractionRings[RNUM];
	for ( jj = 0; jj < RNUM; jj ++ ) {
		fractionRings[jj] = 1.;
	}
	transfer_material_with_rings(p, "BulgeMass", p, "DiskMass", fractionRings, "model_mergers.c", __LINE__);
#else // H2_AND_RINGS
	transfer_material (p, "BulgeMass", p, "DiskMass", 1., "model_mergers.c", __LINE__);
#endif // H2_AND_RINGS
	mass_checks (p, "model_mergers.c", __LINE__);

#ifdef DIBULGESPIN
	// Momentum conservation: As all the disk is moved to the bulge, we must assume that it carries its angular
	// momentum with it.
	for ( j = 0; j < 3; j ++ ) {
		pspin[j] = Gal[p].BulgeSpin[j] * pbulge + Gal[p].DiskSpin[j] * pdisk;
	}
	if ( Gal[p].BulgeMass != 0 ) {
		for ( j = 0; j < 3; j ++ ) {
			Gal[p].BulgeSpin[j] = pspin[j] / Gal[p].BulgeMass;
		}
	}
#endif // DI_BULGESPIN

	// Update the star formation rate.
	Gal[p].SfrBulge = Gal[p].Sfr;
#ifdef H2_AND_RINGS
	for(jj=0;jj<RNUM;jj++) Gal[p].SfrRings[jj]=0;
#endif
}
// end make_bulge_from_burst



















/** @brief If StarBurstModel = 0 (since Croton2006), the Somerville 2001 model of bursts is used. The burst can happen
 * for both major and minor mergers, with a fraction of the added cold gas from the satellite and central being
 * consumed. SN Feedback from starformation is computed and the sizes of bulge and disk followed (not done for the
 * other burst mode).The coefficients in eburst are taken from TJ Cox's PhD thesis and should be more accurate than
 * previous.                                                                                                             */
double
collisional_starburst_recipe (double mass_ratio, int merger_centralgal, int centralgal, double time, double deltaT)
{
	double mstars, eburst, Ggas;
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
	double metallicitySF;
#endif // POST_PROCESS_MAGS
#endif // COMPUTE_SPECPHOT_PROPERTIES
#ifdef H2_AND_RINGS
	double mstarsRings[RNUM];
	int    j;
#endif // H2_AND_RINGS

	Ggas = Gal[merger_centralgal].ColdGas;

	// The bursting fraction given the mass ratio.
	eburst = SfrBurstEfficiency * pow (mass_ratio, SfrBurstSlope);
	mstars = eburst * Gal[merger_centralgal].ColdGas;

	if (mstars < 0.0)
	{
		mstars = 0.0;
	}
#ifdef H2_AND_RINGS
	for ( j    = 0; j < RNUM; j++ ) {
		mstarsRings[j] = eburst * Gal[merger_centralgal].ColdGasRings[j];
		if ( mstarsRings[j] < 0.0 ) {
			mstarsRings[j] = 0.0;
		}
	}
#endif // H2_AND_RINGS

	// Otherwise there is another check inside SN_feedback.
#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
	if ( mstars > Gal[merger_centralgal].ColdGas ) {
		mstars = Gal[merger_centralgal].ColdGas;
	}
#ifdef H2_AND_RINGS
	for ( j    = 0; j < RNUM; j++ ) {
		if ( mstarsRings[j] > Gal[merger_centralgal].ColdGasRings[j] ) {
			mstarsRings[j] = Gal[merger_centralgal].ColdGasRings[j];
		}
	}
#endif // H2_AND_RINGS
#endif // FEEDBACK_COUPLED_WITH_MASS_RETURN

	// Store the value of the metallicity of the cold phase when SF occurs. Used to update luminosities below.
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
	metallicitySF = metals_total(Gal[merger_centralgal].MetalsColdGas) / Gal[merger_centralgal].ColdGas;
#endif // POST_PROCESS_MAGS
#endif // COMPUTE_SPECPHOT_PROPERTIES

	// If FEEDBACK_COUPLED_WITH_MASS_RETURN feedback happens only when stars die, there is no need to balance it with SF
#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
	if (mstars > 0.)
#ifndef H2_AND_RINGS
	{
		update_stars_due_to_reheat (merger_centralgal, centralgal, &mstars);
	}
#else // H2_AND_RINGS
	update_stars_due_to_reheat(merger_centralgal, centralgal, &mstars, mstarsRings);
#endif // H2_AND_RINGS
#endif //FEEDBACK_COUPLED_WITH_MASS_RETURN

	// Update the star formation rate
	Gal[merger_centralgal].Sfr += mstars / deltaT;
#ifdef H2_AND_RINGS
	for ( j = 0; j < RNUM; j++ ) {
		Gal[merger_centralgal].SfrRings[j] += mstarsRings[j] / deltaT;
	}
#endif // H2_AND_RINGS
	mass_checks (merger_centralgal, "model_mergers.c", __LINE__);

	// Update_from_star_formation can only be called after SN_feeedback recipe since stars need to be re_set once the
	// reheated mass is known (star formation and feedback share the same fraction of cold gas).
	int nstep = - 1;
#ifndef H2_AND_RINGS
	if (mstars > 0.)
	{
		update_from_star_formation (merger_centralgal, mstars, "merger", nstep);
	}
#else // H2_AND_RINGS
	if ( mstars > 0. ) {
		update_from_star_formation(merger_centralgal, mstars, mstarsRings, "merger", nstep);
	}
#endif // H2_AND_RINGS

	mass_checks (merger_centralgal, "model_mergers.c", __LINE__);

	update_massweightage (merger_centralgal, mstars, time);

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
#ifndef H2_AND_RINGS
	if (mstars > 0.)
	{
		SN_feedback (merger_centralgal, centralgal, mstars, "ColdGas");
	}
#else // H2_AND_RINGS
	if ( mstars > 0. ) {
		SN_feedback(merger_centralgal, centralgal, mstars, mstarsRings, "ColdGas");
	}
#endif // H2_AND_RINGS
#endif // FEEDBACK_COUPLED_WITH_MASS_RETURN

	// Update the luminosities due to the stars formed.
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
	if ( mstars > 0.0 ) {
		add_to_luminosities(merger_centralgal, mstars, time, deltaT / STEPS, metallicitySF);
	}
#endif // POST_PROCESS_MAGS
#endif // COMPUTE_SPECPHOT_PROPERTIES

	if (Ggas > 0.)
	{
		return mstars / Ggas;
	}
	else
	{
		return 0.0;
	}
}
// end collisional_starburst_recipe



















/** @brief Calculates the bulge size after a merger for any type of merger. We need the gas and stellar mass values
 * before the merger since we have already invoked the add_galaxies_together function.                                  */
void
size_of_merger_remnants (double mass_ratio, int merger_centralgal, double MstarCentral, double MbulgeCentral,
                         double MgasCentral, double MstarSat, double MgasSat, double frac, double HMRCentral,
                         double HMRSat)
{

	double Efinal, Einitial, Eorbital, Eradiative, Rcen;
	float  tiny = 1.0e-8;
	double Crad = 0.1;

	// On a major merger both the progenitors's total stellar masses and stars formed during merger
	// (collisional_starburst_recipe) form the final bulge (aka elliptical galaxy), therefore only these four stellar
	// components are used to compute the radius and masses.
	if (mass_ratio > ThreshMajorMerger)
	{

		// Calculate the energy components.
		Efinal     = (MstarCentral + frac * MgasCentral + MstarSat + frac * MgasSat) *
		             (MstarCentral + frac * MgasCentral + MstarSat + frac * MgasSat);
		Einitial   = (MstarCentral + MgasCentral) * (MstarCentral + MgasCentral) / HMRCentral +
		             (MstarSat + MgasSat) * (MstarSat + MgasSat) / HMRSat;
		Eorbital   = (MstarCentral + MgasCentral) * (MstarSat + MgasSat) / (HMRCentral + HMRSat);
		Eradiative = Crad * Einitial * (MgasCentral + MgasSat) / (MstarCentral + MgasCentral + MstarSat + MgasSat);

		Gal[merger_centralgal].BulgeSize = Efinal / (Einitial + Eorbital + Eradiative);
	}

	else
	{
		// In a minor merger the total stellar mass of the satellite galaxy is moved to the bulge of the central galaxy,
		// therefore only these two stellar components are used to compute the radius and masses.
		Rcen = Gal[merger_centralgal].BulgeSize;

		// Calculate the energy components.
		Efinal     = (MbulgeCentral + MstarSat) * (MbulgeCentral + MstarSat);
		Einitial   = MbulgeCentral * MbulgeCentral / Rcen + MstarSat * MstarSat / HMRSat;
		Eorbital   = MbulgeCentral * MstarSat / (Rcen + HMRSat);
		Eradiative = Crad * Einitial * (MgasCentral + MgasSat) / (MstarCentral + MgasCentral + MstarSat + MgasSat);

		Gal[merger_centralgal].BulgeSize = Efinal / (Einitial + Eorbital + Eradiative);
	}

	if (Gal[merger_centralgal].BulgeSize < tiny)
	{
		Gal[merger_centralgal].BulgeSize = tiny;
	}
}
// end size_of_merger_remnants



















/** @brief Calculates the bulge size after a merger. For any type of merger calculates the new bulge size using Eq. 33
 * in Guo2010:
 *
 *         \f$C\frac{GM^2_{\rm{new,bulge}}}{R_{\rm{new,bulge}}}=
 *             C\frac{GM^2_1}{R_1}+C\frac{GM^2_2}{R_2} + \alpha_{\rm{inter}}\frac{GM_1M_2}{R_1+R_2}\f$.
 *
 * This implementation assumed that the new bulge occupies the same space as the components that formed it.           */
void
bulgesize_from_merger (double mass_ratio, int merger_centralgal, int p, double MstarCentral, double MbulgeCentral,
                       double MgasCentral, double MstarSat, double MbulgeSat, double MgasSat, double frac,
                       double RgasCentral, double RStellarDiskCentral, double RgasSat, double RStellarDiskSat)
{
	double Mcen, Rcen, Msat, Rsat;
	double aint = 0.5;
	double c    = 0.5;

	// In a minor merger only the stars of the satellite galaxy are moved to the bulge of the central galaxy, therefore
	// only stellar components are used to compute radius and masses. In a minor merger only consider the bulge mass of
	// the central galaxy.
	if (mass_ratio < ThreshMajorMerger)
	{
		Mcen = MbulgeCentral;
		Rcen = Gal[merger_centralgal].BulgeSize;
		if (BulgeFormationInMinorMergersOn)
		{
			Msat = MstarSat;
		}
		else
		{
			Msat = MbulgeSat;
		}
		Rsat = (RStellarDiskSat * 1.68 * (MstarSat - MbulgeSat) + Gal[p].BulgeSize * MbulgeSat) / MstarSat;
	}
		// On a major merger both stellar and gas (after a burst) components form the final bulge and need to be
		// considered.
	else
	{
		Mcen = MstarCentral + frac * MgasCentral;
		Rcen = (RStellarDiskCentral * 1.68 * (MstarCentral - MbulgeCentral) +
		        Gal[merger_centralgal].BulgeSize * MbulgeCentral + RgasCentral * frac * MgasCentral * 1.68) /
		       (MgasCentral * frac + MstarCentral);

		Msat = MstarSat + frac * MgasSat;
		Rsat = (RStellarDiskSat * 1.68 * (MstarSat - MbulgeSat) + Gal[p].BulgeSize * MbulgeSat +
		        RgasSat * frac * MgasSat * 1.68) / (MgasSat * frac + MstarSat);
	}

	float tiny;
	tiny = 1.e-8;
	if (Rsat > 0. && Rsat < tiny)
	{
		Rsat = tiny;
	}
	if (Rcen > 0. && Rcen < tiny)
	{
		Rcen = tiny;
	}

	// If both original radius are bigger than 0 then this is Eq. 33 in Guo 2010 solved for R_new,bulge with all terms
	// divided by G and C. The problem with these conditions is that, once we have an erroneously small value of
	// BulgeSize, we will keep it.
	if (Rcen >= tiny && Rsat >= tiny)
	{
		Gal[merger_centralgal].BulgeSize =
				(Msat + Mcen) * (Msat + Mcen) / (Msat * Msat / Rsat + Mcen * Mcen / Rcen +
				                                 aint / c * Msat * Mcen / (Rsat + Rcen));
	}
	else if (Rcen >= tiny)
	{
		Gal[merger_centralgal].BulgeSize =
				(Msat + Mcen) * (Msat + Mcen) / (Mcen * Mcen / Rcen + aint / c * Msat * Mcen / (Rsat + Rcen));
	}
	else if (Rcen <= tiny && Rsat >= tiny)
	{
		Gal[merger_centralgal].BulgeSize =
				(Msat + Mcen) * (Msat + Mcen) / (Msat * Msat / Rsat + aint / c * Msat * Mcen / (Rsat + Rcen));
	}
	else
	{
		Gal[merger_centralgal].BulgeSize = 0.0;
	}

	if ((Msat + Mcen > 0.0 && Gal[merger_centralgal].BulgeSize == 0.0) ||
	    (Msat + Mcen == 0.0 && Gal[merger_centralgal].BulgeSize > 0.0))
	{
		char sbuf[1000];
		printf ("halonr=%d, merger_centralgal %d\n\n", Gal[merger_centralgal].HaloNr, merger_centralgal);
		printf ("New Bulge Mass from Central (Mcen)=%e\n New Bulge Mass from Satellite (Msat)=%e\n NewBulge size=%e\n\n",
		        Msat, Mcen, Gal[merger_centralgal].BulgeSize);
		if (mass_ratio < ThreshMajorMerger)
		{
			printf ("minor merger, new mass=original mass\n");
			printf ("New Bulge Mass From Central (Mcen)   = MbulgeCentral = %f\n", MbulgeCentral);
			printf ("New Bulge Mass From Satellite (Msat) = MstarSat  = %f\n", MstarSat);
		}
		else
		{
			printf ("New Bulge From Central (Mcen)   = MstarCentral+frac*MgasCentral = %f+%f*%f\n", MstarCentral,
			        frac,
			        MgasCentral);
			printf ("New Bulge From Satellite (Msat) = MstarSat+frac*MgasSat = %f+%f*%f\n", MstarSat, frac, MgasSat);
		}
		printf ("BulgeSize from Central (Rcen)=%e\nBulgeSize from Satellite (Rsat)=%e\nmass ratio=%f\n\n", Rcen,
		        Rsat,
		        mass_ratio);
		printf ("the following masses don't tell a lot because things have been merged already!!!\n");
		printf ("    sat: BulgeMass=%0.7f, BulgeSize=%0.7f, GasMass=%0.7f, GasSize=%0.7f, DiskMass=%0.7f StellarSize=%0.7f \n",
		        Gal[p].BulgeMass, Gal[p].BulgeSize, Gal[p].ColdGas, Gal[p].ColdGasRadius, Gal[p].DiskMass,
		        Gal[p].DiskRadius);
		printf ("central: BulgeMass=%0.7f, BulgeSize=%0.7f, GasMass=%0.7f, GasSize=%0.7f, DiskMass=%0.7f StellarSize=%0.7f \n",
		        Gal[merger_centralgal].BulgeMass, Gal[merger_centralgal].BulgeSize, Gal[merger_centralgal].ColdGas,
		        Gal[merger_centralgal].ColdGasRadius, Gal[merger_centralgal].DiskMass,
		        Gal[merger_centralgal].DiskRadius);
		sprintf (sbuf, "\n bulgesize wrong in merger");
		terminate(sbuf);
		exit (0);
	}
}
// end bulgesize_from_merger
