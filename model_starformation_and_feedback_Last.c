#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

/** @file recipe_starformation_and_feedback.c
 *  @brief recipe_starformation_and_feedback.c computes the amount of stars formed from the cold gas, the amount of gas
 *  reheated from cold to hot and the amount of gas ejected from hot to external.
 *
 * The routine is divided in two parts, star formation and SN feedback, with a number of different implementations
 * controlled by input parameters.
 *  0 -\f$M_{\rm{crit}}=3.8\times 10^9
 *     \left(\frac{V_{\rm{max}}}{200\,\rm{km s}^{-1}}\right)
 *     \left(\frac{r_{\rm{disk}}}{10\,\rm{kpc}}\right)M_{\odot}\f$
 *     (Eq. 16 Guo2010) (StarFormationModel = 0), \n
 *        - same as 0 but using \f$V_{\rm{max}}\f$ or \f$V_{\rm{max,infall}}\f$
 *          instead of \f$V_{\rm{vir}}\f$ and allowing SF in satellites. *
 *
 * There are 2 options for the <B>SN Feedback Recipe</B>:
 * 0 - \f$\epsilon_{\rm{disk}}=\epsilon
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_1}\biggr]\f$,
 *     \f$\epsilon_{\rm{halo}}=\eta
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_2}\biggr]\f$
 *     (Eqs. 19 & 21 Guo2010)(FeedbackReheatingModel = 0)
 *
 * Also, Guo2010 alowed for type 1 satellite to have gas cycles and receive gas from their own satellites when these
 * are outside Rvir of the type 0.    																			 	  */



















/** @brief Main recipe, calculates the fraction of cold gas turned into stars due to star formation; the fraction of
 * mass instantaneously recycled and returned to the cold gas; the fraction of gas reheated from cold to hot, ejected
 * from hot to external and returned from ejected to hot due to SN feedback.  										  */
void starformation (int p, int centralgal, double time, double dt, int nstep) {
	double tdyn, strdot = 0., stars, cold_crit;

#ifdef H2_AND_RINGS
	double strdotr[RNUM], starsRings[RNUM], sfe, cold_crit_rate, SigmaGas, SigmaGasratio;
	int j;
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
	double metallicitySF;
#endif
#endif

	double Vmax, gas_radius;

	if (Gal[p].Type == 0) {
		Vmax = Gal[p].Vmax;
	} else {
		Vmax = Gal[p].InfallVmax;
	}
	gas_radius          = Gal[p].ColdGasRadius;
	tdyn                = gas_radius / Vmax;
	cold_crit           = SfrColdCrit * Vmax / 200. * gas_radius * 100.;

	// Standard star formation law (Croton2006, Delucia2007, Guo2010, Henriques2015)
	if (StarFormationModel == 0) {
		if (Gal[p].ColdGas > cold_crit) {
			strdot = SfrEfficiency * (Gal[p].ColdGas - cold_crit) / tdyn;
		} else {
			strdot = 0.0;
		}
	}

#ifdef H2_AND_RINGS
	update_h2fraction(p);

	sfe     = SfrEfficiency * UnitTime_in_years / Hubble_h; // convert from yr-1 into code units of time
	if ( SFRtdyn == 1 ) {
		sfe = (sfe / tdyn) / UnitTime_in_years * Hubble_h; // for star formation rate proportional to 1/t_dyn
	}
	for ( j = 0; j < RNUM; j ++ ) {
		if ( StarFormationModel == 0 ) {
			strdotr[j] = strdot * Gal[p].ColdGasRings[j] / Gal[p].ColdGas;
		}
		else if ( StarFormationModel == 2 ) {
			if ( Gal[p].Type == 0 ) {
				cold_crit_rate = SfrColdCrit * Gal[p].Vmax / 200. * Gal[p].ColdGasRadius / Gal[p].ColdGas * 100.;
			}
			else {
				cold_crit_rate = SfrColdCrit * Gal[p].InfallVmax / 200. * Gal[p].ColdGasRadius / Gal[p].ColdGas * 100.;
			}
			if ( cold_crit_rate < 1 && cold_crit_rate >= 0 ) {
				strdotr[j] = sfe * Gal[p].ColdGasRings[j] * (1 - cold_crit_rate);
			}
			else {
				strdotr[j] = 0.0;
			}
		}
		else if ( StarFormationModel == 3 ) // star formation law in Krumholz et al. 2009
		{
			double SigmaGas0 = 85.0, SF_Law_pow = 0.33;

			if ( j == 0 ) {
				SigmaGas = Gal[p].ColdGasRings[j] / (M_PI * RingRadius[j] * RingRadius[j]) / WARM_PHASE_FACTOR *
						   Clumpingfactor;
			}
			else {
				SigmaGas = Gal[p].ColdGasRings[j] /
						   (M_PI * (RingRadius[j] * RingRadius[j] - RingRadius[j - 1] * RingRadius[j - 1])) /
						   WARM_PHASE_FACTOR * Clumpingfactor;
			}

			SigmaGas      = SigmaGas * 0.01 * Hubble_h; // convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2)
			SigmaGasratio = pow(SigmaGas / SigmaGas0, SF_Law_pow);

			if ( SigmaGasratio < 1.0 && SigmaGasratio > 0.0 ) {
				strdotr[j] =
						sfe / SigmaGasratio * Gal[p].ColdGasRings[j] * Gal[p].H2fractionRings[j] / WARM_PHASE_FACTOR;
			} // only cold H2 component is proportional to star formation rate.
			else {
				if ( SigmaGasratio >= 1.0 ) {
					strdotr[j] =
							sfe * SigmaGasratio * Gal[p].ColdGasRings[j] * Gal[p].H2fractionRings[j] /
							WARM_PHASE_FACTOR;
				}
				else {
					strdotr[j] = 0.0;
				}
			}
		}
		else if ( StarFormationModel == 4 )    // The star formation law in Fu et al. 2010
		{
			if ( Gal[p].H2fractionRings[j] >= 0.0 ) {
				strdotr[j] = sfe * Gal[p].ColdGasRings[j] * Gal[p].H2fractionRings[j] / WARM_PHASE_FACTOR;
			} // only cold H2 component is proportional to star formation rate.
			else {
				strdotr[j] = 0.0;
			}
		}
		else {
			strdotr[j] = 0.0;
		}
	}
	for ( j = 0, strdot = 0; j < RNUM; j ++ ) {
		strdot += strdotr[j];
	}
#endif // H2_AND_RINGS

	stars = strdot * dt; // units of dynamical time are Mpc/Km/s - no conversion on dt needed 3.06e19 to 3.15e19

#ifdef H2_AND_RINGS
	for ( j = 0; j < RNUM; j++ ) {
		if ( strdotr[j] < 0.0 ) {
			strdotr[j]    = 0.;
			starsRings[j] = strdotr[j] * dt;
		}
		if ( starsRings[j] < 0.0 ) {
			starsRings[j] = 0.0;
		}
	}
#endif // H2_AND_RINGS

	// Otherwise this check is done inside update_stars_due_to_reheat for stars + reheated mass !> cold gas.
#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
	if ( stars > Gal[p].ColdGas ) {
		stars = Gal[p].ColdGas;
	}
#ifdef H2_AND_RINGS
	for ( j = 0; j < RNUM; j++ ) {
		if ( starsRings[j] > Gal[p].ColdGasRings[j] ) {
			starsRings[j] = Gal[p].ColdGasRings[j];
		}
	}
#endif // H2_AND_RINGS
#endif // FEEDBACK_COUPLED_WITH_MASS_RETURN

	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
	mass_checks (centralgal, "model_starformation_and_feedback.c", __LINE__);

	// Store the value of the metallicity of the cold phase when SF occurs.
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
	if ( Gal[p].ColdGas > 0. ) {
		metallicitySF = metals_total(Gal[p].MetalsColdGas) / Gal[p].ColdGas;
	}
	else {
		metallicitySF = 0.;
	}
#endif // POST_PROCESS_MAGS
#endif // COMPUTE_SPECPHOT_PROPERTIES

	//if FEEDBACK_COUPLED_WITH_MASS_RETURN feedback happens only when stars die, there is no need to balance it with SF
#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
	if (stars > 0.) {
#ifndef H2_AND_RINGS
		update_stars_due_to_reheat (p, centralgal, &stars);
#else // H2_AND_RINGS
		update_stars_due_to_reheat(p, centralgal, &stars, starsRings);
#endif // H2_AND_RINGS
	}
#endif //FEEDBACK_COUPLED_WITH_MASS_RETURN

	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
	mass_checks (centralgal, "model_starformation_and_feedback.c", __LINE__);

	// Update the star formation rate.
	Gal[p].Sfr += stars / (dt * STEPS);

#ifdef H2_AND_RINGS
	for(j=0;j<RNUM;j++)
	{
	  Gal[p].SfrRings[j] += starsRings[j] / (dt * STEPS);
	}
#endif // H2_AND_RINGS

	// update_from_star_formation can only be called after SD_feeedback recipe since stars need to be re_set once the
	// reheated mass is known (star formation and feedback share the same fraction of cold gas).
	if (stars > 0.0) {
#ifndef H2_AND_RINGS
		update_from_star_formation (p, stars, "insitu", nstep); // false indicates not a burst
#else // H2_AND_RINGS
		update_from_star_formation(p, stars, starsRings, "insitu", nstep); // false indicates not a burst
#endif // H2_AND_RINGS
	}

	update_massweightage (p, stars, time);

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
	if (stars > 0.0) {
#ifndef H2_AND_RINGS
		SN_feedback (p, centralgal, stars, "ColdGas");
#else // H2_AND_RINGS
		SN_feedback(p, centralgal, stars, starsRings, "ColdGas");
#endif // H2_AND_RINGS
	}
#endif // FEEDBACK_COUPLED_WITH_MASS_RETURN

	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
	mass_checks (centralgal, "model_starformation_and_feedback.c", __LINE__);

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS

	//  Update the luminosities due to the stars formed
	if ( stars > 0.0 ) {
		add_to_luminosities(p, stars, time, dt, metallicitySF);
	}
#endif // POST_PROCESS_MAGS
#endif // COMPUTE_SPECPHOT_PROPERTIES

	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
	mass_checks (centralgal, "model_starformation_and_feedback.c", __LINE__);

	if (DiskInstabilityModel == 0) {
		if (Gal[p].DiskMass > 0.0) {
			check_disk_instability (p, dt, time);
		}
	}

	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
	mass_checks (centralgal, "model_starformation_and_feedback.c", __LINE__);
}
// end starformation



















/** @brief 																											  */
#ifndef H2_AND_RINGS

void update_stars_due_to_reheat (int p, int centralgal, double*stars)
#else // H2_AND_RINGS
void update_stars_due_to_reheat( int p, int centralgal, double *stars, double starsRings[] )
#endif // H2_AND_RINGS
{
	double reheated_mass, frac, Radius_low = 0.;
#ifndef H2_AND_RINGS
	reheated_mass = compute_SN_reheat (p, centralgal, *stars, Gal[p].ColdGas, metals_total (Gal[p].MetalsColdGas),
	                                   Radius_low, Gal[p].ColdGasRadius);
	if ((*stars + reheated_mass) > Gal[p].ColdGas) {
		frac = Gal[p].ColdGas / (*stars + reheated_mass);
		*stars *= frac;
	}
#else // H2_AND_RINGS
	int jj;
	for ( jj = 0; jj < RNUM; jj++ ) {
		if ( jj > 0 ) {
			Radius_low = RingRadius[jj - 1];
		}
		reheated_mass = compute_SN_reheat(p, centralgal, starsRings[jj], Gal[p].ColdGasRings[jj],
										  metals_total(Gal[p].MetalsColdGasRings[jj]), Radius_low, RingRadius[jj]);
		if ((starsRings[jj] + reheated_mass) > Gal[p].ColdGasRings[jj] ) {
			frac = Gal[p].ColdGasRings[jj] / (starsRings[jj] + reheated_mass);
			starsRings[jj] *= frac;
		}
		*stars += starsRings[jj];
	}
#endif // H2_AND_RINGS
}
// end update_stars_due_to_reheat



















/** @brief Updates the different components due to star formation: mass and metals in stars and cold gas and stellar
 * spin. 																											  */
#ifndef H2_AND_RINGS

void update_from_star_formation (int p, double stars, char type_of_event[], int nstep)
#else // H2_AND_RINGS
void update_from_star_formation(int p, double stars, double starsRings[], char type_of_event[], int nstep)
#endif // H2_AND_RINGS
{
	int    ii;
	double stars_to_add = 0., NonRecycledFraction = 0.;
#ifdef H2_AND_RINGS
	int jj;
	double stars_to_addr[RNUM], fractionRings[RNUM];
#else // H2_AND_RINGS
	double fraction;
#endif // H2_AND_RINGS
	if (Gal[p].ColdGas <= 0. || stars <= 0.) {
		printf ("Gal[p].ColdGas <= 0. || stars <= 0., Coldgas=%0.5e stars=%0.5e, in function update_from_star_formation, "
		        "model_starformation_and_feedback.c line:%d\n", Gal[p].ColdGas, stars, __LINE__);
		exit (0);
	}

	// If DETAILED_METALS_AND_MASS_RETURN, no longer an assumed instantaneous recycled fraction. Mass is returned over
	// time via SNe and AGB winds. Update the Stellar Spin when forming stars.
#ifndef DETAILED_METALS_AND_MASS_RETURN
	NonRecycledFraction = (1 - RecycleFraction);
#else // DETAILED_METALS_AND_MASS_RETURN
	NonRecycledFraction = 1.;
#endif // DETAILED_METALS_AND_MASS_RETURN

#ifndef H2_AND_RINGS
	stars_to_add = NonRecycledFraction * stars;
#else // H2_AND_RINGS
	for ( jj = 0; jj < RNUM; jj++ ) {
		stars_to_addr[jj] = NonRecycledFraction * starsRings[jj];
		stars_to_add += stars_to_addr[jj];
	}
#endif //H2_AND_RINGS

	if (Gal[p].DiskMass + stars_to_add > 1.e-8) {
		for (ii = 0; ii < 3; ii ++) {
			Gal[p].DiskSpin[ii] = ((Gal[p].DiskSpin[ii]) * (Gal[p].DiskMass) + stars_to_add * Gal[p].ColdGasSpin[ii]) /
			                      (Gal[p].DiskMass + stars_to_add);
		}
	}

	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);

	//  Update Gas and Metals from star formation.
#ifdef H2_AND_RINGS
	for ( jj = 0; jj < RNUM; jj++ ) {
		if ( Gal[p].ColdGasRings[jj] > 0. ) {
			fractionRings[jj] = stars_to_addr[jj] / Gal[p].ColdGasRings[jj];
		}
		else {
			fractionRings[jj] = 0.;
		}
	}
	transfer_material_with_rings(p, "DiskMass", p, "ColdGas", fractionRings, "model_starformation_and_feedback.c",
								 __LINE__);
#else // H2_AND_RINGS
	fraction = stars_to_add / Gal[p].ColdGas;
	transfer_material (p, "DiskMass", p, "ColdGas", fraction, "model_starformation_and_feedback.c", __LINE__);
#endif // H2_AND_RINGS


	// For this calculation we want just the long lived mass and take the instantaneous recycling approximation even
	// for the detailed chemical enrichment because it is not possible to know which component to eject mass from
	// afterwards.
#ifdef TRACK_MASSGROWTH_CHANNELS
	double long_lived_mass;
	long_lived_mass = stars_to_add;
#ifdef DETAILED_METALS_AND_MASS_RETURN
	long_lived_mass *= (1 - RecycleFraction);
#endif // DETAILED_METALS_AND_MASS_RETURN

	if ( strcmp(type_of_event, "insitu") == 0 ) {
		Gal[p].MassFromInSitu += long_lived_mass;
#ifdef STAR_FORMATION_HISTORY
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
		Gal[p].sfh_MassFromInSitu[Gal[p].sfh_ibin] += long_lived_mass;
#endif // STAR_FORMATION_HISTORY
#endif // TRACK_SFH_MASSGROWTH_CHANNELS
	}

	if ( strcmp(type_of_event, "merger") == 0 ) {
		Gal[p].MassFromBursts +=
				long_lived_mass;
#ifdef STAR_FORMATION_HISTORY
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
		Gal[p].sfh_MassFromBursts[Gal[p].sfh_ibin] += long_lived_mass;
#endif // STAR_FORMATION_HISTORY
#endif // TRACK_SFH_MASSGROWTH_CHANNELS
	}
#endif //TRACK_MASSGROWTH_CHANNELS

#ifdef TRACK_BURST
	if ( strcmp(type_of_event, "merger") == 0 ) {
		Gal[p].BurstMass += stars_to_add;
#ifdef STAR_FORMATION_HISTORY
		Gal[p].sfh_BurstMass[Gal[p].sfh_ibin] += stars_to_add;
#endif // STAR_FORMATION_HISTORY
	}
#endif // TRACK_BURST

	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);

	if (FeedbackReheatingModel == 0 || FeedbackReheatingModel == 1) {
		// stars (instead of star_to_add) used because the Yield is defined as a fraction of all stars formed, not just
		// long lived.
#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef METALS_SELF
		Gal[p].MetalsHotGasSelf.type2 += Yield * FracZtoHot * stars;
#endif // METALS_SELF
#else //DETAILED_METALS_AND_MASS_RETURN
		// This part is not used if OPT+=DELAYED_ENRICHMENT_AND MASS_RETURN as yield and recycling fraction are not
		// fixed.
#ifndef H2_AND_RINGS
		Gal[p].MetalsColdGas += Yield * (1. - FracZtoHot) * stars;
#else // H2_AND_RINGS
		for ( jj = 0; jj < RNUM; jj++ ) {
			Gal[p].MetalsColdGasRings[jj] += Yield * (1. - FracZtoHot) * starsRings[jj];
			Gal[p].MetalsColdGas += Yield * (1. - FracZtoHot) * starsRings[jj];
		}
#endif // H2_AND_RINGS
		Gal[Gal[p].CentralGal].MetalsHotGas += Yield * FracZtoHot * stars;
		Gal[Gal[p].CentralGal].HotGas += Yield * FracZtoHot * stars;
#ifdef METALS_SELF
		Gal[p].MetalsHotGasSelf += Yield * FracZtoHot * stars;
#endif // METALS_SELF
#endif //DETAILED_METALS_AND_MASS_RETURN
	}

	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);

	if (DiskRadiusModel == 0) {
		Gal[p].DiskRadius = get_stellar_disk_radius (p);
	}
}
// end update_from_star_formation



















/** @brief  There are two modes for supernova feedback corresponding to when the mass returning by dying stars is
 * returned to the cold gas - reheat and ejection; and when the mass is returned to the hot gas - only ejection.	  */
#ifndef H2_AND_RINGS

void SN_feedback (int p, int centralgal, double stars, char feedback_location[])
#else // H2_AND_RINGS
void SN_feedback(int p, int centralgal, double stars, double starsRings[], char feedback_location[])
#endif // H2_AND_RINGS
{
	double EjectVmax, EjectVvir, SN_Energy, Reheat_Energy, ReScaled_EnergySNcode, Radius_low;
	double reheated_mass = 0., ejected_mass = 0.;
	// SN FEEDBACK RECIPES
#ifdef H2_AND_RINGS
	double reheated_massr[RNUM];
	int jj;
#endif // H2_AND_RINGS

	Radius_low = 0.;
	//REHEAT
#ifndef H2_AND_RINGS
	//when FEEDBACK_COUPLED_WITH_MASS_RETURN some mass goes into HOTGAS and does not produce reheating
	if (strcmp (feedback_location, "HotGas") == 0) {
		reheated_mass = 0;
	} else {
		reheated_mass = compute_SN_reheat (p, centralgal, stars, Gal[p].ColdGas, metals_total (Gal[p].MetalsColdGas),
		                                   Radius_low, Gal[p].ColdGasRadius);
	}
#else // H2_AND_RINGS
	reheated_mass = 0.0;

	for ( jj = 0; jj < RNUM; jj++ ) {
		if ( jj > 0 ) {
			Radius_low = RingRadius[jj - 1];
		}
		//when FEEDBACK_COUPLED_WITH_MASS_RETURN some mass goes into HOTGAS and does not produce reheating
		if ( strcmp(feedback_location, "HotGas") == 0 ) {
			reheated_massr[jj] = 0.;
		}
		else {
			reheated_massr[jj] = compute_SN_reheat(p, centralgal, starsRings[jj], Gal[p].ColdGasRings[jj],
												   metals_total(Gal[p].MetalsColdGasRings[jj]), Radius_low,
												   RingRadius[jj]);
		}
		reheated_mass += reheated_massr[jj];
	}

	//reheated_mass > Gal[p].ColdGas might happen due to precision
	if ( reheated_mass > Gal[p].ColdGas ) {
		reheated_mass = Gal[p].ColdGas;
	}
#endif // H2_AND_RINGS

	// Determine ejection (for FeedbackEjectionModel 0 we have the dependence on Vmax) Guo2010 - eq 22
	// Note that satellites can now retain gas and have their own gas cycle.
	if (FeedbackEagleScaling == 1) {
		Radius_low            = 0.;
		ReScaled_EnergySNcode = EnergySNcode *
		                        EAGLE2015_rescale_of_EnergySN (Gal[p].ColdGas, metals_total (Gal[p].MetalsColdGas),
		                                                       Radius_low, Gal[p].ColdGasRadius);
	} else {
		ReScaled_EnergySNcode = EnergySNcode;
	}

	if (Gal[Gal[p].CentralGal].Type == 0) {
		EjectVmax = Gal[centralgal].Vmax;
		EjectVvir = Gal[centralgal].Vvir; // main halo Vvir
	} else {
		EjectVmax = Gal[Gal[p].CentralGal].InfallVmax;
		EjectVvir = Gal[Gal[p].CentralGal].Vvir; // central subhalo Vvir
	}

	if (FeedbackEjectionModel == 0) {
		ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * ReScaled_EnergySNcode) * stars *
		                min(1. / FeedbackEjectionEfficiency,
		                    .5 + 1 / pow (EjectVmax / EjectPreVelocity, EjectSlope)) -
		                reheated_mass * EjectVvir * EjectVvir) / (EjectVvir * EjectVvir);
	} else if (FeedbackEjectionModel == 1) //the ejected material is assumed to have V_SN
	{
		SN_Energy     = 0.5 * stars * (EtaSNcode * ReScaled_EnergySNcode);
		Reheat_Energy = 0.5 * reheated_mass * EjectVvir * EjectVvir;

		ejected_mass = (SN_Energy - Reheat_Energy) /
		               (0.5 * FeedbackEjectionEfficiency * (EtaSNcode * ReScaled_EnergySNcode));

		// If VSN^2<Vvir^2 nothing is ejected
		if (FeedbackEjectionEfficiency * (EtaSNcode * ReScaled_EnergySNcode) < EjectVvir * EjectVvir) {
			ejected_mass = 0.0;
		}
	}

	// Finished calculating mass exchanges, so just check that none are negative
	if (reheated_mass < 0.0) {
		reheated_mass = 0.0;
	}
	if (ejected_mass < 0.0) {
		ejected_mass = 0.0;
	}


	// Update For Feedback: update cold, hot, ejected gas fractions and respective metallicities there are a number of
	// changes introduced by Guo2010 concerning where the gas ends up.

	if (reheated_mass + ejected_mass > 0.) {
#ifndef H2_AND_RINGS
		update_from_feedback (p, centralgal, reheated_mass, ejected_mass);
#else // H2_AND_RINGS
		update_from_feedback(p, centralgal, reheated_mass, ejected_mass,  reheated_massr);
#endif // H2_AND_RINGS
	}

}
// end SN_feedback



















/** @brief In Guo2010 type 1s can eject, reincorporate gas and get gas from their own satellites (is not sent to the
 * type 0 galaxy as in Delucia2007), for gas flow computations: If satellite is inside Rvir of main halo, Vvir of main
 * halo used If it is outside, the Vvir of its central subhalo is used. 											  */
double compute_SN_reheat (int p, int centralgal, double stars, double ColdGas, double MetalsColdGas, double Radius_low,
                          double Radius_high) {
	double reheated_mass = 0., MergeCentralVvir = 0.;
	double ReScaled_EnergySNcode;

	if (FeedbackEagleScaling == 1) {
		ReScaled_EnergySNcode =
				EnergySNcode * EAGLE2015_rescale_of_EnergySN (ColdGas, MetalsColdGas, Radius_low, Radius_high);
	} else {
		ReScaled_EnergySNcode = EnergySNcode;
	}

	// Reheat
	if (ColdGas > 0.) {
		// Feedback depends on the circular velocity of the host halo Guo2010 - eq 18 & 19
		if (FeedbackReheatingModel == 0) {
			if (Gal[Gal[p].CentralGal].Type == 0) {
				reheated_mass = FeedbackReheatingEpsilon * stars *
				                (.5 + 1. / pow (Gal[Gal[p].CentralGal].Vmax / ReheatPreVelocity, ReheatSlope));
			} else {
				reheated_mass = FeedbackReheatingEpsilon * stars *
				                (.5 + 1. / pow (Gal[Gal[p].CentralGal].InfallVmax / ReheatPreVelocity, ReheatSlope));
			}

			if (FeedbackReheatingDeansityScaling == 1) {
				double SigmaGas;
#ifdef H2_AND_RINGS
				SigmaGas = ColdGas / (M_PI * (Radius_high * Radius_high - Radius_low * Radius_low)) / WARM_PHASE_FACTOR
						   * Clumpingfactor;
#else // H2_AND_RINGS
				SigmaGas = ColdGas / (M_PI * (Radius_high * Radius_high - Radius_low * Radius_low));
#endif // H2_AND_RINGS
				SigmaGas = SigmaGas * 0.01 * Hubble_h; // convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2)

				if (SigmaGas > 0.) {
					reheated_mass /= SigmaGas;
				}
			}

			if (reheated_mass * Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir >
			    stars * (EtaSNcode * ReScaled_EnergySNcode)) {
				reheated_mass = stars * (EtaSNcode * ReScaled_EnergySNcode) /
				                (Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir);
			}
		}
	} else {
		reheated_mass = 0.;
	}

	if (reheated_mass > ColdGas) {
		reheated_mass = ColdGas;
	}
	return reheated_mass;
}
// end compute_SN_reheat



















/** @brief 																											  */
double EAGLE2015_rescale_of_EnergySN (double ColdGas, double MetalsColdGas, double Radius_low, double Radius_high) {
	double fth_min = 0.3, fth_max = 3.0, nH_birth, nH_grams = 1.6737236e-24, nH_0 = 0.67, fth, metallicity_Z;
	double z_solar = 0.0127, n = 2 * log (10);

	metallicity_Z = MetalsColdGas / ColdGas;
	nH_birth      = ColdGas * UnitMass_in_g / nH_grams /
	                (4. / 3. * M_PI * UnitLength_in_cm * (Radius_high * Radius_high - Radius_low * Radius_low));

	fth = fth_min +
	      (fth_max - fth_min) / (1 + pow (metallicity_Z / (0.1 * z_solar), n) * pow (nH_birth / nH_0, - 1. * n));
	return fth;
}
// end EAGLE2015_rescale_of_EnergySN



















/** @brief Updates cold, hot and external gas components due to SN reheating and ejection. 							  */
#ifndef H2_AND_RINGS

void update_from_feedback (int p, int centralgal, double reheated_mass, double ejected_mass)
#else // H2_AND_RINGS
void update_from_feedback( int p, int centralgal, double reheated_mass, double ejected_mass,double reheated_massr[] )
#endif // H2_AND_RINGS
{
	int    merger_centre = 0;
	double MassRemain    = 0., dis = 0., fraction;
#ifdef H2_AND_RINGS
	double fractionRings[RNUM], tmpfractionRings[RNUM], MassRemainRings[RNUM];
	int jj;
#endif // H2_AND_RINGS

	if (Gal[p].ColdGas > 0.) {
		// REHEAT if galaxy is a type 1 or a type 2 orbiting a type 1 with hot gas being stripped, some of the reheated
		// and ejected masses goes to the type 0 and some stays in the type 1
#ifdef H2_AND_RINGS
		for ( jj = 0; jj < RNUM; jj++ ) {
			if ( Gal[p].ColdGasRings[jj] > 0. ) {
				fractionRings[jj] = reheated_massr[jj] / (Gal[p].ColdGasRings[jj]);
			}
			else {
				fractionRings[jj] = 0.;
			}
		}
#endif // H2_AND_RINGS

		if (Gal[p].Type == 0) {
#ifdef H2_AND_RINGS
			transfer_material_with_rings(p, "HotGas", p, "ColdGas", fractionRings, "model_starformation_and_feedback.c",
										 __LINE__);
#else // H2_AND_RINGS
			transfer_material (p, "HotGas", p, "ColdGas", ((float) reheated_mass) / ((float) Gal[p].ColdGas),
			                   "model_starformation_and_feedback.c", __LINE__);
#endif // H2_AND_RINGS
		}

			// For satellite galaxies compute how much gas stays in the galaxy and how much goes to central companion
		else if (Gal[p].Type < 3) {
			if (Gal[p].Type == 1) {
				merger_centre = centralgal;
			} else if (Gal[p].Type == 2) {
				merger_centre = Gal[p].CentralGal;
			}

			// If no hot gas in type 2's, share gas between  0 and 1.
			if (HotGasOnType2Galaxies == 0) {
				dis = separation_gal (centralgal, Gal[p].CentralGal) / (1 + ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);
			}
				// If hot gas in type 2's, share gas between itself and merger centre
			else if (HotGasOnType2Galaxies == 1) {
				dis = separation_gal (merger_centre, p) / (1 + ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);
			}

			// Compute share of reheated mass
			if ((dis < Gal[centralgal].Rvir && Gal[Gal[p].CentralGal].Type == 1 && HotGasOnType2Galaxies == 0) ||
			    (dis < Gal[merger_centre].Rvir && HotGasOnType2Galaxies == 1)) {
				//mass that remains on type1 (the rest goes to type 0) for reheat - MassRemain, for eject - ejected_mass
				MassRemain   = reheated_mass * Gal[p].HotRadius / Gal[p].Rvir;
				ejected_mass = ejected_mass * Gal[p].HotRadius / Gal[p].Rvir;
				if (MassRemain > reheated_mass) {
					MassRemain = reheated_mass;
				}
#ifdef H2_AND_RINGS
				for ( jj     = 0; jj < RNUM; jj++ ) {
					MassRemainRings[jj] = reheated_massr[jj] * Gal[p].HotRadius / Gal[p].Rvir;
					if ( MassRemainRings[jj] > reheated_massr[jj] ) {
						MassRemainRings[jj] = reheated_massr[jj];
					}
				}
#endif // H2_AND_RINGS
			} else {
				MassRemain = reheated_mass;
#ifdef H2_AND_RINGS
				for ( jj   = 0; jj < RNUM; jj++ ) {
					MassRemainRings[jj] = reheated_massr[jj];
				}
#endif // H2_AND_RINGS
			}

			// Needed due to precision issues, since we first remove MassRemain and then (reheated_mass-MassRemain)
			// from the satellite into the type 0 and type 1 the fraction might not add up on the second call since
			// Gal[p].ColdGas is a float and reheated_mass & MassRemain are doubles
			if ((MassRemain + reheated_mass) > Gal[p].ColdGas) {
				MassRemain = Gal[p].ColdGas - reheated_mass;
			}

#ifdef H2_AND_RINGS
			for ( jj = 0; jj < RNUM; jj++ ) {
				if ((MassRemainRings[jj] + reheated_massr[jj]) > Gal[p].ColdGasRings[jj] ) {
					MassRemainRings[jj] = Gal[p].ColdGasRings[jj] - reheated_massr[jj];
				}
			}
#endif // H2_AND_RINGS

			mass_checks (p, "model_starformation_and_feedback.c", __LINE__);

			// Transfer MassRemain
			if (reheated_mass > 0.) {
#ifdef H2_AND_RINGS
				for ( jj = 0; jj < RNUM; jj++ ) {
					if ( Gal[p].ColdGasRings[jj] > 0. ) {
						tmpfractionRings[jj] = MassRemainRings[jj] / (Gal[p].ColdGasRings[jj]);
					}
					else {
						tmpfractionRings[jj] = 0.;
					}
				}
				if ( HotGasOnType2Galaxies == 0 ) {
					// Transfer to itself if type 1, merger centre if type 2
					if ( Gal[p].CentralGal == p ) {
						transfer_material_with_rings(Gal[p].CentralGal, "HotGas", p, "ColdGas", tmpfractionRings,
													 "model_starformation_and_feedback.c", __LINE__);
					}
					else {
						transfer_material_with_rings(Gal[p].CentralGal, "HotGas", p, "ColdGas", tmpfractionRings,
													 "model_starformation_and_feedback.c", __LINE__);
					}
				}
				else if ( HotGasOnType2Galaxies == 1 ) {
					// Transfer to itself.
					transfer_material_with_rings(p, "HotGas", p, "ColdGas", tmpfractionRings,
												 "model_starformation_and_feedback.c", __LINE__);
				}
#else // H2_AND_RINGS
				if (HotGasOnType2Galaxies == 0) {
					// Transfer to itself if type 1, merger centre if type 2
					if (Gal[p].CentralGal == p) {
						transfer_material (Gal[p].CentralGal, "HotGas", p, "ColdGas", MassRemain / Gal[p].ColdGas,
						                   "model_starformation_and_feedback.c", __LINE__);
					} else {
						transfer_material (Gal[p].CentralGal, "HotGas", p, "ColdGas", MassRemain / Gal[p].ColdGas,
						                   "model_starformation_and_feedback.c", __LINE__);
					}
				} else if (HotGasOnType2Galaxies == 1) {
					// Transfer to itself
					transfer_material (p, "HotGas", p, "ColdGas", MassRemain / Gal[p].ColdGas,
					                   "model_starformation_and_feedback.c", __LINE__);
				}
#endif
			}

			mass_checks (p, "model_starformation_and_feedback.c", __LINE__);

			//transfer reheated_mass-MassRemain from galaxy to the type 0
			if (reheated_mass > MassRemain) {
				if (Gal[p].ColdGas > 0.) {
					//if the reheat to itself left cold gas below limit do not reheat to central.
#ifdef H2_AND_RINGS
					// Cannot use tmpfractionRings defined from fractionRings since Gal[p].ColdGasRings has changed
					// from MassRemain above
					for ( jj = 0; jj < RNUM; jj++ ) {
						if ( Gal[p].ColdGasRings[jj] > 0. ) {
							fractionRings[jj] = (reheated_massr[jj] - MassRemainRings[jj]) / Gal[p].ColdGasRings[jj];
						}
						else {
							fractionRings[jj] = 0.;
						}
					}
					// Transfer to type 0
					if ( HotGasOnType2Galaxies == 0 ) {
						transfer_material_with_rings(centralgal, "HotGas", p, "ColdGas", fractionRings,
													 "model_starformation_and_feedback.c", __LINE__);
					}
						// Transfer to merger centre
					else if ( HotGasOnType2Galaxies == 1 ) {
						transfer_material_with_rings(merger_centre, "HotGas", p, "ColdGas", fractionRings,
													 "model_starformation_and_feedback.c", __LINE__);
					}
#else // H2_AND_RINGS
					// With rings ColdGas is a double and using (float) might cause
					// (float)(reheated_mass-MassRemain)/Gal[p].ColdGas to be >1
					// Transfer to type 0
					if (HotGasOnType2Galaxies == 0) {
						transfer_material (centralgal, "HotGas", p, "ColdGas",
						                   (float) (reheated_mass - MassRemain) / Gal[p].ColdGas,
						                   "model_starformation_and_feedback.c", __LINE__);
					}
						// Transfer to merger centre
					else if (HotGasOnType2Galaxies == 1) {
						transfer_material (merger_centre, "HotGas", p, "ColdGas",
						                   (float) (reheated_mass - MassRemain) / Gal[p].ColdGas,
						                   "model_starformation_and_feedback.c", __LINE__);
					}
#endif //H2_AND_RINGS
				}
			}

			mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
		} // else if ( Gal[p].Type < 3 )

	} // if ( Gal[p].ColdGas > 0. )

	// Do ejection of gas
	if ((Gal[Gal[p].CentralGal].HotGas > 0. && HotGasOnType2Galaxies == 0) ||
	    (Gal[p].HotGas > 0. && HotGasOnType2Galaxies == 1)) {
		if (HotGasOnType2Galaxies == 0) {
			if (ejected_mass > Gal[Gal[p].CentralGal].HotGas)
				// Either eject own gas or merger_centre gas for ttype 2's
			{
				ejected_mass = Gal[Gal[p].CentralGal].HotGas;
			}
			fraction = ((float) ejected_mass) / Gal[Gal[p].CentralGal].HotGas;
		} else if (HotGasOnType2Galaxies == 1) {
			if (ejected_mass > Gal[p].HotGas && HotGasOnType2Galaxies == 1) {
				// Always eject own gas
				ejected_mass = Gal[p].HotGas;
			}
			fraction = ((float) ejected_mass) / Gal[p].HotGas;
		}

		if (Gal[Gal[p].CentralGal].Type == 1) {
			// If type 1, or type 2 orbiting type 1 near type 0
			if (FateOfSatellitesGas == 0) {
				if (HotGasOnType2Galaxies == 0) {
					transfer_material (Gal[p].CentralGal, "EjectedMass", Gal[p].CentralGal, "HotGas", fraction,
					                   "model_starformation_and_feedback.c", __LINE__);
				} else if (HotGasOnType2Galaxies == 1) {
					transfer_material (Gal[p].CentralGal, "EjectedMass", p, "HotGas", fraction,
					                   "model_starformation_and_feedback.c", __LINE__);
				}
			} else if (FateOfSatellitesGas == 1) {
				if (dis < Gal[centralgal].Rvir) {
					transfer_material (centralgal, "HotGas", Gal[p].CentralGal, "HotGas", fraction,
					                   "model_starformation_and_feedback.c", __LINE__);
				} else {
					transfer_material (Gal[p].CentralGal, "EjectedMass", Gal[p].CentralGal, "HotGas", fraction,
					                   "model_starformation_and_feedback.c", __LINE__);
				}
			}
		} else {
			// If galaxy type 0 or type 2 merging into type 0
			if (HotGasOnType2Galaxies == 0) {
				transfer_material (centralgal, "EjectedMass", Gal[p].CentralGal, "HotGas", fraction,
				                   "model_starformation_and_feedback.c", __LINE__);
			} else if (HotGasOnType2Galaxies == 1) {
				transfer_material (centralgal, "EjectedMass", p, "HotGas", fraction,
				                   "model_starformation_and_feedback.c", __LINE__);
			}
		}

	} // (Gal[Gal[p].CentralGal].HotGas > 0.)

#ifdef H2_AND_RINGS
	update_h2fraction(p);
#endif // H2_AND_RINGS

	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
} // end update_from_feedback



















/** @brief 																											  */
void update_massweightage (int p, double stars, double time) {
	int    outputbin;
	double age; // age in Mpc/Km/s/h - code units

	for (outputbin = 0; outputbin < NOUT; outputbin ++) {
		age = time - NumToTime (ListOutputSnaps[outputbin]);
#ifdef DETAILED_METALS_AND_MASS_RETURN
		Gal[p].MassWeightAge[outputbin] += age * stars;
#else
		Gal[p].MassWeightAge[outputbin] += age * stars * (1. - RecycleFraction);
#endif
	}
} // end update_massweightage



















/* @brief Calculates the stability of the stellar disk as discussed in Mo, Mao & White (1998). For unstable stars, the
 * required amount is transfered to the bulge to make the disk stable again. Mass, metals and luminosities updated.
 * After Guo2010 the bulge size is followed and needs to be updated. Eq 34 & 35 in Guo2010 are used. 				  */

void check_disk_instability (int p, double dt, double time) {

#ifdef DI_INSTABILITIES
	double Mcrit, UnstableMass, GasMass, GasRd, DiskMass, DiskRd, TotalMass, InverseEpsilonStars, InverseEpsilonGas, InverseEpsilonTotal;
	double fraction, EpsilonGas, EpsilonStars, EpsilonMin, EpsilonTotal, Tester;
#ifdef H2_AND_RINGS
	double mstarsRings[RNUM];
	int j;
#endif // H2_AND_RINGS
#else // DI_INSTABILITIES
	double Mcrit, fraction, stars, diskmass;
#endif // DI_INSTABILITIES

#ifdef H2_AND_RINGS
	double rstar, vmax;
	int j;
#endif // H2_AND_RINGS

#ifdef DI_INSTABILITIES
	GasMass   = Gal[p].ColdGas;
	DiskMass  = Gal[p].DiskMass;
	TotalMass = DiskMass + GasMass;
	DiskRd    = Gal[p].DiskRadius / 3.0;
	GasRd     = Gal[p].ColdGasRadius / 3.0;


	// if the mass of one component is zero then the other should be stable on it's own (i.e. Qcomp = EpsilonTotal = 1)
	if (GasMass > 0.0 && DiskMass > 0.0)
	{
		Tester = 2.0;

		// Calculate epsilon parameters for gaseous and stellar disks.
		if (Gal[p].Type != 0)
		{
			EpsilonGas   = Gal[p].InfallVmax / sqrt (G * GasMass / GasRd);
			EpsilonStars = Gal[p].InfallVmax / sqrt (G * DiskMass / DiskRd);
		}
		else
		{
			EpsilonGas   = Gal[p].Vmax / sqrt (G * GasMass / GasRd);
			EpsilonStars = Gal[p].Vmax / sqrt (G * DiskMass / DiskRd);
		}
	}
	else
	{
		Tester = 1.0;
	}

	//	Calculate the total epsilon parameter.
	if (EpsilonStars >= EpsilonGas && Tester == 2.0)
	{
		EpsilonMin          = (TotalMass + DiskMass) / TotalMass;
		InverseEpsilonTotal = DiskMass / (TotalMass * EpsilonStars) + 1.0 / EpsilonGas;
		EpsilonTotal        = 1.0 / InverseEpsilonTotal;

		// Check the stability of the two-component disk.
		if (EpsilonTotal < 1.0)
		{
			if (EpsilonStars >= EpsilonMin)
			{
				InverseEpsilonGas = (TotalMass * EpsilonStars - DiskMass) / (TotalMass * EpsilonStars);

				// Calculate disk's critical mass.
				if (Gal[p].Type != 0)
				{
					Mcrit = InverseEpsilonGas * InverseEpsilonGas * (Gal[p].InfallVmax * Gal[p].InfallVmax * GasRd / G);
				}
				else
				{
					Mcrit = InverseEpsilonGas * InverseEpsilonGas * (Gal[p].Vmax * Gal[p].Vmax * GasRd / G);
				}
				// Calculate the unstable gaseous mass.
				UnstableMass      = GasMass - Mcrit;

				// Invoke a function to deal with the unstable gas.
				if (UnstableMass > 0.0)
				{
					fix_gaseous_disk_instability (p, GasMass, UnstableMass, dt, time);
				}
			} // if (EpsilonStars > EpsilonMin)

			else if (EpsilonStars < EpsilonMin)
			{
				InverseEpsilonGas = TotalMass / (TotalMass + DiskMass);

				// Calculate disk's critical mass.
				if (Gal[p].Type != 0)
				{
					Mcrit = InverseEpsilonGas * InverseEpsilonGas * (Gal[p].InfallVmax * Gal[p].InfallVmax * GasRd / G);
				}
				else
				{
					Mcrit = InverseEpsilonGas * InverseEpsilonGas * (Gal[p].Vmax * Gal[p].Vmax * GasRd / G);
				}
				// Calculate the unstable gaseous mass.
				UnstableMass      = GasMass - Mcrit;

				// Invoke a function to deal with the unstable gas.
				if (UnstableMass > 0.0)
				{
					fix_gaseous_disk_instability (p, GasMass, UnstableMass, dt, time);
				}

				GasMass   = Gal[p].ColdGas;
				DiskMass  = Gal[p].DiskMass;
				TotalMass = DiskMass + GasMass;
				DiskRd    = Gal[p].DiskRadius / 3.0;

				InverseEpsilonStars = (1.0 - InverseEpsilonGas) * TotalMass / DiskMass; // New masses after the formation of stars from unstable gas

				// Calculate disk's critical mass.
				if (Gal[p].Type != 0)
				{
					Mcrit = InverseEpsilonStars * InverseEpsilonStars * (Gal[p].InfallVmax * Gal[p].InfallVmax * DiskRd / G);
				}
				else
				{
					Mcrit = InverseEpsilonStars * InverseEpsilonStars * (Gal[p].Vmax * Gal[p].Vmax * DiskRd / G);
				}
				// Calculate the unstable gaseous mass.
				UnstableMass        = DiskMass - Mcrit;
				fraction            = UnstableMass / DiskMass;

				// Invoke a function to deal with the unstable gas.
				if (UnstableMass > 0.0)
				{
					fix_stellar_disk_instability (p, UnstableMass, fraction);
				}
			} // else if (EpsilonStars < 1.0 + EpsilonMin)
		} // if (EpsilonTotal < 1.0 && Tester == 2.0)
	} // if (EpsilonStars >= EpsilonGas)
	else if (EpsilonStars < EpsilonGas && Tester == 2.0)
	{
		EpsilonMin          = (TotalMass + GasMass) / TotalMass;
		InverseEpsilonTotal = 1.0 / EpsilonStars + GasMass / (TotalMass * EpsilonGas);
		EpsilonTotal        = 1.0 / InverseEpsilonTotal;

		// Check the stability of the two-component disk.
		if (EpsilonTotal < 1.0)
		{
			if (EpsilonGas >= 1.0 + EpsilonMin)
			{
				InverseEpsilonStars = (TotalMass * EpsilonGas - GasMass) / (TotalMass * EpsilonGas);

				// Calculate disk's critical mass.
				if (Gal[p].Type != 0)
				{
					Mcrit = InverseEpsilonStars * InverseEpsilonStars * (Gal[p].InfallVmax * Gal[p].InfallVmax * DiskRd / G);
				}
				else
				{
					Mcrit = InverseEpsilonStars * InverseEpsilonStars * (Gal[p].Vmax * Gal[p].Vmax * DiskRd / G);
				}
				// Calculate the unstable gaseous mass.
				UnstableMass        = DiskMass - Mcrit;
				fraction            = UnstableMass / DiskMass;

				// Invoke a function to deal with the unstable gas.
				if (UnstableMass > 0.0)
				{
					fix_stellar_disk_instability (p, UnstableMass, fraction);
				}
			} // if (EpsilonGas > 1.0 + EpsilonMin)

			else if (EpsilonGas < 1.0 + EpsilonMin)
			{
				InverseEpsilonStars = TotalMass / (TotalMass + GasMass);

				// Calculate disk's critical mass.
				if (Gal[p].Type != 0)
				{
					Mcrit = InverseEpsilonStars * InverseEpsilonStars * (Gal[p].InfallVmax * Gal[p].InfallVmax * DiskRd / G);
				}
				else
				{
					Mcrit = InverseEpsilonStars * InverseEpsilonStars * (Gal[p].Vmax * Gal[p].Vmax * DiskRd / G);
				}
				// Calculate the unstable gaseous mass.
				UnstableMass        = DiskRd - Mcrit;
				fraction            = UnstableMass / DiskMass;

				// Invoke a function to deal with the unstable gas.
				if (UnstableMass > 0.0)
				{
					fix_stellar_disk_instability (p, UnstableMass, fraction);
				}

				GasMass   = Gal[p].ColdGas;
				DiskMass  = Gal[p].DiskMass;
				TotalMass = DiskMass + GasMass;
				GasRd     = Gal[p].ColdGasRadius / 3.0;

				InverseEpsilonGas = (1.0 - InverseEpsilonStars) * TotalMass / GasMass; // New masses after the formation of stars from unstable gas

				// Calculate disk's critical mass.
				if (Gal[p].Type != 0)
				{
					Mcrit = InverseEpsilonGas * InverseEpsilonGas * (Gal[p].InfallVmax * Gal[p].InfallVmax * GasRd / G);
				}
				else
				{
					Mcrit = InverseEpsilonGas * InverseEpsilonGas * (Gal[p].Vmax * Gal[p].Vmax * GasRd / G);
				}
				// Calculate the unstable gaseous mass.
				UnstableMass      = GasMass - Mcrit;

				// Invoke a function to deal with the unstable gas.
				if (UnstableMass > 0.0)
				{
					fix_gaseous_disk_instability (p, GasMass, UnstableMass, dt, time);
				}
			} // else if (EpsilonStars < 1.0 + EpsilonMin)
		} // if (EpsilonTotal < 1.0)
	} // else if (EpsilonStars < EpsilonGas && Tester == 2.0)
	else if (Tester == 1.0)
	{
		if (GasMass <= 0.0)
		{
			if (Gal[p].Type != 0)
			{
				Mcrit = (Gal[p].InfallVmax * Gal[p].InfallVmax * (Gal[p].DiskRadius / 3.0)) / G;
			}
			else
			{
				Mcrit = (Gal[p].Vmax * Gal[p].Vmax * (Gal[p].DiskRadius / 3.0)) / G;
			}
			// Calculate the unstable stellar mass.
			UnstableMass = Gal[p].DiskMass - Mcrit;
			fraction     = UnstableMass / Gal[p].DiskMass;

			if (UnstableMass > 0.0)
			{
				fix_stellar_disk_instability (p, UnstableMass, fraction);
			}
		}
		else if (DiskMass <= 0.0)
		{
			// Calculate disk's critical mass.
			if (Gal[p].Type != 0)
			{
				Mcrit = (Gal[p].InfallVmax * Gal[p].InfallVmax * (Gal[p].ColdGasRadius / 3.0)) / G;
			}
			else
			{
				Mcrit = (Gal[p].Vmax * Gal[p].Vmax * (Gal[p].ColdGasRadius / 3.0)) / G;
			}
			// Calculate the unstable gaseous mass.
			UnstableMass = GasMass - Mcrit;

			// Invoke a function to deal with the unstable gas.
			if (UnstableMass > 0.0)
			{
				fix_gaseous_disk_instability (p, GasMass, UnstableMass, dt, time);
			}
		}
	} // else if (EpsilonTotal < 1.0 && EpsilonMin == 1.0)
#else // DI_INSTABILITIES

	diskmass = Gal[p].DiskMass;
#ifndef H2_AND_RINGS
// Check stellar disk -> eq 34 Guo2010
	if (Gal[p].Type != 0) {
		Mcrit = Gal[p].InfallVmax * Gal[p].InfallVmax * Gal[p].DiskRadius / G;
	} else {
		Mcrit = Gal[p].Vmax * Gal[p].Vmax * Gal[p].DiskRadius / G;
	}
#else // H2_AND_RINGS
	if ( diskmass < 1.0e-6 ) {
		rstar = 0.5 * RingRadius[0];
	}
	else {
		rstar   = 0.5 * RingRadius[0] * Gal[p].DiskMassRings[0];
		for ( j = 1; j < RNUM; j++ ) {
			rstar += 0.5 * (RingRadius[j - 1] + RingRadius[j]) * Gal[p].DiskMassRings[j];
		}
		rstar   = rstar / diskmass / 2.0;      //2.0=mean radius/scale length for exponential disk
	}

	if ( Gal[p].Type != 0 ) {
		vmax = Gal[p].InfallVmax;
	}
	else {
		vmax = Gal[p].Vmax;
	}

	Mcrit = vmax * vmax * rstar / G;
#endif // H2_AND_RINGS
	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);

	stars    = diskmass - Mcrit;
	fraction = stars / diskmass;

	// Add excess stars to the bulge.
	if (stars > 0.0) {
		// Calculate the bulge size and update the disk size.
		update_bulgesize_from_stellardisk_instability (p, stars);
		//The bulge will be formed in the same place as the disk was, so the disk rings are transferred directly into bulge rings.
#ifdef H2_AND_RINGS
		double dmass, fractionRings[RNUM];
		for ( j = 0; j < RNUM; j++ ) {
			fractionRings[j] = 0.;
		}

		dmass = stars;
		j     = 0; //avoid non-definded j if dmass<1e-6
		if ( dmass > 1.0e-6 ) {
			for ( j = 0; j < RNUM; j++ ) {
				// Mass is transferred first from the inner rings until the necessary mass is achieved.
				if ( dmass > Gal[p].DiskMassRings[j] ) {
					dmass -= Gal[p].DiskMassRings[j];
					fractionRings[j] = 1.;
				}
				else
					break;
			}
		}

		// Check needed in case there is a ring with 0 mass in the middle
		if ( Gal[p].DiskMassRings[j] > 0 ) {
			fractionRings[j] = dmass / Gal[p].DiskMassRings[j];
		}
		else {
			fractionRings[j] = 0.;
		}
		transfer_material_with_rings(p, "BulgeMass", p, "DiskMass", fractionRings, "model_starformation_and_feedback.c",
									 __LINE__);
#else // H2_AND_RINGS

// Transfer the unstable stars directly to the bulge. The disk has to regain dynamical equilibrium by
// transporting mass to the bulge. We store the mass transferred (stars) from the disk to the bulge which
// represents the instability-driven bulge mass.
		Gal[p].PBMass += stars;
		transfer_material (p, "BulgeMass", p, "DiskMass", fraction, "model_starformation_and_feedback.c", __LINE__);
#endif // H2_AND_RINGS

		if (BHGrowthInDiskInstabilityModel == 1) {
			if (Gal[p].ColdGas > 0.) {
#ifdef H2_AND_RINGS
				for (j = 0; j < RNUM; j ++)
				{
					fractionRings[j] *=
							0.1 * BlackHoleGrowthRate / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir)));
					fractionRings[j] = min(1.0, fractionRings[j]);
				}

				transfer_material_with_rings (p, "BlackHoleMass", p, "ColdGas", fractionRings,
						"model_starformation_and_feedback.c", __LINE__);
#else // H2_AND_RINGS
				fraction *= 0.1 * BlackHoleGrowthRate / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir)));
				transfer_material (p, "BlackHoleMass", p, "ColdGas", fraction, "model_starformation_and_feedback.c", __LINE__);
#endif // H2_AND_RINGS
				Gal[p].QuasarAccretionRate += fraction * Gal[p].ColdGas / (dt * STEPS);
			}
		}
#ifdef BULGESIZE_DEBUG
		mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

		if ((Gal[p].BulgeMass > TINY_MASS && Gal[p].BulgeSize < TINY_LENGTH)||
		(Gal[p].BulgeMass < TINY_MASS && Gal[p].BulgeSize > TINY_LENGTH)) {
		printf("BulgeMass=%g BulgeSize=%g\n",Gal[p].BulgeMass,Gal[p].BulgeSize);
		terminate("bulgesize wrong in disk instablility\n");
		}
#endif // BULGESIZE_DEBUG
		mass_checks (p, "model_starformation_and_feedback.c", __LINE__);

		if ((Gal[p].BulgeMass > 1e-9 && Gal[p].BulgeSize == 0.0) ||
		    (Gal[p].BulgeMass == 0.0 && Gal[p].BulgeSize > 1e-9)) {
			char sbuf[1000];
			sprintf (sbuf,
			         "bulgesize wrong in disk instablility.c \n");
			printf ("BulgeMass=%g BulgeSize=%g\n", Gal[p].BulgeMass, Gal[p].BulgeSize);
			terminate(sbuf);
		}
	}// if(stars > 0.0)
#endif // DI_INSTABILITIES
} //end check_disk_instability


















/* Calculate the bulge size and update the disk size. Transfer the unstable stars directly to the bulge.*/
void fix_stellar_disk_instability (int p, double UnstableMass, double fraction) {
	// The bulge will be formed in the same place as the disk was.
	update_bulgesize_from_stellardisk_instability (p, UnstableMass);
	// We store the mass transferred (UnstableMass) from the disk to the bulge which represents the pseudo-bulge mass.
	Gal[p].PBMass += UnstableMass;
	// The disk has to regain its dynamical equilibrium by transporting mass to the bulge.
	transfer_material (p, "BulgeMass", p, "DiskMass", fraction, "model_starformation_and_feedback.c", __LINE__);
} // fix_stellar_disk_instability

















/* Gaseous disk instabilities can lead to black hole growth and star formation.*/
void fix_gaseous_disk_instability (int p, double GasMass, double UnstableMass, double dt, double time) {
	double BHGrowth, BHfraction, starburst;
	int    j;
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
	double metallicitySF;
#endif // POST_PROCESS_MAGS
#endif // COMPUTE_SPECPHOT_PROPERTIES

	// A fraction of the unstable gas is accreted into the black hole. We use here the same formula as in Kauffmann & Haehnelt (2000)(Quasar Mode).
	if (Gal[p].BlackHoleMass > 0.0) {
		BHGrowth   = BlackHoleGrowthRate * UnstableMass / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir)));
		BHfraction = BHGrowth / GasMass;

		transfer_material (p, "BlackHoleMass", p, "ColdGas", BHfraction, "model_starformation_and_feedback.c", __LINE__);
		// Update gas disk spin and radius.
		for (j = 0; j < 3; j ++) {
			Gal[p].ColdGasSpin[j] = Gal[p].ColdGasSpin[j] / (1 - BHfraction);
		}
		if (DiskRadiusModel == 0) {
			Gal[p].ColdGasRadius = get_gas_disk_radius (p);
		}
		// Update the quasar accretion rate.
		Gal[p].QuasarAccretionRate += BHGrowth / (dt * STEPS);

		// The remaining gas is converted into stars that are left in the stellar disk.
		starburst = UnstableMass - BHGrowth;
	} else {
		// All unstable gas is converted into stars that are left in the stellar disk.
		starburst = UnstableMass;
	}
	if (starburst > 0.0) {
		// Deal with the disk instability-induced starburst feedback.
		// Store the value of the metallicity of the cold phase when SF occurs. Used to update luminosities below.
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
		metallicitySF = metals_total(Gal[p].MetalsColdGas) / Gal[p].ColdGas;
#endif // POST_PROCESS_MAGS
#endif // COMPUTE_SPECPHOT_PROPERTIES
		// If FEEDBACK_COUPLED_WITH_MASS_RETURN feedback happens only when stars die, no need to balance it with SF.
#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
		update_stars_due_to_reheat (p, Gal[p].CentralGal, &starburst);
#endif //FEEDBACK_COUPLED_WITH_MASS_RETURN
		// Update the star formation rate.
		Gal[p].Sfr += starburst / (dt * STEPS);
		mass_checks (p, "model_starformation_and_feedback.c", __LINE__);

		// Update_from_star_formation can only be called after SN_feeedback recipe since stars need to be re_set once
		// the reheated mass is known (star formation and feedback share the same fraction of cold gas).
		int nstep = - 1;
		update_from_star_formation (p, starburst, "model_starformation_and_feedback", nstep);
		mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
		update_massweightage (p, starburst, time);

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
		SN_feedback (p, Gal[p].CentralGal, starburst, "ColdGas");
#endif // FEEDBACK_COUPLED_WITH_MASS_RETURN

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
		// Update the luminosities due to the stars formed.
		add_to_luminosities (p, starburst, time, (dt * STEPS) / STEPS, metallicitySF);
#endif // POST_PROCESS_MAGS
#endif // COMPUTE_SPECPHOT_PROPERTIES
	} // if (starburst > 0.0)
} // fix_gaseous_disk_instability



















/* @brief Updates bulge from disk instability -> stars represents the mass transferred to the bulge, which occupies a
 * size in the bulge equal to that occupied in the disk. 															  */
void update_bulgesize_from_stellardisk_instability (int p, double stars) {
	double bulgemass, bulgesize, diskmass, fint, massfrac;
	int    j;

	// alpha_inter=2.0/C=0.5 (alpha larger than in mergers since the existing and newly formed bulges are concentric).
	fint = 4.0;

	// Disk instabilities are assumed to NOT remove angular momentum. Since the mass of the disk changes, the spin is
	// changed by the same amount to keep angular momentum constant.

#ifdef DI_BULGESPIN
	bulgemass = Gal[p].BulgeMass;
	massfrac  = stars / bulgemass;

	for ( j   = 0; j < 3; j ++ ) {
		Gal[p].BulgeSpin[j] = Gal[p].BulgeSpin[j] / (1 + massfrac);
	}
#endif // DI_BULGESPIN

	diskmass = Gal[p].DiskMass;
	massfrac = stars / diskmass;

	for (j = 0; j < 3; j ++) {
		// if everything is transferred to the bulge
		if (massfrac == 1) {
			Gal[p].DiskSpin[j] = 0;
		} else {
			Gal[p].DiskSpin[j] = Gal[p].DiskSpin[j] / (1 - massfrac);
		}
	}
	if (DiskRadiusModel == 0) {
		Gal[p].DiskRadius = get_stellar_disk_radius (p);
	}

#ifndef H2_AND_RINGS
#ifdef BULGESIZE_DEBUG
	double orisize;
	orisize = Gal[p].BulgeSize;
#endif // BULGESIZE_DEBUG

	// Size of newly formed bulge, which consists of the stellar mass transferred from the disk. This is calculated
	// using bulge_from_disk which receives Delta_M/DiskMass and returns Rb/Rd. From eq 35 and since
	// DiskMass=2PISigma(Rd)^2 we see that Delta_M/DiskMass=1-(1+Rb/Rd)*exp(-Rb/Rd), so function bulge_from_disk
	// avoids calculating the slow "ln" function.
	bulgesize = bulge_from_disk (massfrac) * Gal[p].DiskRadius / 3.;
	if (Gal[p].BulgeMass < TINY_MASS) {
		// If previous Bulge Mass = 0 then the bulge size is given directly from newly formed bulge
		Gal[p].BulgeSize = bulgesize;
	} else {
		// Combine the old with newly formed bulge and calculate the bulge size assuming energy conservation as for
		// mergers but using alpha=2. - eq 33 *
		Gal[p].BulgeSize = (Gal[p].BulgeMass + stars) * (Gal[p].BulgeMass + stars) /
		                   (Gal[p].BulgeMass * Gal[p].BulgeMass / Gal[p].BulgeSize +
		                    stars * stars / bulgesize +
		                    fint * Gal[p].BulgeMass * stars / (Gal[p].BulgeSize + bulgesize));
	}
	// Added by PAT to see if the cause of low bulge sizes could be propagation of tightly-bound bulges. These originate
	// in unfeasibly small disks.
#define BULGESIZE_MIN 1e-4
	//	Gal[p].BulgeSize = max(Gal[p].BulgeSize, BULGESIZE_MIN);

#ifdef BULGESIZE_DEBUG
	if ((Gal[p].BulgeMass + stars > TINY_MASS && Gal[p].BulgeSize < TINY_LENGTH)
		|| (Gal[p].BulgeMass + stars < TINY_MASS && Gal[p].BulgeSize > TINY_LENGTH)) {
		printf("Original DiskMass=%e, DiskSize=%e\nOriginal BulgeMass=%e, BulgeSize=%e\nTransferred stars=%e,"
					   "bulgesize=%e\nFinal BulgeMass=%e, BulgeSize=%e\n", Gal[p].DiskMass, Gal[p].DiskRadius,
			   Gal[p].BulgeMass, orisize, stars, bulgesize, Gal[p].BulgeMass + stars, Gal[p].BulgeSize);
		terminate("bulgesize or mass wrong in disk instablility");
	}
#endif // BULGESIZE_DEBUG
#else //H2_AND_RINGS

	// Size of new formed bulge, which consist of the stellar mass trasfered from the disk combine the old bulge with
	// the new materials and caculate the bulge size assuming energy conservation
	diskmass = stars;
	j        = 0; // avoid non-definded j if dmass<1e-6
	if ( diskmass > 1.0e-6 ) {
		for ( j = 0; j < RNUM; j++ ) {
			// Mass is transfered first from the inner rings first until the necessary mass is achieved.
			if ( diskmass > Gal[p].DiskMassRings[j] ) {
				diskmass -= Gal[p].DiskMassRings[j];
			}
			else
				break;
		}
		if ( j == RNUM ) {
			bulgesize = RingRadius[RNUM - 1];
		}
		else {
			if ( j == 0 ) {
				bulgesize = diskmass / Gal[p].DiskMassRings[j] * RingRadius[j];
			}
			else {
				bulgesize = diskmass / Gal[p].DiskMassRings[j] * RingRadius[j] +
							(1 - diskmass / Gal[p].DiskMassRings[j]) * RingRadius[j - 1];
			}
		}
	}
	else {
		bulgesize = 0.5 * RingRadius[0];
	}

	if ( Gal[p].BulgeMass < 1.e-9 ) {
		Gal[p].BulgeSize = bulgesize;
	}
	else {
		Gal[p].BulgeSize = (Gal[p].BulgeMass + stars) * (Gal[p].BulgeMass + stars) /
						   (Gal[p].BulgeMass * Gal[p].BulgeMass / Gal[p].BulgeSize + stars * stars / bulgesize +
							fint * Gal[p].BulgeMass * stars / (Gal[p].BulgeSize + bulgesize));
	}
#endif //H2_AND_RINGS
}
// end update_bulgesize_from_stellardisk_instability



















/** @brief Calculates the size of the disk that contains the mass transferred to the bulge. The bulge is assumed to
 * form with the same size. avoid doing "ln" from eq 35																  */
double bulge_from_disk (double frac) {
	double x1, x2, x0, value;

	x1    = 0.0;
	x2    = 1.;
	while ((func_size (x2, frac) * func_size (x1, frac)) > 0) {
		x1 = x2;
		x2 = x2 * 2;
	}
	x0    = x1 + (x2 - x1) / 2.;
	value = func_size (x0, frac);
	if (value < 0) {
		value = - value;
	}

	while (value > 0.00001) {
		if (func_size (x0, frac) * func_size (x2, frac) > 0) {
			x2 = x0;
		} else {
			x1 = x0;
		}
		x0    = x1 + (x2 - x1) / 2.;
		value = func_size (x0, frac);
		if (value < 0) {
			value = - value;
		}
	}

	return x0;
}
// end bulge_from_disk



















/** @brief 																											  */
double func_size (double x, double a) {
	return exp (- x) * (1 + x) - (1 - a);
}
// end func_size



















/** @brief 																											  */
#ifdef H2_AND_RINGS
void update_h2fraction( int p ) {
	int    j;
	// The central stellar surface density converted from (10^10M_sun/h)/(Mpc/h)^2 to (M_sun/pc^2)
	double SigmaHRings;

	Gal[p].H2fraction = 0.;
	for ( j = 0; j < RNUM; j++ ) {
		//KMT09 or Krumholz et al. 2008
		if ( H2FractionRecipe == 0 || H2FractionRecipe == 1 ) {
			double metallicityr;

			if ( Gal[p].ColdGasRings[j] < 1.0e-8 ) {
				metallicityr = 0.0;
			}
			else {
				metallicityr = metals_total(Gal[p].MetalsColdGasRings[j]) / (Gal[p].ColdGasRings[j] * 0.0134);
			}

			if ( metallicityr < 0.01 ) {
				metallicityr = 0.01;
			}
			if ( j == 0 ) {
				SigmaHRings = Gal[p].ColdGasRings[j] / ((M_PI * RingRadius[j] * RingRadius[j]) * WARM_PHASE_FACTOR);
			}
			else {
				SigmaHRings = Gal[p].ColdGasRings[j] /
							  ((M_PI * (RingRadius[j] * RingRadius[j] - RingRadius[j - 1] * RingRadius[j - 1])) *
							   WARM_PHASE_FACTOR);
			}

			// Now convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2) -> *(1.e10/h)/(1.e6/h*1.e6/h)
			SigmaHRings = SigmaHRings * 0.01 * Hubble_h;

			// Update to clumping factor, as in Fu2013, to solve problems with low-z galaxies
			if ( metallicityr < 1.0 ) {
				SigmaHRings = SigmaHRings * Clumpingfactor * pow((1.0 / metallicityr), 0.7);
			}

			// KMT09 - updated in Fu2013 (eq 11 and 12)
			if ( H2FractionRecipe == 0 ) {
				double tau, khi, s;
				khi = 3.1 * (1. + 3.1 * pow(metallicityr, 0.365)) / 4.1;
				tau = 0.066 * SigmaHRings * metallicityr;
				s   = log(1. + 0.6 * khi + 0.01 * khi * khi) / (0.6 * tau);

				if ( s < 2.0 ) {
					Gal[p].H2fractionRings[j] = 1 - 0.75 * s / (1 + 0.25 * s);
				}
				else {
					Gal[p].H2fractionRings[j] = 0.0;
				}
			}
				// Krumholz et al. 2008
			else if ( H2FractionRecipe == 1 ) {
				// Convert to log10
				metallicityr = log10(metallicityr);
				SigmaHRings  = log10(SigmaHRings);
				Gal[p].H2fractionRings[j] = update_H2fraction_KMT08(SigmaHRings, metallicityr);
			}
		}
			//Blitz & Rosolowsky 2006, pressure recipe
		else if ( H2FractionRecipe == 2 ) {
			double SigmaStarRings, alpha_p = 0.92;
			double SigmaStar0              =
						   (Gal[p].DiskMassRings[0] / (RingRadius[0] * RingRadius[0] * M_PI)) * 0.01 * Hubble_h;
			if ( j == 0 ) {
				SigmaHRings    = Gal[p].ColdGasRings[j] / (M_PI * RingRadius[j] * RingRadius[j]) / WARM_PHASE_FACTOR;
				SigmaStarRings = Gal[p].DiskMassRings[j] / (M_PI * RingRadius[j] * RingRadius[j]);
			}
			else {
				SigmaHRings    = Gal[p].ColdGasRings[j] /
								 (M_PI * (RingRadius[j] * RingRadius[j] - RingRadius[j - 1] * RingRadius[j - 1])) /
								 WARM_PHASE_FACTOR;
				SigmaStarRings = Gal[p].DiskMassRings[j] /
								 (M_PI * (RingRadius[j] * RingRadius[j] - RingRadius[j - 1] * RingRadius[j - 1]));
			}
			SigmaHRings *= (0.01 * Hubble_h);
			SigmaStarRings *= (0.01 * Hubble_h);    //from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2) */
			Gal[p].H2fractionRings[j] =
					1.38e-3 * pow(SigmaHRings * (SigmaHRings + 0.1 * sqrt(SigmaStar0 * SigmaStarRings)), alpha_p);
			if ( Gal[p].H2fractionRings[j] < 1.0e-8 ) {
				Gal[p].H2fractionRings[j] = 0.0;
			}
			else {
				Gal[p].H2fractionRings[j] = 1 / (1 + 1 / Gal[p].H2fractionRings[j]);
			}
		}
		else {
			Gal[p].H2fractionRings[j] = 0;
		}

		Gal[p].H2fraction += Gal[p].H2fractionRings[j] * Gal[p].ColdGasRings[j] / Gal[p].ColdGas;

	}

	if ( Gal[p].ColdGas < 1.0e-7 ) {
		Gal[p].H2fraction = 0.0;
	}
}
// end update_h2fraction



















/** @brief Only used if H2FractionRecipe=1 																			  */
double update_H2fraction_KMT08( double logsigmah, double metallicity ) {
	int    i, j;
	double logNHtot[LENSIGMAH], lgZ[LENZ], mf, mf1, mf2;
	for ( i = 0, logNHtot[0] = -1; i < (LENSIGMAH - 1); i++ ) {
		logNHtot[i + 1] = logNHtot[i] + 0.05;
	}
	for ( j = 0, lgZ[0] = -2; j < (LENZ - 1); j++ ) {
		lgZ[j + 1] = lgZ[j] + 0.25;
	}

	if ( logsigmah < logNHtot[0] ) {
		logsigmah = logNHtot[0];
	}
	if ( logsigmah > logNHtot[i - 1] ) {
		logsigmah = logNHtot[i - 1];
	}
	for ( i             = 0; logsigmah > logNHtot[i + 1]; i++ );

	if ( metallicity < lgZ[0] ) {
		metallicity = lgZ[0];
	}
	if ( metallicity > lgZ[j - 1] ) {
		metallicity = lgZ[j - 1];
	}
	for ( j             = 0; metallicity > lgZ[j + 1]; j++ );

	mf1 = h2frac[i][j] + (h2frac[i][j + 1] - h2frac[i][j]) * (metallicity - lgZ[j]) / (lgZ[j + 1] - lgZ[j]);
	mf2 = h2frac[i + 1][j] + (h2frac[i + 1][j + 1] - h2frac[i + 1][j]) * (metallicity - lgZ[j]) / (lgZ[j + 1] - lgZ[j]);
	mf  = mf1 + (mf2 - mf1) * (logsigmah - logNHtot[i]) / (logNHtot[i + 1] - logNHtot[i]);

	return (mf);
}
// end update_H2fraction_KMT08
#endif // H2_AND_RINGS
