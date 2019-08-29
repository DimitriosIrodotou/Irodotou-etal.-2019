#include "allvars.h"
#include "proto.h"
#include "hdf5.h"
#include "hdf5_hl.h"

#define NRECORDS (hsize_t) 0
#define HDF5_STRING_SIZE 2048

/* File to set all the HDF5 table properties */
int     ifield;
int     rowsize    = 0;
hsize_t chunk_size = CHUNK_SIZE;
int*fill_data = NULL;
hid_t file_id;
size_t*output_offsets;
hid_t *field_types;
size_t*output_sizes;
size_t output_size;

// Write the datatype for HDF5 to write out the data
char types[] = {
		'i',
		'i',
		'f',
		'f',
		'f',
		'3',
		'3',
		'f',
		'f',
		'3',
		'i',
		'f',
		'f',
		'f',
		'3',
		'3',
		'3',
		'i',
		'f',
		'f',
		'f',
		'f',
		'3',
		'3',
		'3',
		'f',
		'f',
		'i',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'f',
		'i',
		'i',
		'o',
		'o',
		'o',
		'f',
};

const char*field_names[] = {
		"Type",
		"HaloIndex",
		"HaloM_Mean200",
		"HaloM_Crit200",
		"HaloM_TopHat",
		"HaloPos",
		"HaloVel",
		"HaloVelDisp",
		"HaloVmax",
		"HaloSpin",
		"SnapNum",
		"LookBackTimeToSnap",
		"CentralMvir",
		"CentralRvir",
		"DistanceToCentralGal",
		"Pos",
		"Vel",
		"Len",
		"Mvir",
		"Rvir",
		"Vvir",
		"Vmax",
		"ColdGasSpin",
		"DiskSpin",
		"BulgeSpin",
		"InfallVmax",
		"InfallVmaxPeak",
		"InfallSnap",
		"InfallHotGas",
		"HotRadius",
		"OriMergTime",
		"MergTime",
		"NMajorMergers",
		"NMinorMergers",
		"ColdGas",
		"StellarMass",
		"DiskMass",
		"BulgeMass",
		"PBMass",
		"CBMassMajor",
		"CBMassMinor",
		"CBMass",
		"HotGas",
		"ReheatedGas",
		"EjectedMass",
		"BlackHoleMass",
		"ICM",
		"MetalsColdGas",
		"MetalsStellarMass",
		"MetalsDiskMass",
		"MetalsBulgeMass",
		"MetalsHotGas",
		"MetalsEjectedMass",
		"MetalsICM",
		"PrimordialAccretionRate",
		"CoolingRadius",
		"CoolingRate",
		"CoolingRate_beforeAGN",
		"QuasarAccretionRate",
		"RadioAccretionRate",
		"Sfr",
		"SfrBulge",
		"XrayLum",
		"BulgeSize",
		"DiskRadius",
		"ColdGasRadius",
		"StellarHalfMassRadius",
		"StellarHalfLightRadius",
		"StellarNinetyLightRadius",
		"CosInclination",
		"DisruptOn",
		"MergeOn",
		"MagDust",
		"Mag",
		"MagBulge",
		"MassWeightAge",
};

//The number of fields in the data
int nfields = 76;

// Define the dimensions for the HDF5 table
hsize_t float_dims[1] = {3};
#ifdef COMPUTE_SPECPHOT_PROPERTIES
hsize_t    nmag_dims[1]={NMAG};
#endif
#ifdef STAR_FORMATION_HISTORY
hsize_t    sfh_dims[1]={SFH_NBIN};
hsize_t    sfh_3_dims[1]={3*SFH_NBIN};
#ifdef INDIVIDUAL_ELEMENTS
hsize_t    numele_dims[1]={NUM_ELEMENTS};
hsize_t    numele_sfh_dims[2]={NUM_ELEMENTS,SFH_NBIN};
hsize_t    mag_sfh_dims[2]={3,SFH_NBIN};
#endif //INDIVIDUAL_ELEMENTS
#endif //STAR_FORMATION_HISTORY
