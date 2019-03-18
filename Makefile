EXEC  = L-Galaxies

# Default object files (others may be added with -D options)
OBJS  = /Users/Bam/ClionProjects/ITH18/main.o \
	/Users/Bam/ClionProjects/ITH18/io_tree.o \
	/Users/Bam/ClionProjects/ITH18/init.o \
	/Users/Bam/ClionProjects/ITH18/cool_func.o \
	/Users/Bam/ClionProjects/ITH18/save.o \
	/Users/Bam/ClionProjects/ITH18/save_galtree.o \
	/Users/Bam/ClionProjects/ITH18/mymalloc.o \
	/Users/Bam/ClionProjects/ITH18/read_parameters.o \
	/Users/Bam/ClionProjects/ITH18/peano.o \
	/Users/Bam/ClionProjects/ITH18/allvars.o \
	/Users/Bam/ClionProjects/ITH18/age.o \
	/Users/Bam/ClionProjects/ITH18/update_type_two.o \
	/Users/Bam/ClionProjects/ITH18/metals.o \
	/Users/Bam/ClionProjects/ITH18/model_infall.o \
	/Users/Bam/ClionProjects/ITH18/model_cooling.o \
	/Users/Bam/ClionProjects/ITH18/model_starformation_and_feedback.o \
	/Users/Bam/ClionProjects/ITH18/model_reincorporation.o \
	/Users/Bam/ClionProjects/ITH18/model_mergers.o \
	/Users/Bam/ClionProjects/ITH18/model_dust.o \
	/Users/Bam/ClionProjects/ITH18/model_misc.o \
	/Users/Bam/ClionProjects/ITH18/model_disrupt.o \
	/Users/Bam/ClionProjects/ITH18/model_stripping.o \
	/Users/Bam/ClionProjects/ITH18/scale_cosmology.o

# The following is used only to set dependencies
INCL  = Makefile \
    /Users/Bam/ClionProjects/ITH18/allvars.h \
	/Users/Bam/ClionProjects/ITH18/h_funcs.h \
	/Users/Bam/ClionProjects/ITH18/h_params.h \
	/Users/Bam/ClionProjects/ITH18/h_metals.h \
	/Users/Bam/ClionProjects/ITH18/h_galaxy_output.h \
	/Users/Bam/ClionProjects/ITH18/h_galaxy_tree_data.h \
	/Users/Bam/ClionProjects/ITH18/h_galaxy.h \
	/Users/Bam/ClionProjects/ITH18/h_halo_data.h \
	/Users/Bam/ClionProjects/ITH18/h_halo_ids_data.h \
	/Users/Bam/ClionProjects/ITH18/h_halo_aux_data.h \
	/Users/Bam/ClionProjects/ITH18/h_lightcone.h \
	/Users/Bam/ClionProjects/ITH18/h_variables.h \
	/Users/Bam/ClionProjects/ITH18/proto.h
ifeq (ALL_SKY_LIGHTCONE,$(findstring ALL_SKY_LIGHTCONE,$(OPT)))
INCL  += /Users/Bam/ClionProjects/ITH18/lightcone.h
endif

# Either include the default set of Makefile options, or define your own
#include Makefile_options
include My_Makefile_options/My_Makefile_options_Hen15_MR

# Choose your system type (needs to match an entry in Makefile_compilers)
#SYSTYPE = "ETH"
#include Makefile_compilers
# Alternatively, My_Makefile_compilers is an extract from Makefile_compilers
include My_Makefile_compilers

#LIBS   =   -g $(LDFLAGS) -lm  $(GSL_LIBS)  $(RLIBS) -lgsl -lgslcblas $(HDF5_LIBS) -lhdf5_serial -lhdf5_serial_hl
LIBS   =   -g $(LDFLAGS) -lm  $(GSL_LIBS)  $(RLIBS) -lgsl -lgslcblas $(HDF5_LIBS) -lhdf5 -lhdf5_hl

CFLAGS =   -g $(OPTIONS) $(OPT) -DCOMPILETIMESETTINGS=\""$(OPT)"\" $(OPTIMIZE) $(GSL_INCL) $(HDF5_INCL)

all: metadata $(EXEC)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) Makefile Makefile_options My_Makefile_options My_Makefile_options_MCMC My_Makefile_options_MCMC_Halo_Model Makefile_compilers ##My_Makefile_compilers

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) .$(EXEC)

# use next target to generate metadata about the result files
# uses -E compiler option to preprocess the allvars.h file, stores result in allvars.i
# uses -CC compiler option to save comments, needed for HDF5 output
# then calls awk scripts from ./awk/ folder to extract cleand-up version of GALAXY_OUTPUT struct
# and generate different representations of use for post-processing the result 	
metadata:
	${CC_MD} ${OPT} ${CFLAGS} -E -CC /Users/Bam/ClionProjects/ITH18/h_galaxy_output.h -o /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i
	${CC_MD} ${OPT} ${CFLAGS} -E -CC /Users/Bam/ClionProjects/ITH18/h_metals.h -o /Users/Bam/ClionProjects/ITH18/h_metals.i
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_TypeString.awk > ./AuxCode/awk/output/L-Galaxies_Types.txt
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_DDL.awk > ./AuxCode/awk/output/L-Galaxies_DDL.sql
ifeq (NORMALIZEDDB,$(findstring NORMALIZEDDB,$(OPT)))
	awk -f ./AuxCode/awk/extract_SFH_BIN.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/SFH_BIN_2_DDL.awk >> ./AuxCode/awk/output/L-Galaxies_DDL.sql
else
	awk -f ./AuxCode/awk/extract_SFH_Time.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/SFH_Time_2_DDL.awk >> ./AuxCode/awk/output/L-Galaxies_DDL.sql
endif	
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_struct.awk >  ./AuxCode/awk/output/idl/LGalaxy.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_hists.awk > ./AuxCode/awk/output/idl/LGalaxy_plot.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_testfloats.awk > ./AuxCode/awk/output/idl/LGalaxy_testfloats.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_zerofloats.awk > ./AuxCode/awk/output/idl/LGalaxy_zerofloats.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_LGalaxy.awk > ./AuxCode/awk/output/L-Galaxies.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_FileFormat.awk > ./AuxCode/awk/output/L-Galaxies_FileFormat.csv
	awk -f ./AuxCode/awk/extract_SFH_BIN.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/MOMAF_INPUT_2_MoMaFGalaxy.awk >> ./AuxCode/awk/output/L-Galaxies.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_python_struct.awk >  ./AuxCode/awk/output/python/LGalaxy.py
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_HDF5.awk > /Users/Bam/ClionProjects/ITH18/io_hdf5.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT_props.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_prop_2_HDF5_proptable.awk > ./input/hdf5_field_props.txt

	awk -f ./AuxCode/awk/extract_struct_metals.awk /Users/Bam/ClionProjects/ITH18/h_metals.i > ./AuxCode/awk/output/structs.dat
	awk -f ./AuxCode/awk/extract_struct_elements.awk /Users/Bam/ClionProjects/ITH18/h_metals.i >> ./AuxCode/awk/output/structs.dat
	awk -f ./AuxCode/awk/extract_struct_GALAXY_OUTPUT.awk /Users/Bam/ClionProjects/ITH18/h_galaxy_output.i >> ./AuxCode/awk/output/structs.dat