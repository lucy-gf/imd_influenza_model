
###### INTENDED OUTPUTS ########################################################

default: localdef

localdef: all_dummy

###### SUPPORT DEFINITIONS #####################################################

# if need to override directories e.g.
-include local.makefile

# convenience make definitions
R = $(strip Rscript $^ $(1) $@)

# analysis directories + build rules
CODEDIR ?= scripts
SETUPDIR ?= ${CODEDIR}/setup
DUMMYDIR ?= ${CODEDIR}/dummy_data
MCMCDIR ?= ${CODEDIR}/mcmc
DATADIR ?= data
INPUTDIR ?= ${DATADIR}/inputs
DUMMYDAT ?= ${DATADIR}/dummy_data
CMDIR ?= ${DATADIR}/contact_matrix
POPDIR ?= ${DATADIR}/population
UKHSADIR ?= ${DATADIR}/ukhsa
OUTDIR ?= output
FIGDIR ?= ${OUTDIR}/figures
DATDIR ?= ${OUTDIR}/data

${OUTDIR} ${DATDIR} ${FIGDIR}:
	mkdir -p $@

RENV = .Rprofile

# build renv/library & other renv infrastructure
${RENV}: install.R 
	 Rscript --vanilla $^

##### INPUTS ###################################################################

${INPUTDIR}/contact_matrix.rds: ${SETUPDIR}/load_contact_data.R ${CMDIR}/fitted_matrs_balanced.csv
	$(call R)

${INPUTDIR}/imd_age_pop.rds: ${SETUPDIR}/load_pop_data.R ${POPDIR}/imd_2025.xlsx ${POPDIR}/lsoa_to_region.csv
	$(call R)

${INPUTDIR}/subtype_years.rds: ${SETUPDIR}/subtype_setup.R ${UKHSADIR}/annual_influenza_2025_2026.ods
	$(call R)

all_inputs: ${INPUTDIR}/contact_matrix.rds ${INPUTDIR}/imd_age_pop.rds ${INPUTDIR}/subtype_years.rds

##### PARAMETERS ###################################################################

${DUMMYDAT}/known_parameters.rds: ${DUMMYDIR}/produce_known_parameters.R ${INPUTDIR}/imd_age_pop.rds ${INPUTDIR}/subtype_years.rds ${POPDIR}/risk_group_population_data.rds
	$(call R)

${DUMMYDAT}/unknown_parameters.rds: ${DUMMYDIR}/produce_unknown_parameters.R ${INPUTDIR}/imd_age_pop.rds ${INPUTDIR}/subtype_years.rds
	$(call R)

all_pars: ${DUMMYDAT}/known_parameters.rds ${DUMMYDAT}/unknown_parameters.rds

##### DUMMY DATA ###################################################################

${DUMMYDAT}/dummy_infections.rds: ${DUMMYDIR}/dummy_infections.R ${INPUTDIR}/imd_age_pop.rds ${INPUTDIR}/contact_matrix.rds ${DUMMYDAT}/known_parameters.rds ${DUMMYDAT}/unknown_parameters.rds
	$(call R)

${DUMMYDAT}/dummy_surveillance.rds: ${DUMMYDIR}/dummy_surveillance.R ${DUMMYDAT}/dummy_infections.rds ${DUMMYDAT}/known_parameters.rds ${DUMMYDAT}/unknown_parameters.rds
	$(call R)

dummy_inf: ${DUMMYDAT}/dummy_infections.rds
all_dummy: ${DUMMYDAT}/dummy_infections.rds ${DUMMYDAT}/dummy_surveillance.rds

















