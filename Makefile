
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
SETUP ?= ${CODEDIR}/setup
DATADIR ?= data
INPUTDIR ?= ${DATADIR}/inputs
CMDIR ?= ${DATADIR}/contact_matrix
POPDIR ?= ${DATADIR}/population
OUTDIR ?= output
FIGDIR ?= ${OUTDIR}/figures
DATDIR ?= ${OUTDIR}/data

${OUTDIR} ${DATDIR} ${FIGDIR}:
	mkdir -p $@

RENV = .Rprofile

# build renv/library & other renv infrastructure
${RENV}: install.R 
	 Rscript --vanilla $^

# ages 
ALLAGES ?= 0-4 5-9 10-14 15-19 20-24 25-29 30-34 35-39 40-44 45-49 50-54 55-59 60-64 65-69 70-74 75+
NHSAGES ?= 0-4 5-11 12-17 18-25 26-34 35-49 50-69 70-79 80+

# sensitivity analyses
SENS_ANALYSES ?= base regional

##### INPUTS ###################################################################

${INPUTDIR}/contact_matrix.rds: ${CODEDIR}/load_contact_data.R ${CMDIR}/fitted_matrs_balanced.csv
	$(call R)

${INPUTDIR}/imd_age_pop.rds: ${CODEDIR}/load_pop_data.R ${POPDIR}/imd_2025.xlsx ${POPDIR}/lsoa_to_region.csv
	$(call R)

all_inputs: ${INPUTDIR}/contact_matrix.rds ${INPUTDIR}/imd_age_pop.rds

##### DUMMY DATA ###################################################################

${INPUTDIR}/dummy_flu_data.rds: ${CODEDIR}/dummy_infections.R ${INPUTDIR}/imd_age_pop.rds ${INPUTDIR}/contact_matrix.rds
	$(call R)

all_dummy: ${INPUTDIR}/dummy_flu_data.rds


