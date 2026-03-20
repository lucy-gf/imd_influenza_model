
select <- dplyr::select
wday <- lubridate::wday
map <- purrr::map

# Suppress summarise info
options(readr.show_col_types = FALSE)
options(dplyr.summarise.inform = FALSE)

na_to_0 <- function(x){
  x[is.na(x)] <- 0
  x
}

susc_vector <- function(x,
                        rep1 = 3,
                        rep2 = 4,
                        rep3 = 2){
  
  x <- unname(unlist(x))
  
  if(length(x) %notin% 2:3){stop('Length not 2 or 3')}
  if(length(x) == 2){x <- c(x[1], 1, x[2])}

  c(rep(x[1],rep1), rep(x[2],rep2), rep(x[3],rep3))
  
}

imd_spline <- function(imd_pars){
  
  if(length(imd_pars) != 2){stop('Input length not 2')}
  
  # IMD 1 is exp(imd_pars[1]) times the reporting rate of IMD 3,
  # IMD 5 is exp(imd_pars[1]) times the reporting rate of IMD 3
  
  x <- exp(imd_pars)
  
  c(x[1], 1 + (x[1] - 1)/2, 1, 1 + (x[2] - 1)/2, x[2])
  
}

simp_labels <- function(string){
  if(string == 'p_sec_input'){return('nssec')}
  string <- gsub('_group','',string)
  string <- gsub('_grp','',string)
  string <- gsub('nssec','',string)
  string <- gsub('p_sec_input_','nsseccode_',string)
  string <- gsub('_p_sec_input','_nsseccode',string)
  string <- gsub('c_','',string)
  string <- gsub('p_','',string)
  string <- gsub('_cd','',string)
  string <- gsub('_nm','',string)
  string <- gsub('_input','',string)
  string <- gsub('_short','',string)
  string <- gsub('_6','',string)
  string <- gsub('_8','',string)
  return(string)
}

variables_from_name <- function(varname){
  
  if(! varname %in% names(variables_and_names)){stop('varname not found')}
  
  out <- variables_and_names[[varname]]
  
  out
}

name_from_variables <- function(vars){
  
  i <- 0; j <- 0
  
  for(k in 1:length(variables_and_names)){
    if(length(variables_and_names[[k]]) == length(vars)){
      if(sum(variables_and_names[[k]] == vars) == length(vars)){
        i <- k
        j <- j + 1
      }
    }
  }
  
  if(j > 1){stop('vars match more than one list entry')}
  if(j == 0){stop("vars don't match any list entries")}
  
  out <- names(variables_and_names)[i]
  
  out
}

variables_and_names <- list(
  'engreg' = c('eng_reg'),
  'pcd' = c('pcd1'),
  'pcdage' = c('pcd1','age_grp'),
  'ageethn' = c('p_age_group','p_ethnicity'),
  'pcdageethn' = c('pcd1','p_age_group','p_ethnicity'),
  'utlaageethn' = c('utla','p_age_group','p_ethnicity'),
  'pcdagehiqualnssec' = c('pcd1','age_grp_8','p_hiqual','p_sec_input'),
  'pcdhousehold' = c('pcd1','hh_size_nm','hh_tenure_nm'),
  'pcdhhsize' = c('pcd1','hh_size_nm'),
  'pcdhhtenure' = c('pcd1','hh_tenure_nm'),
  'pcdagehiqual' = c('pcd1','age_grp_8','p_hiqual'),
  'pcdagenssec' = c('pcd1','age_grp_8','p_sec_input'),
  'pcdethn' = c('pcd1','p_ethnicity'),
  'pcdethntenure' = c('pcd1','p_ethnicity','p_tenure_short'),
  'pcdethnhiqual' = c('pcd1','p_ethnicity','p_hiqual'),
  'pcdtenurenssec' = c('pcd1','p_tenure_short','p_sec_input'),
  'ethnnssec' = c('p_ethnicity','p_sec_input'),
  'pcdethnnssec' = c('pcd1','p_ethnicity','p_sec_input'),
  'utlaethnnssec' = c('utla','p_ethnicity','p_sec_input'),
  'pcdageethnnssec' = c('pcd1','p_age_group','p_ethnicity','p_sec_input'),
  'utlaageethnnssec' = c('utla','p_age_group','p_ethnicity','p_sec_input'),
  'ageethnnssec' = c('p_age_group','p_ethnicity','p_sec_input')
)

## COLORS FOR IMD ANALYSIS ##

# var_cols <- c('#9B7EDE','#832161','#F7B801','#52050A','#BCD2EE','#F35B04')
# names(var_cols) <- fitted_pars

eng_reg_colors <- c("London"="#31688EFF", "North West" = '#CC4678FF',
                    'Yorkshire and The Humber' = '#65156EFF', "North East" = '#006837', 
                    "West Midlands" = '#F89441FF', 'East Midlands' = '#e78ac3',
                    'South West' = '#253494', 'East of England' = '#bd0026', 'South East' = '#c2e699')

colors_p_engreg <- c('East Midlands' = '#e78ac3', 'East of England' = '#bd0026', 'Greater London' = '#31688EFF',
                     'North East' = '#006837', 'North West' = '#CC4678FF', 'South East' = '#c2e699', 
                     'South West' = '#253494', 'West Midlands' = '#F89441FF', 'Yorkshire and the Humber' = '#65156EFF')

urban_rural_colors <- c('Urban' = '#756bb1', 'Rural' = '#e6550d')

imd_quintile_colors <- c('1' = '#7a0177', '2' = '#c51b8a', '3' = '#f768a1', '4' = '#fa9fb5', '5' = '#fcc5c0')

colors_p_urban_rural <- c('Urban' = 'dodgerblue3', 'Rural' = 'green4')

colors_p_sec_input <- c('1' = '#d9f0a3', '2' = '#addd8e', '3' = '#78c679',
                        '4' = '#41ab5d', '5' = '#238443', '6' = '#006837', '7' = '#004529')

colors_p_hiqual <- c('No qualifications' = '#08589e', '1-4 GCSEs' = '#2b8cbe', '5+ GCSEs' = '#4eb3d3',
                        '2+ A levels' = '#7bccc4', 'Degree' = '#a8ddb5', 'Apprentice/vocational' = '#41ae76')

colors_p_ethnicity <- c('Asian' = '#67a9cf', 'Black' = '#ef8a62', 'Mixed' = '#af8dc3',
                        'Other' = '#2166ac', 'White' = '#b2182b')

colors_hh_size_nm <- c('1 person in household' = '#dadaeb', '2 people in household' = '#bcbddc', 
                       '3 people in household' = '#9e9ac8', '4 people in household' = '#807dba', 
                       '5 people in household' = '#6a51a3', '6 or more people in household' = '#4a1486')

colors_hh_tenure_nm <- c('Owned: Owns outright' = 'darkorange2', 
                         'Owned: Owns with a mortgage or loan or shared ownership' = 'orangered2', 
                         'Private rented: Other private rented or lives rent free' = 'lightskyblue2', 
                         'Private rented: Private landlord or letting agency' = 'steelblue2', 
                         'Social rented: Other social rented' = 'violet', 
                         'Social rented: Rents from council or Local Authority' = 'mediumpurple',
                         'Owned: Owns with a mortgage or\nloan or shared ownership' = 'orangered2', 
                         'Private rented: Other private\nrented or lives rent free' = 'lightskyblue2', 
                         'Private rented: Private landlord\nor letting agency' = 'steelblue2', 
                         'Social rented: Rents from council\nor Local Authority' = 'mediumpurple')

colors_p_age_group <- c('0-4' = '#4d004b', '5-9' = '#7a0177',
                        '10-14' = '#ae017e', '15-19' = '#dd3497',
                        '20-24' = '#f768a1', '25-29' = '#fa9fb5',
                        '30-34' = '#fcc5c0', '35-39' = '#fec44f',
                        '40-44' = '#fe9929', '45-49' = '#ec7014',
                        '50-54' = '#ef6548', '55-59' = '#d7301f',
                        '60-64' = '#990000', '65-69' = '#08589e',
                        '70-74' = '#3690c0', '75+' = '#a6bddb')

colors_p_age_group_nhs <- c('0-4' = '#4d004b',
                           '5-11' = '#ae017e', 
                           '12-17' = '#dd3497',
                           '18-25' = '#f768a1', 
                           '26-34' = '#fcc5c0',
                           '35-49' = '#fec44f',
                           '50-69' = '#fe9929',
                           '70-79' = '#ef6548',
                           '80+' = '#990000')

age_colors <- colors_p_age_group_nhs

colors_age_grp <- c('Aged 4 years and under' = '#4d004b', 'Aged 5 to 9 years' = '#7a0177',
                        'Aged 10 to 14 years' = '#ae017e', 'Aged 15 to 19 years' = '#dd3497',
                        'Aged 20 to 24 years' = '#f768a1', 'Aged 25 to 29 years' = '#fa9fb5',
                        'Aged 30 to 34 years' = '#fcc5c0', 'Aged 35 to 39 years' = '#fec44f',
                        'Aged 40 to 44 years' = '#fe9929', 'Aged 45 to 49 years' = '#ec7014',
                        'Aged 50 to 54 years' = '#ef6548', 'Aged 55 to 59 years' = '#d7301f',
                        'Aged 60 to 64 years' = '#990000', 'Aged 65 to 69 years' = '#08589e',
                        'Aged 70 to 74 years' = '#3690c0', 'Aged 75+' = '#a6bddb')

model_colors <- c('utlaageethnnssec' = '#31a354', 'pcd' = '#7a0177', 
                  'pcdage' = '#fe9929', 'pcdageethn' = '#bdc9e1', 
                  'pcdagehiqualnssec' = '#a50f15', 'pcdhousehold' = '#bae4b3',
                  'pcdagenssec' = '#006d2c', 'pcdagehiqual' = '#2b8cbe', 
                  'pcdethn' = '#de2d26', 'pcdhhsize' = '#f768a1', 'pcdhhtenure' = '#d7b5d8',
                  'pcdethntenure' = 'blue4', 'pcdethnhiqual' = '#fec44f',
                  'pcdethnnssec' = '#ec7014', 'pcdageethnnssec' = '#8c96c6',
                  'ageethnnssec' = '#fb9a99', 'engreg' = '#72874EFF') 

model_names <- unlist(unname(lapply(sapply(lapply(names(model_colors), variables_from_name),
                             paste, collapse = '_'),
                             simp_labels)))
for(k in 1:length(model_names)){model_names[k] <- gsub('nsseccode','nssec',model_names[k])}
for(k in 1:length(model_names)){model_names[k] <- gsub('pcd1','pcd',model_names[k])}
names(model_names) <- names(model_colors)

method_shapes <- c('det' = 1, 'prob' = 19)  
method_names <- c('Deterministic','Probabilistic')
names(method_names) <- c('det','prob')

variable_colors <- c('age_grp' = '#31a354', 'hh_size_nm' = '#7a0177', 
                  'hh_tenure_nm' = '#fd8d3c', 'p_ethnicity' = '#bdc9e1', 
                  'p_hiqual' = '#f768a1', 'p_sec_input' = '#2b8cbe', 
                  'p_urban_rural' = '#de2d26')

variable_colors <- c('nssec' = '#66c2a4', 'age' = '#7a0177', 
                     'hh_size' = '#fd8d3c', 'hh_tenure' = '#bdc9e1', 
                     'region' = '#f768a1', 'hiqual' = '#2b8cbe', 
                     'ethnicity' = '#de2d26', 'urban_rural' = '#006d2c')

gender_colors <- c('Female' = '#31a354', 'Male' = '#7a0177')


imd_model_colors <- c(imd_quintile_colors,rev(c('#c6dbef','#6baed6','#4292c6','#2171b5','#08519c')))
names(imd_model_colors) <- paste0(rep(1:5,2), '.', c(rep('National-level',5),rep('Regional-level',5)))

imd_model_labels <- paste0(c(rep('National-level',5),rep('Regional-level',5)),
                           ', IMD ', rep(1:5,2))





