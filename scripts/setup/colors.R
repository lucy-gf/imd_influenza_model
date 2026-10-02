## FUNCTIONS, COLORS, VARIABLE NAMES, ETC. ##

select <- dplyr::select
wday <- lubridate::wday
map <- purrr::map

# Suppress summarise info
options(readr.show_col_types = FALSE)
options(dplyr.summarise.inform = FALSE)

## COLORS FOR IMD ANALYSIS ##

# var_cols <- c('#9B7EDE','#832161','#F7B801','#52050A','#BCD2EE','#F35B04')
# names(var_cols) <- fitted_pars

flu_strain_colors <- c('flu_a' = '#81559B', 'flu_b' = '#6699CC')#'#B2EF9B')
flu_strain_names <- c('flu_a' = 'Flu A', 'flu_b' = 'Flu B')
subtype_colors <- c('A' = '#81559B', 'AH1N1' = '#FF3C38', 'AH3N2' = '#FF8C42',
                   'AUNSUB' = '#A23E48', 'B' = '#6699CC')
subtype_names <- c('A' = 'Flu A', 'AH1N1' = 'Flu A (H1N1)', 'AH3N2' = 'Flu A (H3N2)', 'AUNSUB' = 'Flu A (unsubtyped)', 'B' = 'Flu B')

flu_subtype_colors <- c('flu_a' = '#81559B', 'flu_a_h1n1pdm09' = '#FF3C38', 'flu_a_h3n2' = '#FF8C42',
                        'flu_a_unsubtyped' = '#A23E48', 'flu_b' = '#6699CC')
flu_subtype_names <- c('flu_a' = 'Flu A', 'flu_a_h1n1pdm09' = 'Flu A H1N1pmd09', 'flu_a_h3n2' = 'Flu A H3N2',
                       'flu_a_unsubtyped' = 'Flu A Unsubtyped', 'flu_b' = 'Flu B')

risk_colors <- c('low_risk' = '#D90368', 'high_risk' = '#F75C03')

season_colors <- c('2023' = '#5D576B', '2024' = '#99E1D9', '2025' = '#F7567C')

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





