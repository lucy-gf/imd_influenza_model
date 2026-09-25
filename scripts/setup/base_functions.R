#### BASE FUNCTIONS FOR ANALYSIS ####

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

## GGPLOT THEME ##

theme_lg <- function(gridline_x = FALSE, gridline_y = TRUE) {
  
  gridline <- element_line(
    linewidth = 0.15,
    color = "#999999"
  )
  
  gridline_x <- if (isTRUE(gridline_x)) {
    gridline
  } else {
    element_blank()
  }
  
  gridline_y <- if (isTRUE(gridline_y)) {
    gridline
  } else {
    element_blank()
  }
  
  theme_bw(
  ) +
    theme(
      text = element_text(family = "Helvetica"),
      plot.title = element_text(
        size = 18,
        face = "bold",
        color = "black",
        margin = margin(b = 10, l = 15)
      ),
      plot.subtitle = element_text(
        size = 14,
        color = "#333333",
        margin = margin(b = 10)
      ),
      plot.caption = element_text(
        size = 13,
        color = "#777777",
        margin = margin(t = 15),
        hjust = 0
      ),
      axis.text = element_text(
        size = 12,
        color = "#1F1F1F"
      ),
      axis.title = element_text(
        size = 14,
        color = "#1F1F1F"
      ),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = gridline_x,
      panel.grid.major.y = gridline_y,
      axis.ticks.x = element_line(
        linetype = "solid",
        linewidth = 0.25,
        color = "#333333"
      ),
      axis.ticks.y = element_blank(),
      axis.ticks.length.x = unit(4, units = "pt"),
      strip.background = element_rect(
        fill = "#F0F0F0",
        color = "#777777"
      ),
      strip.text = element_text(
        size = 12,
        # face = "semibold",
        color = "black"
      )
    ) 
  
}

# library(datasauRus)
# 
# datasaurus_dozen %>%
#   filter(dataset %in% c('dino', 'star', 'bullseye', 'slant_down')) %>%
#   ggplot() +
#   geom_point(aes(x, y), size = 2) +
#   scale_y_continuous(limits = c(0, NA), expand = expansion(c(0.00075, 0.025))) +
#   theme_lg() + ggtitle('Datasaurus') + labs(x = '', y = '') +
#   facet_wrap(dataset ~ .)

## FUNCTION TO CONVERT BETWEEN SUBTYPE SYNTAX ##

convert_to_subtype <- function(string_vec){
  
  for(k in 1:length(string_vec)){
    
    text <- string_vec[k]
    
    updated_text <- gsub('flu|pdm09|_', '', text)
    
    string_vec[k] <- toupper(updated_text)
    
  }
  
  string_vec
  
}


## FUNCTION TO ASSIGN AGE GROUPS TO BROAD AGE GROUPS ##

# PURPOSE
# 1. risk: Risk groups (if not using UKHSA data)
# 2. vacc_cov: Vaccination coverage
# 3. vacc_eff: Vaccine efficacy
# 4. susc: Population susceptibility
# 5. reporting rates: Vaccine efficacy

fcn_assign_ages <- function(
    child_value,
    adult_value,
    older_adult_value,
    age_labels_in,
    purpose = NULL
    ){
  
  n_c <- 3
  n_a <- 4
  n_o_a <- length(age_labels_in) - n_c - n_a
  
  vector <- c(rep(child_value, n_c),
              rep(adult_value, n_a),
              rep(older_adult_value, n_o_a))

  return(vector)
  
  }


