#### BASE FUNCTIONS FOR ANALYSIS ####

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


