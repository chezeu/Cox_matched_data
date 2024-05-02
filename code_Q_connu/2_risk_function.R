

###1. Finding the risk set R(t) given some time_of_interest

GetRiskSet <- function(time_of_interest, time_vector) {
  
  return( which(((time_vector >= time_of_interest) )))
  
}

#####
#position of times before time tilde_d

observe <- function(time_of_interest, time_vector, event_vector) {
  return(which( (time_vector <= time_of_interest) & (event_vector == 1) ))
}
####

#number of risk at time tilde_d

risk_nber <- function(time_of_interest, time_vector) {
  return( sum( as.numeric( (time_vector >= time_of_interest)  )))
}
####
#####
#position of times before time time_of_interest

observe_Nd <- function(time_of_interest, time_vector) {
  return(which( (time_vector <= time_of_interest )  ))
}
####
