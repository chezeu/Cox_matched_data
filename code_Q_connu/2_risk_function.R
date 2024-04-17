

###1. Finding the risk set R(t) given some time t

GetRiskSet <- function(time_of_interest, time_vector, event_vector) {
  
  return(which(((time_vector == time_of_interest & event_vector == 1) | (time_vector > time_of_interest))))
  
}


#####
#number of observations before time t

observe <- function(time_of_interest, time_vector, event_vector) {
  return(which( (time_vector <= time_of_interest) ))
}
####
