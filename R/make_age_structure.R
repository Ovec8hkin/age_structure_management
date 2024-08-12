make_age_structure <- function(abi, reference, threshold){
    nages <- length(reference)
    a_MSY <- min(which(cumsum(reference/sum(reference)) > threshold))+1
    naa <- rep(NA, nages)

    refnaa_prop <- reference/sum(reference)

    old_mass <- sum(refnaa_prop[(a_MSY+1):nages])

    new_mass <- old_mass*abi
    mass_change <- ifelse(new_mass > 0, new_mass - old_mass, old_mass + new_mass)

    naa <- refnaa_prop
    naa[(a_MSY+1):nages] <- naa[(a_MSY+1):nages]+(naa[(a_MSY+1):nages]/sum(naa[(a_MSY+1):nages]))*mass_change
    naa[2:a_MSY] <- naa[2:a_MSY]-(naa[2:a_MSY]/sum(naa[2:a_MSY]))*mass_change
    naa <- naa/sum(naa)
    return(naa)
}
