
compute_biomass_RP <- function(nyears, dem_params, ...){

    dp <- subset_dem_params(dem_params, 1:nyears, d=1, drop=FALSE)

    naas <- run_simulation(nyears=nyears, dem_params=dp, ...)    

    ssb <- apply(
        apply(naas, 5, \(x) compute_ssb(x[1:nyears,,,,drop=FALSE], dp)),
        1,
        mean
    )
    Bref <- mean(ssb[(nyears*0.8):nyears])
    return(Bref)

}

compute_refnaa <- function(nyears, dem_params, ...){

    dp <- subset_dem_params(dem_params, 1:nyears, d=1, drop=FALSE)

    naas <- run_simulation(nyears=nyears, dem_params=dp, ...)    

    ref_naa <- apply(apply(naas, 5, \(z) z[100,,,]), 1, mean)
    return(ref_naa)

}

beverton_holt <- function(naa, dem_params, h, R0, S0, sigR){
    # set.seed(seed)
    ssb <- compute_ssb(naa, dem_params)[1,1]
    bh <- (4*R0*h*ssb)/((1-h)*R0*(S0/R0) + (5*h - 1)*ssb)
    rec <- rlnorm(1, meanlog = log(bh)-sigR*sigR/2, sdlog=sigR)
    return(rec)
}

compute_maturity_curve <- function(dp, F, sp){

    get_naapr <- function(dp, f){

        mort <- dp$mort
        sel <- subset_matrix(dp$sel, 1, d=5, drop=TRUE)
        ret <- subset_matrix(dp$ret, 1, d=5, drop=TRUE)
        F <- f

        naa <- rep(NA, nages)
        naa[1] <- 1
        zaa <- mort + sel*ret*F
        for(a in 2:(nages)){
            naa[a] <- naa[a-1]*exp(-zaa[a-1])
        }
        return(naa)
    }

    mat_curve <- function(par1, par2){
        x <- 2:(nages+1)
        return(exp(par1+par2*x)/(1+exp(par1+par2*x)))
    }

    spawning_potential <- function(pars, naa, waa, sp){
        p1 <- pars[1]
        p2 <- pars[2]
        est_mat <- array(mat_curve(p1, p2), dim=c(1, nages, 1, 1))

        sbaa <- array(naa, dim=c(1, nages, 1, 1))*waa*est_mat
        sbpr <- sum(sbaa)

        A_ref <- min(which(cumsum(naa/sum(naa)) > 0.9))

        old_sp <- sum(sbaa[(A_ref+1):nages])/sbpr
        return(abs(old_sp-sp))
    }

    dp_y <- subset_dem_params(dp, 1, d=1, drop=FALSE)
    naapr <- get_naapr(dp=dp_y, f=F)

    ests <- nlminb(c(0, 0), spawning_potential, naa=naapr, waa=dp_y$waa, sp=sp, upper=c(-4, 10), lower=c(-10, 0.4))
    est_mat <- array(mat_curve(ests$par[1], ests$par[2]), dim=c(1, nages, 1, 1))
    return(est_mat)

}

