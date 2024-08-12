run_simulation <- function(dims, nyears, nsims, dem_params, init_naa, model_options, f, steepness, sigR, seed=1120){
    
    model_options$recruitment_pars <- list(
        h = steepness,
        R0 = 25,
        S0 = 300,
        sigR = sigR
    )

    f_timeseries <- array(f, dim=c(nyears, 1, 1, 1, 1))

    set.seed(seed)
    seeds <- sample(1:1e6, nsims)

    naas <- array(NA, dim=c(nyears+1, dims$nages, dims$nsexes, dims$nregions, nsims))
    for(i in 1:nsims){
        set.seed(seeds[i])
        # recruitment_timeseries <- rlnorm(nyears+1, meanlog = log(mean(hist_recruits/2)), sdlog=0.1)
        # model_options$recruitment_devs <- rnorm(nyears+1, mean=0, sd=sigR)
        sim <- project_multi(
            init_naa = init_naa, 
            removals_timeseries = f_timeseries, 
            recruitment = beverton_holt, 
            dem_params = dem_params, 
            nyears = nyears, 
            model_options = model_options
        )
        naas[,,,,i] <- sim$naa
    }
    return(naas)
}

