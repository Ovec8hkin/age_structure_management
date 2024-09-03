#' Project Population Forward Multiple Years
#' 
#' A wrapper function around the base `project` function that handles
#' projecting forward multiple years given a removals timeseries and a 
#' recruitment timeseries.
#'
#' @param init_naa numbers-at-age matrix in starting year ([1, nages, nsexes, nregions])
#' @param removals_timeseries vector of removals (either catch of F) of length nyears
#' @param recruitment vector of recruitmene of length nyears
#' @param dem_params list of demographic parameter matrices
#' @param nyears number of projection yeas
#' @param model_option list of additional model options
#'
#' @export project_multi
#'
#' @example
#'
project_multi2 <- function(init_naa, recruitment, dem_params, nyears, model_options){

    model_dimensions <- get_model_dimensions(dem_params$sel)
    nyears <- model_dimensions$nyears
    nages <- model_dimensions$nages
    nsexes <- model_dimensions$nsexes
    nregions <- model_dimensions$nregions
    nfleets <- model_dimensions$nfleets
    # nsurveys <- get_model_dimensions(dem_params$surv_sel)$nsurveys
    nsurveys <- ifelse(model_options$simulate_observations, get_model_dimensions(dem_params$surv_sel)$nfleets, 0)

    caa         = array(NA, dim=c(nyears, nages, nsexes, nregions))
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions))
    naa[1,,,] = init_naa

    recruits    = array(NA, dim=c(nyears+1, 1, 1, nregions))
    recruits[1,,,] = sum(init_naa[1,1,,])

    removals_timeseries <- array(NA, dim=c(nyears, 1, 1, nregions, nfleets))

    naa_proj <- naa[1,,,,drop=FALSE]

    set.seed(model_options$seed)
    for(y in 1:nyears){

        # Subset the demographic parameters list to only the current year
        # and DO NOT drop lost dimensions.
        dp.y <- subset_dem_params(dem_params = dem_params, y, d=1, drop=FALSE)

        # DO HARVEST CONTROL RULE
        # Solve for reference points, F from the HCR,
        # and compute TAC for the next year. Note that
        # selectivity for RP calculations is weighted
        # by terminal year F.

        mp <- model_options$hcr
        ssb <- compute_ssb(naa[y,,,,drop=FALSE], dp.y)
        hcr_parameters <- list(x=ssb[1,1])
        if(is.list(mp$hcr$extra_pars)){
            hcr_parameters <- c(hcr_parameters, mp$hcr$extra_pars)
        }

        hcr_out <- do.call(mp$hcr$func, hcr_parameters)
        removals_timeseries[y,,,,] <- hcr_out

        # mgmt_out <- simulate_TAC(
        #     hcr_F = hcr_out, 
        #     naa = naa_proj, 
        #     recruitment = mean(recruits, na.rm=TRUE), 
        #     joint_sel = subset_matrix(dp.y$sel, 1, d=5, drop=TRUE), 
        #     dem_params = dp.y,
        #     hist_tac = 0,
        #     hcr_options = mp$hcr$extra_options,
        #     options = mp$management
        # )

        
        removals_input <- subset_matrix(removals_timeseries, y, d=1, drop=FALSE)
        # removals_input <- subset_matrix(removals_timeseries, y, d=1, drop=FALSE)
        # region_props <- subset_matrix(model_options$region_apportionment, y, d=1, drop=FALSE)
        # fleet_props <- subset_matrix(model_options$fleet_apportionment, y, d=1, drop=FALSE)

        log_dev <- model_options$recruitment_devs[y+1]
        if(is.function(recruitment)){
            r_y <- do.call(recruitment, c(list(naa=naa[y,,,,drop=FALSE], dem_params=dp.y, logdev=log_dev), model_options$recruitment_pars))
        }else{
            rs <- array(recruitment, dim=c(nyears+1, 1))
            r_y <- subset_matrix(rs, y+1, d=1, drop=FALSE)
        }

        # r_y <- ifelse(!is.null(model_options$recruitment_devs), as.vector(exp(log(r_y)+model_options$recruitment_devs[y+1])), as.vector(r_y))
        r_y <- as.vector(r_y)

        r <- apportion_recruitment_single(
            recruits = as.vector(r_y),
            apportionment = subset_matrix(model_options$recruit_apportionment, y+1, d=1, drop=FALSE),
            nregions = nregions
        )

        rec <- get_annual_recruitment(
            recruitment = r$full_recruitment,
            apportionment = r$rec_props,
            apportion_random = model_options$random_apportion_recruits,
            apportionment_pars = model_options$recruit_apportionment_pars,
            nregions = nregions,
            list(naa=naa[y,,,,drop=FALSE], dem_params=dp.y)
        )
        
        out_vars <- project(
            removals = removals_input,
            dem_params=dp.y,
            prev_naa=naa[y,,,, drop = FALSE],
            recruitment=rec,
            options=model_options
        )

        # update state
        naa[y+1,,,] <- out_vars$naa_tmp
        caa[y,,,] <- out_vars$caa_tmp

        # f[y,,,,] <- out_vars$F_f_tmp
        recruits[y+1,,,] <- rec

        naa_proj <- out_vars$naa_tmp

    }

    return(listN(naa, caa, removals_timeseries))

}
