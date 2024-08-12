rm(list=ls())
library(devtools)
library(tidyverse)
library(doParallel)
library(viridis)
library(ggdist)
library(reshape2)

afscOM_dir <- "~/Desktop/Projects/afscOM"
sablefishmse_dir <- "~/Desktop/Projects/SablefishMSE"

devtools::load_all(afscOM_dir)
devtools::load_all(sablefishmse_dir)

lapply(list.files("R", full.names = TRUE), source)

sable_om <- readRDS(file.path(sablefishmse_dir, "data/sablefish_om_big.RDS"))
joint_selret <- calculate_joint_selret(
    sable_om$dem_params$sel[450,,,,,drop=FALSE],
    sable_om$dem_params$ret[450,,,,,drop=FALSE],
    c(0.80, 0.20)
)
# Define recruitment to occur via historical resampling
assessment <- dget(file.path(sablefishmse_dir, "data/sablefish_assessment_2023.rdat"))
hist_recruits <- assessment$natage.female[,1]*2

nyears <- 500
nages <- 30
nsexes <- 1
nregions <- 1
nfleets <- 1

dimension_names <- list(
    "time" = 1:nyears,
    "age" = 1:nages,
    "sex" = "F",
    "region" = 1,
    "fleet" = 1
)

mortality <- generate_param_matrix(0.15, dimension_names=dimension_names)
maturity <- generate_param_matrix(
    sable_om$dem_params$mat[450,,1,],
    dimension_names=dimension_names,
    by = c("age")
)

weight <- generate_param_matrix(
    sable_om$dem_params$waa[450,,1,],
    dimension_names=dimension_names,
    by = c("age")
)
sex_ratio <- generate_param_matrix(1, dimension_names=dimension_names)
selectivity <- generate_param_matrix(
    joint_selret$sel[,,1,],
    dimension_names=dimension_names,
    by = c("age"),
    include_fleet_dim = TRUE
)
retention <- generate_param_matrix(
    joint_selret$ret[,,1,],
    dimension_names=dimension_names,
    by = c("age"),
    include_fleet_dim = TRUE
)
dmr <- generate_param_matrix(1, dimension_names=dimension_names, include_fleet_dim=TRUE)

dem_params <- list(
    waa=weight,
    mat=maturity,
    mort=mortality,
    sexrat=sex_ratio,
    sel=selectivity,
    ret=retention,
    dmr=dmr,
    surv_sel=retention,
    movement=NA
)

model_dims <- get_model_dimensions(dem_params$sel)

dp_y <- subset_dem_params(dem_params, nyears, d=1, drop=FALSE)

mort <- dp_y$mort
sel <- subset_matrix(dp_y$sel, 1, d=5, drop=TRUE)
ret <- subset_matrix(dp_y$ret, 1, d=5, drop=TRUE)
# F <- 0.107

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

naapr_ref <- get_naapr(dp_y, f=0.1078)
naapr1 <- get_naapr(dp_y, f=0.1312857)
abi2(naapr1, ref=naapr_ref, threshold=0.90)

solve_naa_abi <- function(par, abi, reference, dp){
    f <- par[1]
    naapr <- get_naapr(dp, f=f)
    abi_out <- abi2(naapr, ref=reference, threshold=0.90)
    return(abs(abi-abi_out))
}

nlminb(c(0.1), solve_naa_abi, abi=1.2, reference=naapr_ref, dp=dp_y)

sum(naapr1*25*dp_y$waa*dp_y$mat)
