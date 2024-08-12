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

dp_y <- subset_dem_params(dem_params, nyears, d=1, drop=FALSE)

mort <- dp_y$mort
sel <- subset_matrix(dp_y$sel, 1, d=5, drop=TRUE)
ret <- subset_matrix(dp_y$ret, 1, d=5, drop=TRUE)
F <- 0.107

naa <- rep(NA, nages)
naa[1] <- 1
zaa <- mort + sel*ret*F
for(a in 2:(nages)){
    naa[a] <- naa[a-1]*exp(-zaa[a-1])
}

sbaa <- array(naa, dim=c(1, nages, 1, 1))*dp_y$waa*dp_y$mat
sbpr <- sum(sbaa)

A_ref <- min(which(cumsum(naa/sum(naa)) > 0.9))

sum(sbaa[(A_ref+1):nages])/sbpr

logistic <- function(k, mid){
    x <- 1:nages
    return(1/(1+exp(-k*(x-mid))))
}

mat_curve <- function(par1, par2){
    x <- 2:(nages+1)
    return(exp(par1+par2*x)/(1+exp(par1+par2*x)))
}

mat_curve(-5.1560, 0.7331)
as.vector(dp_y$mat)

get_maturity_pars <- function(pars, true_mat){
    p1 <- pars[1]
    p2 <- pars[2]
    est_mat <- array(mat_curve(p1, p2), dim=c(1, nages, 1, 1))
    return(sum(abs(est_mat-true_mat)))
}

ests <- nlminb(c(0 ,0), get_maturity_pars, true_mat=dp_y$mat)

ests$par

get_maturity_pars(ests$par, dp_y$mat)

est_mat <- logistic(ests$par[1], ests$par[2])
as.vector(dp_y$mat)

sbaa <- array(naa, dim=c(1, nages, 1, 1))*dp_y$waa*est_mat
sbpr <- sum(sbaa)


spawning_potential <- function(pars, naa, waa, sp){
    k <- pars[1]
    mid <- pars[2]
    est_mat <- array(mat_curve(k, mid), dim=c(1, nages, 1, 1))

    sbaa <- array(naa, dim=c(1, nages, 1, 1))*waa*est_mat
    sbpr <- sum(sbaa)

    A_ref <- min(which(cumsum(naa/sum(naa)) > 0.9))

    old_sp <- sum(sbaa[(A_ref+1):nages])/sbpr
    return(abs(old_sp-sp))
}

ests <- nlminb(c(0, 0), spawning_potential, naa=naa, waa=dp_y$waa, sp=0.35987, upper=c(-4, 10), lower=c(-10, 0.4))
est_mat <- array(mat_curve(ests$par[1], ests$par[2]), dim=c(1, nages, 1, 1))
sbaa <- array(naa, dim=c(1, nages, 1, 1))*dp_y$waa*est_mat
sbpr <- sum(sbaa)
sum(sbaa[(A_ref+1):nages])/sbpr

plot(1:nages, est_mat, type="l")
lines(1:nages, dp_y$mat, col="red")


est_mat2 <- compute_maturity_curve(dp=dp_y, F=0.107, sp=0.90)
plot(1:nages, est_mat2, type="l")
lines(1:nages, dp_y$mat, col="red")

