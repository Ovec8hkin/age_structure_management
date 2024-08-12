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

reference_points <- calculate_ref_points(
    nages=nages,
    mort = dp_y$mort,
    mat = dp_y$mat,
    waa = dp_y$waa,
    sel = dp_y$sel[,,,,1],
    ret = dp_y$ret[,,,,1],
    avg_rec = mean(hist_recruits),
    spr_target = 0.40
)


Fref <- reference_points$Fref
Bref <- reference_points$Bref


model_options <- list()
model_options$removals_input = "F"
model_options$fleet_apportionment = NULL
model_options$recruit_apportionment = NULL
model_options$random_apportion_recruits = FALSE
model_options$do_recruits_move = FALSE
model_options$simulate_observations = FALSE
model_options$seed = 1120

init_naa <- sable_om$init_naa[,,1,]


refnaa <- compute_refnaa(dims=model_dims, nyears=100, nsims=100, dem_params=dem_params, init_naa=init_naa, model_options=model_options, f=Fref, steepness=1, sigR=0, seed=1120)

# refnaa <- apply(model_grid, 1, function(x){ print(x); compute_refnaa(dims=model_dims, nyears=100, nsims=100, dem_params=dem_params, init_naa=init_naa, model_options=model_options, f=x[1], steepness=x[2], sigR=x[3], seed=1120)})
nsims <- 20
# recruitment strengths: h=0.5, 0.6, 0.7, 0.8, 0.9, 1.0
#steepness <- c(0.50, 0.60, 0.70, 0.80, 0.90, 1.0)
steepness <- c(0.75, 0.90, 1.0)
# recruitment variability: sigR=0.4, sigR=0.8
sigma_r <- c(0.10, 0.40, 0.80)
# fishing mortality F = 0, 0.25F_MSY, 0.5FMSY, 0.75F_MSY, 1.0F_MSY, 1.5F_MSY
#fishing_mort <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5)
fishing_mort <- c(0.5, 1.0, 1.5)
# population scale: B/B_MSY = 0.25, 0.5, 0.75, 1.0, 1.5 
#depletion = c(0.25, 0.50, 0.75, 1, 1.25, 1.50)
depletion = c(0.20, 0.5, 0.75, 1.0, 1.25, 1.50)
# starting population structure: abi=0.1, 0.2, ..., 3.0
abis <- seq(0.25, 2, 0.25)#seq(0.25, 3, 0.25)
# old spawning potential
sp <- seq(0.3, 0.9, 0.2) 

# sp_refs <- apply(as.matrix(sp), 1, function(x){
#     mat <- as.vector(compute_maturity_curve(dp_y, F=Fref, sp=x))
#     d <- dem_params
#     d$mat <- generate_param_matrix(mat, dimension_names=dimension_names, by="age")
#     ref_naa <- compute_refnaa(dims=model_dims, nyears=100, nsims=50, dem_params=d, init_naa=init_naa, model_options=model_options, f=Fref, steepness=1, sigR=0, seed=1120)
#     return(ref_naa)
# })

# simulations: nsims=20
seeds <- sample(1:1e6, nsims)

model_grid <- expand.grid(steepness, fishing_mort, depletion, abis, sp, seeds)
nyears <- 100
sim_objects = list()
for(i in 1:nrow(model_grid)){
    print(paste0(i, "/", nrow(model_grid)))
    h <- model_grid[i, 1]
    model_options$recruitment_pars <- list(
        h = h,
        R0 = 25,
        S0 = 300,
        sigR = 0.1
    )

    # Set maturity curve
    mat <- as.vector(compute_maturity_curve(dp_y, F=Fref, sp=model_grid[i, 5]))
    dem_params$mat <- generate_param_matrix(mat, dimension_names=dimension_names, by="age")

    dp_y <- subset_dem_params(dem_params, nyears, d=1, drop=FALSE)
    reference_points <- calculate_ref_points(
        nages=nages,
        mort = dp_y$mort,
        mat = dp_y$mat,
        waa = dp_y$waa,
        sel = dp_y$sel[,,,,1],
        ret = dp_y$ret[,,,,1],
        avg_rec = mean(hist_recruits),
        spr_target = 0.40
    )


    Fref <- reference_points$Fref
    Bref <- reference_points$Bref

    f <- Fref*model_grid[i,2]
    # ref_naa <- refnaa[,(grid$F == f & grid$h == h & grid$sigR == 0.1)]
    ref_naa <- refnaa

    abi <- model_grid[i,4]
    init_naa <- make_age_structure(abi=abi, reference = ref_naa, threshold=0.90)

    #Bref <- grid$B[grid$F == f & grid$h == h & grid$sigR == 0.1]

    desired_ssb <- Bref*model_grid[i,3]
    multiplier <- desired_ssb/sum(init_naa*dem_params$waa[1,,,]*dem_params$mat[1,,,])

    set.seed(model_grid[i, 6])
    model_options$recruitment_devs <- rnorm(nyears+1, mean=0, sd=0.1)

    model_options$hcr <- setup_mp_options()
    model_options$hcr$hcr <- list(
        func = threshold_f,
        extra_pars = list(
            f_min = 0,
            f_max = f,
            lrp = 0.05*Bref,
            urp = Bref
        ),
        extra_options = list(
            max_stability = NA,
            harvest_cap = NA
        ),
        units = "F"
    )
    
    dp <- subset_dem_params(dem_params, 1:nyears, d=1, drop=FALSE)
    sim <- project_multi2(
        init_naa = array(multiplier*init_naa, dim=c(1, nages, nsexes, nregions)),
        recruitment = beverton_holt, 
        dem_params = dp, 
        nyears = nyears,
        model_options = model_options
    )
    sim_objects[[i]] <- list(naa=sim$naa)

}

extra_columns <- expand.grid(
    steepness=steepness, 
    fishing=fishing_mort, 
    scale=depletion, 
    structure=abis,
    sp=sp,
    sim=seeds
) %>% as.data.frame

sim_naas <- bind_mse_outputs(sim_objects, c("naa"), extra_columns) 

# write_csv(sim_naas, file=file.path("data/sim_naas1.csv"))

mats <- apply(as.matrix(sp), 1, function(x){
    mat <- as.vector(compute_maturity_curve(dp_y, F=Fref, sp=x))
    return(mat)
})

d <- dem_params
d$mat <- generate_param_matrix(mats, dimension_names=c(dimension_names[c("time", "age", "sex", "region")], list("fleet"=sp)), by=c("age", "fleet"))

dp <- subset_dem_params(d, 1:nyears, d=1, drop=FALSE)

dimnames(dp$waa) <- list("time"=1:nyears, "age"=1:nages, "sex"="fem", "region"=1)
dimnames(dp$mat) <- list("time"=1:nyears, "age"=1:nages, "sex"="fem", "region"=1, "sp"=sp)
dimnames(dp$mort) <- list("time"=1:nyears, "age"=1:nages, "sex"="fem", "region"=1)
dimnames(dp$sel) <- list("time"=1:nyears, "age"=1:nages, "sex"="fem", "region"=1, "fleet"=1)
dimnames(dp$ret) <- list("time"=1:nyears, "age"=1:nages, "sex"="fem", "region"=1, "fleet"=1)




ssb_data <- sim_naas %>% as_tibble() %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    # select(time, age, sex, region, sim, steepness, fishing, scale, structure, L1, value) %>%
    mutate(sex = case_when(sex == 1 ~ "fem", TRUE ~ "mal")) %>%
    filter(time <= nyears) %>%
    left_join(reshape2::melt(dp$waa, value.name="waa"), by=c("time", "age", "sex", "region")) %>%
    left_join(reshape2::melt(dp$mat, value.name="mat"), by=c("time", "age", "sex", "region", "sp")) %>%
    select(time, age, sex, region, sim, steepness, fishing, scale, structure, sp, L1, value, waa, mat) %>%
    mutate(
        ssbaa = value*waa*mat
    ) %>%
    select(-c(L1, value, waa, mat)) %>%
    group_by(time, sim, steepness, fishing, scale, structure, sp) %>%
    summarise(ssb=sum(ssbaa))

abi_data <- sim_naas %>% as_tibble() %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    # select(time, age, sex, region, sim, steepness, fishing, scale, structure, L1, value) %>%
    mutate(sex = case_when(sex == 1 ~ "fem", TRUE ~ "mal")) %>%
    filter(time <= nyears) %>%
    group_by(time, sex, region, L1, steepness, fishing, scale, structure, sp, sim) %>%
    summarise(
        abi=abi2(value, ref=refnaa, threshold=0.9)
    ) %>%
    select(time, steepness, fishing, scale, structure, sim, sp, abi)


# Recovery Time: number of years required for the population to
# recover to a reference spawning biomasslevel
recovery_time <- ssb_data %>%
    filter(time > 5) %>%
    ungroup() %>%
    group_by(steepness, sp) %>%
    mutate(Bref = median(ssb[which(time == max(time))])) %>%
    group_by(sim, steepness, fishing, scale, structure, sp) %>%
    mutate(
        # recovered = case_when(
        #     scale <= 1 ~ ((ssb >= Bref) & (ssb[which(time == max(time))] > Bref | (abs(ssb[which(time == max(time))]-Bref)/Bref) < 0.12)),
        #     scale > 1 ~ ((ssb <= Bref) & (ssb[which(time == max(time))] < Bref | (abs(ssb[which(time == max(time))]-Bref)/Bref) < 0.12)),
        # ),
        recovered = abs(ssb-Bref)/Bref < 0.05,
        recovery_time = min(min(time[which(recovered == TRUE)]), max(time)),
        std_recovery_time = recovery_time/6
    ) %>%
    group_by(steepness, fishing, scale, structure, sp) %>%
    median_qi(recovery_time, std_recovery_time, .width=c(0.95)) %>%
    mutate(
        steepness = factor(steepness),
        fishing = factor(fishing),
        scale = factor(scale),
        structure = factor(structure),
        sp = factor(sp)
    )

ggplot(recovery_time)+
    geom_point(aes(x=fishing, y=recovery_time,color=structure, group=structure))+
    geom_line(aes(x=fishing, y=recovery_time, color=structure, group=structure))+
    facet_grid(rows=vars(steepness), cols=vars(scale), scales="free_y")+
    scale_y_continuous(limits=c(0, 24), expand=c(0, 0))+
    theme_bw()

ggplot(recovery_time)+
    geom_point(aes(x=fishing, y=std_recovery_time,color=structure, group=structure))+
    geom_line(aes(x=fishing, y=std_recovery_time, color=structure, group=structure))+
    facet_grid(rows=vars(steepness), cols=vars(scale), scales="free_y")+
    scale_y_continuous(limits=c(0, 4), expand=c(0, 0))+
    theme_bw()


ggplot(recovery_time %>% filter(fishing == 1))+
    geom_point(aes(x=structure, y=recovery_time, color=scale, group=scale))+
    geom_line(aes(x=structure, y=recovery_time, color=scale, group=scale))+
    facet_grid(cols=vars(steepness), rows=vars(sp), scales="free_y")+
    scale_y_continuous(limits=c(0, 70), expand=c(0, 0))+
    theme_bw()


ggsave("~/Desktop/sablefish_age_structure_analyses.png", width=11*1.5, height=9*1.5, units=c("in"))

# SSB Timeseries: timeseries plots of spawning biomass
ssb_tseries <- ssb_data %>%
    mutate(
        steepness = factor(steepness),
        fishing = factor(fishing),
        scale = factor(scale),
        structure = factor(structure),
        sp = factor(sp)
    ) %>%
    group_by(steepness, sp) %>%
    mutate(Bref = median(ssb[which(time == max(time))])) %>%
    ungroup() %>%
    group_by(time, steepness, fishing, scale, structure, sp) %>%
    median_qi(ssb, Bref, .width=c(0.50))

ggplot(ssb_tseries %>% filter(steepness == 0.75), aes(x=time, y=ssb))+
    geom_line(aes(color=structure, group=interaction(fishing, structure)))+
    geom_hline(data=ssb_tseries %>% filter(steepness == 0.75), aes(yintercept=Bref))+
    # geom_hline(yintercept=Bref*0.935, linetype="dashed")+
    # geom_hline(yintercept=Bref*1.065, linetype="dashed")+
    facet_grid(row=vars(sp), cols=vars(scale), scales="free_y")+
    coord_cartesian(xlim=c(0, 70))


# Spawning Biomass Stability: average annual variation in spawning
# biomass across productivity, maximum fishing pressure, and stock
# status relative to B_MSY
ssb_stability <- ssb_data %>%
    group_by(sim, steepness, fishing, structure, scale, sp) %>%
    summarise(ssb_aav = aav(ssb)) %>%
    group_by(steepness, fishing, structure, scale, sp) %>%
    median_qi(ssb_aav, .width=c(0.95)) %>%
    mutate(
        steepness = factor(steepness),
        fishing = factor(fishing),
        scale = factor(scale),
        structure = factor(structure),
        sp = factor(sp)
    )

ggplot(ssb_stability)+
    geom_point(aes(x=fishing, y=ssb_aav, color=structure, group=structure))+
    geom_line(aes(x=fishing, y=ssb_aav, color=structure, group=structure))+
    facet_grid(rows=vars(steepness), cols=vars(scale), scales="free_y")+
    scale_y_continuous(limits=c(0, 0.05), expand=c(0, 0))+
    theme_bw()

ggplot(ssb_stability %>% filter(fishing == 1))+
    geom_point(aes(x=structure, y=ssb_aav, color=scale, group=scale))+
    geom_line(aes(x=structure, y=ssb_aav, color=scale, group=scale))+
    facet_grid(cols=vars(steepness), rows=vars(sp), scales="free_y")+
    scale_y_continuous(limits=c(0, 0.05), expand=c(0, 0))+
    theme_bw()

# ABI Stability: average annual variation in ABI
# across productivity, maximum fishing pressure, 
# and stock status relative to B_MSY
abi_stability <- abi_data %>%
    group_by(sim, steepness, fishing, structure, scale, sp) %>%
    summarise(abi_aav = aav(abi)) %>%
    group_by(steepness, fishing, structure, scale, sp) %>%
    median_qi(abi_aav, .width=c(0.95)) %>%
    mutate(
        steepness = factor(steepness),
        fishing = factor(fishing),
        scale = factor(scale),
        structure = factor(structure),
        sp = factor(sp)
    )

ggplot(abi_stability)+
    geom_point(aes(x=fishing, y=abi_aav, color=structure, group=structure))+
    geom_line(aes(x=fishing, y=abi_aav, color=structure, group=structure))+
    facet_grid(rows=vars(steepness), cols=vars(scale), scales="free_y")+
    scale_y_continuous(limits=c(0, 0.1), expand=c(0, 0))+
    theme_bw()

ggplot(abi_stability %>% filter(fishing == 1))+
    geom_point(aes(x=structure, y=abi_aav, color=scale, group=scale))+
    geom_line(aes(x=structure, y=abi_aav, color=scale, group=scale))+
    facet_grid(cols=vars(steepness), rows=vars(sp), scales="free_y")+
    scale_y_continuous(limits=c(0, 0.05), expand=c(0, 0))+
    theme_bw()


##########

recs <- sim_naas %>% as_tibble() %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    # select(time, age, sex, region, sim, steepness, fishing, scale, structure, L1, value) %>%
    mutate(sex = case_when(sex == 1 ~ "fem", TRUE ~ "mal")) %>%
    filter(age == 1) %>% 
    group_by(time, steepness, scale, structure, sp) %>%
    median_qi(value) 

ggplot(recs %>% filter(steepness == 0.75), aes(x=time, y=value))+
    geom_line(aes(color=structure))+
    facet_grid(cols=vars(scale), rows=vars(sp))+
    coord_cartesian(ylim=c(0, 100))
