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

mortality <- generate_param_matrix(0.11, dimension_names=dimension_names)
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
    avg_rec = mean(hist_recruits)/2,
    spr_target = 0.40
)
F40 <- reference_points$Fref
B40 <- reference_points$Bref

sp <- c(0.3, 0.5, 0.7)

B0s <- sapply(sp, function(x){
    
    d <- dem_params
    mat <- as.vector(compute_maturity_curve(dp_y, F=F40, sp=x))
    d$mat <- generate_param_matrix(mat, dimension_names=dimension_names, by="age")
    dp_y <- subset_dem_params(d, 1, d=1, drop=FALSE)

    reference_points <- calculate_ref_points(
        nages=nages,
        mort = dp_y$mort,
        mat = dp_y$mat,
        waa = dp_y$waa,
        sel = dp_y$sel[,,,,1],
        ret = dp_y$ret[,,,,1],
        avg_rec = mean(hist_recruits)/2,
        spr_target = 1
    )
    return(reference_points$Bref)
})

B40s <- sapply(sp, function(x){
    
    d <- dem_params
    mat <- as.vector(compute_maturity_curve(dp_y, F=F40, sp=x))
    d$mat <- generate_param_matrix(mat, dimension_names=dimension_names, by="age")
    dp_y <- subset_dem_params(d, 1, d=1, drop=FALSE)

    reference_points <- calculate_ref_points(
        nages=nages,
        mort = dp_y$mort,
        mat = dp_y$mat,
        waa = dp_y$waa,
        sel = dp_y$sel[,,,,1],
        ret = dp_y$ret[,,,,1],
        avg_rec = mean(hist_recruits)/2,
        spr_target = 0.4
    )
    return(reference_points$Bref)
})

F40s <- sapply(sp, function(x){
    
    d <- dem_params
    mat <- as.vector(compute_maturity_curve(dp_y, F=F40, sp=x))
    d$mat <- generate_param_matrix(mat, dimension_names=dimension_names, by="age")
    dp_y <- subset_dem_params(d, 1, d=1, drop=FALSE)

    reference_points <- calculate_ref_points(
        nages=nages,
        mort = dp_y$mort,
        mat = dp_y$mat,
        waa = dp_y$waa,
        sel = dp_y$sel[,,,,1],
        ret = dp_y$ret[,,,,1],
        avg_rec = mean(hist_recruits)/2,
        spr_target = 0.4
    )
    return(reference_points$Fref)
})

model_options <- list()
model_options$removals_input = "F"
model_options$fleet_apportionment = NULL
model_options$recruit_apportionment = NULL
model_options$random_apportion_recruits = FALSE
model_options$do_recruits_move = FALSE
model_options$simulate_observations = FALSE
model_options$seed = 1120

init_naa <- sable_om$init_naa[,,1,]


naa40 <- compute_refnaa(dims=model_dims, nyears=60, nsims=30, dem_params=dem_params, init_naa=init_naa, model_options=model_options, f=F40, steepness=0.75, sigR=0, seed=1120)
naa00 <- compute_refnaa(dims=model_dims, nyears=100, nsims=100, dem_params=dem_params, init_naa=init_naa, model_options=model_options, f=0, steepness=0.75, sigR=0, seed=1120)


# refnaa <- apply(model_grid, 1, function(x){ print(x); compute_refnaa(dims=model_dims, nyears=100, nsims=100, dem_params=dem_params, init_naa=init_naa, model_options=model_options, f=x[1], steepness=x[2], sigR=x[3], seed=1120)})
nsims <- 20
# recruitment strengths: h=0.5, 0.6, 0.7, 0.8, 0.9, 1.0
steepness <- c(0.75)#, 0.90, 1.0)
# steepness <- c(0.75, 0.90, 1.0)
# recruitment variability: sigR=0.4, sigR=0.8
# sigma_r <- c(0.10, 0.40, 0.80)
sigma_r <- c(0.4, 0.8, 1.2)
# fishing mortality F = 0, 0.25F_MSY, 0.5FMSY, 0.75F_MSY, 1.0F_MSY, 1.5F_MSY
#fishing_mort <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5)
fishing_mort <- c(1.0)
# population scale: B/B_MSY = 0.25, 0.5, 0.75, 1.0, 1.5 
depletion = c(1)
# depletion = c(0.20, 0.5, 0.75, 1.0, 1.25, 1.50)
# starting population structure: abi=0.1, 0.2, ..., 3.0
abis <- seq(0.2, 2, 0.1)
# old spawning potential
sp <- seq(0.3, 0.7, 0.2)
# operating models
oms <- c(0) 


get_reference_naa <- function(s, h, sigma_r){
    F40 <- F40s[match(s, sp)]
    ref_naa <- compute_refnaa(
        dims=model_dims, 
        nyears=60, 
        nsims=30, 
        dem_params=dem_params, 
        init_naa=init_naa, 
        model_options=model_options, 
        f=F40, 
        steepness=h, 
        sigR=sigma_r, 
        seed=1120
    )
    return(list(ref_naa))
}

ref_naas <- expand.grid(sp, steepness, sigma_r) %>% as_tibble() %>%
    rowwise() %>%
    mutate(ref_naa = get_reference_naa(Var1, Var2, Var3))

# simulations: nsims=20
seeds <- sample(1:1e6, nsims)


model_grid <- expand.grid(steepness, fishing_mort, depletion, abis, sp, seeds, sigma_r, oms)
nyears <- 350
sim_objects = list()

for(i in 1:nrow(model_grid)){
    print(paste0(i, "/", nrow(model_grid)))

    d <- dem_params
    dp_y <- subset_dem_params(d, 1, d=1, drop=FALSE)

    # Get SPR40 reference points, since they vary with the maturity curve
    # reference_points <- calculate_ref_points(
    #     nages=nages,
    #     mort = dp_y$mort,
    #     mat = dp_y$mat,
    #     waa = dp_y$waa,
    #     sel = dp_y$sel[,,,,1],
    #     ret = dp_y$ret[,,,,1],
    #     avg_rec = mean(hist_recruits),
    #     spr_target = 0.40
    # )

    # Fref <- reference_points$Fref
    # Bref <- reference_points$Bref

    h <- model_grid[i, 1]
    model_options$recruitment_pars <- list(
        h = h,
        R0 = 25/2,
        S0 = B0s[match(model_grid[i, 5], sp)],
        sigR = model_grid[i, 7]
    )

    F40 <- F40s[match(model_grid[i, 5], sp)]
    B40 <- B40s[match(model_grid[i, 5], sp)]

    Fref <- compute_F_for_abi(abi=model_grid[i,4], F40, dp=dp_y)
    f <- Fref

    # Set maturity curve
    mat <- as.vector(compute_maturity_curve(dp_y, F=f, sp=model_grid[i, 5]))
    d$mat <- generate_param_matrix(mat, dimension_names=dimension_names, by="age")

    # ref_naa <- compute_refnaa(dims=model_dims, nyears=60, nsims=30, dem_params=dem_params, init_naa=init_naa, model_options=model_options, f=F40, steepness=h, sigR=model_grid[i, 7], seed=1120)
    # ref_naa <- ref_naas %>% filter(Var1 == model_grid[i, 5], Var2 == h, Var3 == model_grid[i, 7])  %>% pull(ref_naa) %>% unlist

    abi <- model_grid[i,4]

    set.seed(model_grid[i, 6])
    model_options$recruitment_devs <- rep(NA, nyears)
    # if(model_grid[i, 8] == 1){
    #     model_options$recruitment_devs[75:90] <- rnorm(15, mean=-2.3, sd=0.1)
    # }else if(model_grid[i, 8] == 2){
    #     d$mort[100:110,,,] <- d$mort[100:110,,,]*2
    # }

    model_options$recruitment_devs[175:250] <- rnorm(75, mean=-3, sd=sigma_r) 
    model_options$recruitment_devs[c(200, 225)] <- rnorm(2, mean = 1, sd=sigma_r) 
    # d$mort[100:105,,,] <- d$mort[100:105,,,]*3

    desired_ssb <- B40
    multiplier <- desired_ssb/sum(init_naa*d$waa[1,,,]*d$mat[1,,,])

    model_options$seed <- model_grid[i, 6]
    model_options$hcr <- setup_mp_options()
    model_options$hcr$hcr <- list(
        func = threshold_f,
        extra_pars = list(
            f_min = 0,
            f_max = f,
            lrp = 0.05*B40,
            urp = B40
        ),
        extra_options = list(
            max_stability = NA,
            harvest_cap = NA
        ),
        units = "F"
    )
    
    dp <- subset_dem_params(d, 1:nyears, d=1, drop=FALSE)
    # tic()
    sim <- project_multi2(
        init_naa = array(multiplier*init_naa, dim=c(1, nages, 1, 1)),
        recruitment = beverton_holt, 
        dem_params = dp, 
        nyears = nyears,
        model_options = model_options
    )
    # toc()
    sim_objects[[i]] <- list(naa=sim$naa, caa=sim$caa, f=sim$removals_timeseries)

}

expand.grid(a=abis, s=sp) %>% as_tibble() %>%
    rowwise() %>%
    mutate(
        F = compute_F_for_abi(abi=a, F40s[match(s, sp)], dp=dp_y),
        spr = compute_spr(
            nages=nages,
            mort = dp_y$mort,
            mat = compute_maturity_curve(dp_y, F=F, sp=s),
            waa = dp_y$waa,
            sel = dp_y$sel[,,,,1],
            ret = dp_y$ret[,,,,1],
            F = F
        )
    )

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

extra_columns <- expand.grid(
    steepness=steepness, 
    structure=abis,
    sp=sp,
    sim=seeds,
    sigma_r=sigma_r,
    om = oms
) %>% as.data.frame

sim_naas <- bind_mse_outputs(sim_objects, c("naa"), extra_columns) %>% as_tibble()
sim_fs <- bind_mse_outputs(sim_objects, c("f"), extra_columns) %>% as_tibble() %>%
    group_by(structure, sp, steepness, sigma_r) %>%
    summarise(f=mean(value))


sim_fs %>%
    mutate(
        spr = compute_spr(
            nages=nages,
            mort = dp_y$mort,
            mat = dp_y$mat,
            waa = dp_y$waa,
            sel = dp_y$sel[,,,,1],
            ret = dp_y$ret[,,,,1],
            F = f
        ),
        frel = f/F40
    ) %>%
    group_by(sigma_r, sp, steepness, structure) %>%
    median_qi(spr, f, frel)


sim_fs %>% group_by(Var1, structure, sp, steepness, sigma_r) %>%
    summarise(f = mean(f))

# sim_naas <- sim_naas %>% filter(Var1 < 150, om != 2)
sim_caas <- bind_mse_outputs(sim_objects, c("caa"), extra_columns) %>% as_tibble()

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


fdat <- sim_fs %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    filter(time == 100, sigma_r == 0.4) %>%
    group_by(steepness, structure, sp, sigma_r) %>%
    mutate(
        F = compute_F_for_abi(abi=structure, F40s[match(sp, c(0.3, 0.5, 0.7))], dp=dp_y),
        spr = compute_spr(
            nages=nages,
            mort = dp_y$mort,
            mat = compute_maturity_curve(dp_y, F=F, sp=sp),
            waa = dp_y$waa,
            sel = dp_y$sel[,,,,1],
            ret = dp_y$ret[,,,,1],
            F = F
        )
    ) %>%
    group_by(sigma_r, sp, steepness, structure) %>%
    median_qi(spr, F) %>%
    select(sigma_r, sp, structure, spr, F)


lm1.out <- lm(F ~ structure*sp, data=fdat)
lm1.pred <- predict(lm1.out, structure=abi)
lm2.out <- lm(log(F) ~ structure*sp, data=fdat, offset=sp)
lm2.pred <- exp(predict(lm2.out, structure=abi))

fdat %>% mutate(pred1=lm1.pred, pred2=lm2.pred, pred3=exp(-(1.78304)*structure) + (-0.03238854)*structure) %>%
    ggplot(aes(group=sp))+
        geom_point(aes(x=structure, y=F, color=factor(sp)))+
        geom_line(aes(x=structure, y=F, color=factor(sp)))+
        geom_line(aes(x=structure, y=pred1), color="red")+
        geom_line(aes(x=structure, y=pred2), color="blue")
        # scale_y_reverse(limits=c(1, 0), expand=c(0, 0), breaks=seq(1, 0, -0.2))

lm1.out <- lm(spr ~ structure*sp, data=fdat)
lm1.pred <- predict(lm1.out, structure=abi)
lm2.out <- lm(log(spr) ~ structure*sp, data=fdat)
lm2.pred <- exp(predict(lm2.out, structure=abi))
fdat %>% mutate(pred1=lm1.pred, pred2=lm2.pred) %>%
    ggplot(aes(group=sp))+
        geom_point(aes(x=structure, y=spr, color=factor(sp)))+
        geom_line(aes(x=structure, y=spr, color=factor(sp)))+
        geom_line(aes(x=structure, y=pred1), color="red")+
        geom_line(aes(x=structure, y=pred2), color="blue")

ssq <- function(par, datarange, data){
    # if(is.nan(par)) par <- 0
    a <- exp(par[1])
    b <- par[2]
    out <- exp(-a*datarange) + b*datarange
    return(sum((out-data)^2))
}


optim(c(-5,-0.1), ssq, datarange=seq(0.2, 2, 0.1), data=fdat %>% filter(sp == 0.5) %>% pull(F))




spr_plot <- ggplot(fdat, aes(x = factor(structure), y=spr, group=factor(sp))) +
    geom_point(aes(color=factor(sp)))+
    geom_line(aes(color=factor(sp)))+
    # geom_smooth(method="lm")+
    annotate("segment", x = 8.5, y = 0.1, xend = 1, yend = 0.1,
        arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
    annotate("segment", x = 9.5, y = 0.1, xend = 19, yend = 0.1,
        arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
    geom_text(data=data.frame(), aes(x=6.5, y=0.0825, label="More Young Fish", size=7))+
    geom_text(data=data.frame(), aes(x=11.5, y=0.0825, label="More Old Fish", size=7))+
    # geom_hline(yintercept=0.40, linetype="dashed")+
    # geom_smooth(method="lm")+
    scale_y_reverse(limits=c(1, 0), expand=c(0, 0), breaks=seq(1, 0, -0.2))+
    scale_x_discrete()+
    guides(size="none")+
    labs(x="ABI40%", y="Spawning Potential Ratio", color="Old-Age Spawning Proportion")+
    theme_bw()+
    theme(panel.grid.minor = element_blank())

f_plot <- ggplot(fdat, aes(x = factor(structure), y=F, group=factor(sp))) +
    geom_point(aes(color=factor(sp)))+
    geom_line(aes(color=factor(sp)))+
    annotate("segment", x = 8.5, y = 0.225, xend = 1, yend = 0.225,
        arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
    annotate("segment", x = 9.5, y = 0.225, xend = 19, yend = 0.225,
        arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
    geom_text(data=data.frame(), aes(x=6.5, y=0.23, label="More Young Fish", size=7))+
    geom_text(data=data.frame(), aes(x=11.5, y=0.23, label="More Old Fish", size=7))+
    # geom_hline(yintercept=0.40, linetype="dashed")+
    # geom_smooth(method="lm")+
    scale_y_continuous(limits=c(0, 0.25), expand=c(0, 0), breaks=seq(0, 0.25, 0.05),)+
    scale_x_discrete()+
    guides(size="none")+
    labs(x="ABI40%", y="Fishing Mortality Rate", color="Old-Age Spawning Proportion")+
    theme_bw()+
    theme(panel.grid.minor = element_blank())

fspr_plot <- ggplot(fdat, aes(x = F, y=spr, group=factor(sp))) +
    geom_point(aes(color=factor(sp)))+
    geom_line(aes(color=factor(sp)))+
    geom_text(data=fdat %>% filter(structure %in% c(0.2, 2.0), sigma_r == 0.4), aes(x=F+0.005, y=spr, label=structure, color=factor(sp)))+
    guides(color="none")+
    scale_y_reverse(limits=c(1, 0), expand=c(0, 0), breaks=seq(1, 0, -0.2))+
    labs(y="ABI40%", x="Fishing Mortality Rate", color="Old-Age Spawning Proportion")+
    theme_bw()+
    theme(panel.grid.minor = element_blank())

library(patchwork)
((f_plot + spr_plot)/fspr_plot) + plot_layout(guides="collect") & custom_theme+theme(legend.position = "bottom")
ggsave("~/Desktop/fspr.png", width=16, height=12, units="in")



sim_fs %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    filter(time == 100, sigma_r == 0.4) %>%
    group_by(steepness, structure, sp, sigma_r) %>%
    mutate(
        F = compute_F_for_abi(abi=structure, F40s[match(sp, c(0.3, 0.5, 0.7))], dp=dp_y)
    ) %>%
    ungroup() %>%
    filter(sp == 0.3, sigma_r == 0.4, sim == seeds[1]) %>%
    select(structure, F)

expand.grid(biomass=seq(0, 2, 0.01), structure=abis) %>%
    as_tibble() %>%
    left_join(
        sim_fs %>%
            rename(
                "time"="Var1",
                "age"="Var2",
                "sex"="Var3",
                "region"="Var4"
            ) %>%
            filter(time == 100, sigma_r == 0.4) %>%
            group_by(steepness, structure, sp, sigma_r) %>%
            mutate(
                Fmax = compute_F_for_abi(abi=structure, F40s[match(sp, c(0.3, 0.5, 0.7))], dp=dp_y)
            ) %>%
            ungroup() %>%
            filter(sp == 0.3, sigma_r == 0.4, sim == seeds[1]) %>%
            select(structure, Fmax),
        by="structure"
    ) %>%
    rowwise() %>%
    mutate(
        F = threshold_f(
            x=biomass,
            f_min = 0,
            f_max = Fmax,
            lrp = 0.05,
            urp = 1
        )
    ) %>%

    ggplot(aes(x=biomass, y=F, color=factor(structure), group=factor(structure)))+
        geom_line()+
        scale_x_continuous(expand=c(0,0))+
        theme_bw()+
        labs(x="SSB/B40", color="Target ABI", title="Harvest Strategies")+
        guides(color=guide_legend(nrow=3))+
        custom_theme+
        theme(plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"))

ggsave("~/Desktop/hcr.png", width=12, height=12, units="in")    







ssb_data <- sim_naas %>% as_tibble() %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    # select(time, age, sex, region, sim, steepness, fishing, scale, structure, L1, value) %>%
    mutate(sex = case_when(sex == 1 ~ "fem", TRUE ~ "mal")) %>%
    dplyr::filter(time <= nyears) %>%
    left_join(reshape2::melt(dp$waa, value.name="waa"), by=c("time", "age", "sex", "region")) %>%
    left_join(reshape2::melt(dp$mat, value.name="mat"), by=c("time", "age", "sex", "region", "sp")) %>%
    select(time, age, sex, region, sim, steepness, structure, sp, om, sigma_r, L1, value, waa, mat) %>%
    mutate(
        ssbaa = value*waa*mat
    ) %>%
    select(-c(L1, value, waa, mat)) %>%
    group_by(time, sim, steepness, structure, sp, sigma_r, om) %>%
    summarise(ssb=sum(ssbaa)) %>%
    mutate(
        steepness = factor(steepness),
        structure = factor(structure),
        sp = factor(sp),
        om = factor(om, levels=c(0, 1, 2), labels=c("Equilibrium", "Recruit_Crash", "Mortality_Event"))
    )

ssb_data %>% write_csv("data/ssb_data_full2.csv")

# # SSB Distribution: distribution of SSB when at equilibrium
# ssb_dist <- ssb_data %>%
#     group_by(sp) %>%
#     mutate(
#         Bref = median(ssb[which(time > 100 & structure == 1)]),
#         depletion = ssb/Bref
#     ) %>%
#     filter(time > 100, om == "Equilibrium")

    
# ggplot(ssb_dist, aes(y=structure, x=depletion, fill=structure))+
#     stat_halfeye(aes(fill = after_stat(level)), .width = c(.50, 0.80, .95, 1))+
#     scale_fill_brewer(na.translate = FALSE) +
#     facet_wrap(~sp)+
#     ggtitle("SSB Distribution at Equilibrium")

# # SSB Timeseries: timeseries plots of spawning biomass
# ssb_tseries <- ssb_data %>% as_tibble() %>%
#     group_by(steepness, sigma_r, sp) %>%
#     mutate(Bref = median(ssb[which(time > 20 & time < 100 & structure == 1)])) %>%
#     group_by(time, steepness, fishing, scale, structure, sp, sigma_r, om) %>%
#     median_qi(ssb, Bref, .width=c(0.50))

# ggplot(ssb_tseries)+
#     geom_line(aes(x=time, y=ssb, color=structure, group=structure), linewidth=0.9)+
#     geom_hline(aes(yintercept=Bref), linetype="dashed")+
#     geom_hline(aes(yintercept=0.05*Bref), linetype="dashed")+
#     facet_grid(sp ~ sigma_r)+
#     coord_cartesian(ylim=c(0, 500), expand=0)+
#     ggtitle("Spawning Biomass Timeseries")+
#     theme_bw()


# # SSB Decline: depletion relative to pre-crash conditions 
# ssb_decline <- ssb_data %>% as_tibble() %>%
#     group_by(steepness, sigma_r, sp) %>%
#     mutate(Bref = median(ssb[which(time > 130 & structure == 1)])) %>%
#     group_by(sim, structure, sigma_r, steepness, sp, om) %>%
#     mutate(
#         recfail_biomass = mean(ssb[which(time > 75 & time < 91)]),
#         # min_bio = min(ssb[which(time > 75 & time < 126)]),
#         is_nearmin = ssb < Bref,
#         time_nearmin = sum(is_nearmin)
#     ) %>%
#     mutate(
#         sp = factor(sp, levels=c(0.9, 0.7, 0.5, 0.3))
#     )

# ggplot(ssb_decline, aes(y=structure, x=time_nearmin, color=sp, group=sp))+
#     stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
#     # scale_fill_brewer(na.translate = FALSE) +
#     facet_grid(steepness ~ sigma_r, scales="free_x")+
#     ggtitle("SSB Decline")+
#     theme_bw()


# # Spawning Biomass Stability: average annual variation in SSB
# ssb_stability <- ssb_data %>%
#     group_by(sim, steepness, fishing, structure, scale, sp, sigma_r, om) %>%
#     summarise(ssb_aav = aav(ssb))

# ggplot(ssb_stability, aes(y=structure, x=ssb_aav, color=sp, group=sp))+
#     stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
#     # scale_fill_brewer(na.translate = FALSE) +
#     facet_grid(steepness ~ sigma_r)+
#     ggtitle("SSB AAV")+
#     theme_bw()


# ssb_data %>%
#     group_by(steepness, sigma_r, sp) %>%
#     mutate(
#         Bref = median(ssb[which(time > 125 & structure == 1)]),
#         below = ssb < 0.2*reference_points$B0
#     ) %>%
#     group_by(steepness, sigma_r, sp, structure, sim, om) %>%
#     summarise(
#         prob_below = sum(below)/max(time)
#     ) %>%
#     group_by(steepness, sigma_r, sp, structure, om) %>%
#     median_qi(prob_below)







# # Recruitment: timeseries plots of recruitment
# recs <- sim_naas %>% as_tibble() %>%
#     rename(
#         "time"="Var1",
#         "age"="Var2",
#         "sex"="Var3",
#         "region"="Var4"
#     ) %>%
#     mutate(
#         steepness = factor(steepness),
#         fishing = factor(fishing),
#         scale = factor(scale),
#         structure = factor(structure),
#         sp = factor(sp),
#         om = factor(om, levels=c(0, 1, 2), labels=c("Equilibrium", "Recruit_Crash", "Mortality_Event"))
#     ) %>%
#     mutate(sex = case_when(sex == 1 ~ "fem", TRUE ~ "mal")) %>%
#     filter(age == 1)

# avg_rec <- recs %>% group_by(sim, steepness, fishing, structure, scale, sp, om, sigma_r) %>%
#     summarise(mu = mean(value))

# ggplot(avg_rec, aes(y=structure, x=mu, color=structure, group=sp))+
#     stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
#     # scale_fill_brewer(na.translate = FALSE) +
#     facet_grid(sp ~ sigma_r)+
#     ggtitle("Recruitment")+
#     theme_bw()

# recs %>% as_tibble() %>%
#     group_by(steepness, sigma_r, sp) %>%
#     mutate(Rref = median(value[which(time > 50 & time < 100 & structure == 1)])) %>%
#     group_by(time, steepness, fishing, scale, structure, sp, sigma_r, om) %>%
#     median_qi(value, Rref, .width=c(0.50)) %>%

#     ggplot()+
#         geom_line(aes(x=time, y=value, color=structure, group=structure), size=0.9)+
#         geom_hline(aes(yintercept=Rref), linetype="dashed")+
#         # geom_hline(aes(yintercept=0.80*Cref), linetype="dotted")+
#         facet_grid(sp ~ sigma_r)+
#         coord_cartesian(ylim=c(0, 100))+
#         ggtitle("Recruitment Timeseries")+
#         theme_bw()

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################


rel_catch_value <- sim_caas %>% as_tibble() %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    mutate(sex = case_when(sex == 1 ~ "fem", TRUE ~ "mal")) %>%
    dplyr::filter(time <= nyears) %>%
    left_join(reshape2::melt(dp$waa/max(dp$waa), value.name="waa"), by=c("time", "age", "sex", "region")) %>%
    mutate(rel_value = value*waa) %>%
    group_by(time, sim, steepness, structure, sp, sigma_r, om) %>%
    summarise(rel_value=sum(rel_value)) %>%
    mutate(
        steepness = factor(steepness),
        structure = factor(structure),
        sp = factor(sp),
        om = factor(om, levels=c(0, 1, 2), labels=c("Equilibrium", "Recruit_Crash", "Mortality_Event"))
    )

rel_catch_value %>% write_csv("data/catch_value_data_full2.csv")




sim_caas %>% as_tibble() %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    filter(time <= nyears) %>%
    group_by(time, sim, steepness, structure, sp, sigma_r, om) %>%
    summarise(catch=sum(value)) %>%
    mutate(
        steepness = factor(steepness),
        structure = factor(structure),
        sp = factor(sp),
        om = factor(om, levels=c(0, 1, 2), labels=c("Equilibrium", "Recruit_Crash", "Mortality_Event"))
    )


# # Catch Timeseries: timeseries plots of spawning biomass
# catch_tseries <- catch_data %>% as_tibble() %>%
#     group_by(steepness, sigma_r, sp) %>%
#     mutate(Cref = median(catch[which(time > 130 & structure == 1)])) %>%
#     group_by(time, steepness, fishing, scale, structure, sp, sigma_r, om) %>%
#     median_qi(catch, Cref, .width=c(0.50))

# ggplot(catch_tseries %>% filter(sigma_r == 1.2))+
#     geom_line(aes(x=time, y=catch, color=structure, group=structure), size=0.9)+
#     geom_hline(aes(yintercept=Cref), linetype="dashed")+
#     geom_vline(xintercept=75)+
#     geom_vline(xintercept=90)+
#     geom_hline(aes(yintercept=0.05*Cref))+
#     facet_grid(sp ~ steepness)+
#     coord_cartesian(ylim=c(0, 100))+
#     ggtitle("Catch Timeseries")+
#     theme_bw()


# # SSB Decline: depletion relative to pre-crash conditions 
# catch_decline <- catch_data %>% as_tibble() %>%
#     group_by(steepness, sigma_r, sp) %>%
#     mutate(Cref = median(catch[which(time > 130 & structure == 1)])) %>%
#     group_by(sim, structure, steepness, sigma_r, sp, om) %>%
#     mutate(
#         is_nearmin = catch < 0.50*Cref,
#         time_nearmin = sum(is_nearmin)
#     ) %>%
#     mutate(
#         sp = factor(sp, levels=c(0.9, 0.7, 0.5, 0.3))
#     )

# ggplot(catch_decline, aes(y=structure, x=time_nearmin, color=sp, group=sp))+
#     stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
#     # scale_fill_brewer(na.translate = FALSE) +
#     facet_grid(steepness ~ sigma_r, scales="free_x")+
#     ggtitle("Catch Decline")+
#     theme_bw()


# # Spawning Biomass Stability: average annual variation in SSB
# catch_stability <- catch_data %>%
#     group_by(sim, steepness, fishing, structure, scale, sp, om, sigma_r) %>%
#     summarise(catch_aav = aav(catch))

# ggplot(catch_stability, aes(y=structure, x=catch_aav, color=sp, group=sp))+
#     stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
#     # scale_fill_brewer(na.translate = FALSE) +
#     facet_grid(steepness ~ sigma_r)+
#     ggtitle("Catch AAV")+
#     theme_bw()




reference_age_structures <- sim_caas %>% as_tibble() %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    # select(time, age, sex, region, sim, steepness, fishing, scale, structure, L1, value) %>%
    mutate(sex = case_when(sex == 1 ~ "fem", TRUE ~ "mal")) %>%
    filter(time <= nyears) %>%
    filter(time > 50 & time < 174) %>%
    group_by(sigma_r, sp, age) %>%
    summarise(value = mean(value)) %>%
    summarise(ref_naa = list(c(value)))

reference_age_structures %>%
    filter(sp == 0.7) %>%
    pull(ref_naa) %>% lapply(mean)


catch_abi_data <- sim_caas %>% as_tibble() %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    # select(time, age, sex, region, sim, steepness, fishing, scale, structure, L1, value) %>%
    mutate(sex = case_when(sex == 1 ~ "fem", TRUE ~ "mal")) %>%
    filter(time <= nyears) %>%
    filter(time > 50 & time < 174) %>%
    left_join(reference_age_structures, by=c("sp", "sigma_r")) %>%
    group_by(time, steepness, structure, sp, sigma_r, sim, om) %>%
    summarise(
        abi=abi2(value, ref=unlist(ref_naa), threshold=0.9)
    ) %>%
    select(time, steepness, structure, sim, sp, abi, sigma_r, om) %>%
    mutate(
        steepness = factor(steepness),
        structure = factor(structure),
        sp = factor(sp),
        om = factor(om, levels=c(0, 1, 2), labels=c("Equilibrium", "Recruit_Crash", "Mortality_Event"))
    )

catch_abi_data %>% write_csv("data/catch_abi_data_full2.csv")












#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
reference_age_structures <- sim_naas %>% as_tibble() %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    # select(time, age, sex, region, sim, steepness, fishing, scale, structure, L1, value) %>%
    mutate(sex = case_when(sex == 1 ~ "fem", TRUE ~ "mal")) %>%
    filter(time <= nyears) %>%
    filter(time > 50 & time < 174) %>%
    group_by(sigma_r, sp, age) %>%
    summarise(value = mean(value)) %>%
    summarise(ref_naa = list(c(value)))

reference_age_structures %>%
    filter(sp == 0.7) %>%
    pull(ref_naa) %>% lapply(mean)


abi_data2 <- sim_naas %>% as_tibble() %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    # select(time, age, sex, region, sim, steepness, fishing, scale, structure, L1, value) %>%
    filter(time <= nyears) %>%
    filter(time > 50 & time < 174) %>%
    left_join(reference_age_structures, by=c("sp", "sigma_r")) %>%
    group_by(time, steepness, structure, sp, sigma_r, sim, om) %>%
    summarise(
        abi=abi2(value, ref=unlist(ref_naa), threshold=0.9)
    ) %>%
    select(time, steepness, structure, sim, sp, abi, sigma_r, om) %>%
    mutate(
        steepness = factor(steepness),
        structure = factor(structure),
        sp = factor(sp),
        om = factor(om, levels=c(0, 1, 2), labels=c("Equilibrium", "Recruit_Crash", "Mortality_Event"))
    )

abi_data2 %>% write_csv("data/abi_data_full2.csv")

# # ABI Stability: average annual variation in ABI
# abi_stability <- abi_data %>%
#     group_by(sim, steepness, fishing, structure, scale, sp, sigma_r, om) %>%
#     summarise(abi_aav = aav(abi))

# ggplot(abi_stability, aes(y=structure, x=abi_aav, color=sp, group=sp))+
#     stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
#     # scale_fill_brewer(na.translate = FALSE) +
#     facet_grid(steepness ~ sigma_r, scales="free_x")+
#     ggtitle("ABI AAV")+
#     theme_bw()

# # ABI Distribution: distribution of ABI levels when
# # stock is near equilibrium biomass levels
# abi_dist <- abi_data %>%
#     filter(time > 130, om == "Equilibrium")

    
# ggplot(abi_dist, aes(y=structure, x=abi, color=sp, group=sp))+
#     stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
#     # scale_fill_brewer(na.translate = FALSE) +
#     facet_grid(steepness ~ sigma_r, scales="free_x")+
#     ggtitle("ABI Distribution")+
#     theme_bw()


popdiversity_data <- sim_naas %>% as_tibble() %>%
    rename(
        "time"="Var1",
        "age"="Var2",
        "sex"="Var3",
        "region"="Var4"
    ) %>%
    # select(time, age, sex, region, sim, steepness, fishing, scale, structure, L1, value) %>%
    filter(time <= nyears) %>%
    filter(time > 50, time < 175) %>%
    left_join(reference_age_structures, by=c("sp", "sigma_r")) %>%
    group_by(time, steepness, structure, sp, sigma_r, sim, om) %>%
    summarise(
        div=shannon_diversity(value)
    ) %>%
    select(time, steepness, structure, sim, sp, div, sigma_r, om) %>%
    mutate(
        steepness = factor(steepness),
        structure = factor(structure),
        sp = factor(sp),
        om = factor(om, levels=c(0, 1, 2), labels=c("Equilibrium", "Recruit_Crash", "Mortality_Event"))
    )

popdiversity_data %>% write_csv("data/popdiv_data_full2.csv")

# avg_abi <- abi_dist %>% 
#     filter(sp == 0.7) %>%
#     ungroup() %>% 
#     group_by(structure) %>%
#     reframe(
#         mu = mean(abi)
#     ) %>%
#     column_to_rownames("structure")

# hc <- hclust(dist(avg_abi), method="average")
# hcdata <- dendro_data(hc)

# ggplot() + 
#     geom_segment(data = hcdata$segments, 
#                     aes(x = x, y = y, xend = xend, yend = yend)) + 
#     geom_text(data = hcdata$labels, 
#                 aes(x = x, y = y, label = label), size = 5, vjust = 3)+
#     theme_bw()+
#     theme(
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.line.x = element_blank(),
#         panel.border = element_blank(),
#         panel.grid = element_blank()
#     )




# #######################################################################
# #######################################################################
# #######################################################################
# #######################################################################
# #######################################################################
# #######################################################################

# sim_data1 <- ssb_data %>%
#     group_by(sim, sigma_r, steepness, sp, structure) %>%
#     summarise(ssb = mean(ssb)) %>%
#     left_join(
#         ssb_decline %>%
#             group_by(sim, sigma_r, steepness, sp, structure) %>%
#             summarise(collapsed_time = mean(time_nearmin)),
#         by = c("sim", "sp", "structure", "sigma_r", "steepness")
#     ) %>%
#     left_join(
#         ssb_stability %>%
#             group_by(sim, sigma_r, steepness, sp, structure) %>%
#             summarise(ssb_aav = mean(ssb_aav)),
#         by = c("sim", "sp", "structure", "sigma_r", "steepness")
#     ) %>%
#     left_join(
#         abi_stability %>%
#             group_by(sim, sigma_r, steepness, sp, structure) %>%
#             summarise(abi_aav = mean(abi_aav)),
#         by = c("sim", "sp", "structure", "sigma_r", "steepness")
#     ) %>%
#     left_join(
#         catch_stability %>%
#             group_by(sim, sp, structure) %>%
#             summarise(catch_aav = mean(catch_aav)),
#         by = c("sim", "sp", "structure")
#     ) %>%
#     left_join(
#         catch_data %>%
#             group_by(sim, sp, structure) %>%
#             summarise(catch = mean(catch)),
#         by = c("sim", "sp", "structure")
#     )
 
# s <- sim_data1 %>%
#     select(-collapsed_time) %>%
#     group_by(sigma_r, steepness, sp, sim) %>%
#     mutate(
#         across(ssb:catch, ~ (. - .[which(structure == 1)])/ .[which(structure == 1)])
#     ) %>%
#     mutate_at(vars(ends_with("aav")), ~ -1*.) %>%
#     mutate(util = ssb+ssb_aav+abi_aav+catch_aav+catch) %>%
#     group_by(sigma_r, steepness, sp, structure) %>%
#     median_qi(ssb, ssb_aav, abi_aav, catch, catch_aav, util, .width=c(0.50, 0.80, 0.95)) %>%  
#     reformat_ggdist_long(n=4)

# ggplot(s %>% filter(steepness == 0.75), aes(x=median, xmin=lower, xmax=upper, y=structure, color=sp, group=steepness, shape=steepness))+
#     geom_pointrange(position=position_dodge(width = 1))+
#     facet_grid(sigma_r ~ name, scales="free_x")

# sim_data_agg <- sim_data1 %>% group_by(structure, sp) %>%
#     summarise(across(everything(), ~ mean(.))) %>%
#     filter(sp == 0.3) %>%
#     column_to_rownames("structure")

# hc <- hclust(dist(sim_data_agg), method="average")
# hcdata <- dendro_data(hc)

# ggplot() + 
#     geom_segment(data = hcdata$segments, 
#                     aes(x = x, y = y, xend = xend, yend = yend)) + 
#     geom_text(data = hcdata$labels, 
#                 aes(x = x, y = y, label = label), size = 5, vjust=2)+
#     theme_bw()+
#     theme(
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.line.x = element_blank(),
#         panel.border = element_blank(),
#         panel.grid = element_blank()
#     )

# cluster_groups <- data.frame(structure = factor(abis), structure_groups = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5))



# sim_data_agg_grouped <- sim_data1 %>% 
#     left_join(
#         cluster_groups,
#         by = "structure"
#     ) %>%
#     mutate(
#         structure_groups = factor(
#             structure_groups,
#             labels = c("Very Young","Young","MSY","Older","Very Old")
#         )
#     ) %>%
#     group_by(sp, structure_groups) %>%
#     median_qi(ssb, collapsed_time, ssb_aav, abi_aav, catch_aav, catch, .width=c(0.50, 0.80, 0.95)) %>%
#     reformat_ggdist_long(n=2) %>%
#     mutate(
#         name = factor(
#             name,
#             levels = c("ssb", "catch", "collapsed_time", "ssb_aav", "catch_aav", "abi_aav")
#         )
#     )

#     ggplot(sim_data_agg_grouped, aes(x=median, xmin=lower, xmax=upper, y=structure_groups, color=structure_groups))+
#         geom_pointrange(position=position_dodge(width = 1))+
#         geom_vline(data = sim_data_agg_grouped %>% filter(structure_groups == "MSY"), aes(xintercept=median))+
#         facet_grid(rows=vars(sp), cols=vars(name), scales="free_x")+
#         theme_bw()

# sim_data_agg_grouped$name
 

# # K-Means Clustering to find possible breakpoints
# kmeans_data <- sim_data1 %>%
#     select(-collapsed_time) %>%
#     group_by(sigma_r, steepness, sp, sim) %>%
#     mutate(
#         across(ssb:catch, ~ (. - .[which(structure == 1)])/ .[which(structure == 1)])
#     ) %>%
#     mutate_at(vars(ends_with("aav")), ~ -1*.) %>%
#     mutate(util = ssb+ssb_aav+abi_aav+catch_aav+catch) %>% 
#     ungroup() %>% 
#     select(-c(structure, ssb, ssb_aav, abi_aav, catch_aav, catch)) %>%
#     nest(util) %>%
#     mutate(
#         model = map(data, kmeans, 5),    # iterate model over list column of input data
#         assignments = map(model, \(x) data.frame(x["cluster"]) %>% mutate(groupid = data.table::rleid(cluster)) %>% pull(groupid)),
#         # separation = map(model, \(x) x["betweenss"]/x["totss"])
#     ) %>% unnest(assignments) %>% unnest(assignments) %>% mutate(structure = factor(rep(abis, length.out=10260)))

# sim_data1 %>% ungroup() %>% 
#     arrange(sim, sigma_r, steepness, sp, structure) %>% 
#     left_join(
#         kmeans_data %>% select(-c(data, model)),
#         by = c("sim", "sigma_r", "steepness", "sp", "structure")
#     ) %>% 
#     filter(sigma_r == 0.4) %>%

#     ggplot(aes(x=ssb, y=catch, color=assignments))+
#         geom_point()+
#         stat_ellipse(aes(x=ssb, y=catch,color=assignments, group=assignments),type = "norm", level=0.95)+
#         stat_ellipse(level = 0.0001, geom = "point", color = "black")+
#         facet_grid(steepness ~ sp)
#         # geom_polygon(data=hull, alpha=0.1)
        
# sim_data1 %>% ungroup() %>% 
#     arrange(sim, sigma_r, steepness, sp, structure) %>% 
#     left_join(
#         kmeans_data %>% select(-c(data, model)),
#         by = c("sim", "sigma_r", "steepness", "sp", "structure")
#     ) %>%
#     group_by(sigma_r, sp, structure) %>%
#     median_qi(assignments) %>%

#     ggplot(aes(x=assignments, xmin=.lower, xmax=.upper, y=structure, color=structure)) +
#         geom_vline(xintercept = 1)+
#         geom_vline(xintercept = 2)+
#         geom_vline(xintercept = 3)+
#         geom_vline(xintercept = 4)+
#         geom_vline(xintercept = 5)+
#         geom_vline(xintercept = 1.5, linetype="dashed")+
#         geom_vline(xintercept = 2.5, linetype="dashed")+
#         geom_vline(xintercept = 3.5, linetype="dashed")+
#         geom_vline(xintercept = 4.5, linetype="dashed")+
#         geom_pointrange()+
#         facet_grid(sigma_r ~ sp)







# abi_struct_matrix <- abi_dist %>% 
#     ungroup() %>% 
#     filter(sp==0.3) %>%
#     group_by(structure) %>%
#     slice_sample(n=250)

# abi_80tiles <- quantile(abi_struct_matrix$abi[abi_struct_matrix$structure ==1], c(0.10, 0.90))




# percentiles <- abi_dist %>% group_by(sp) %>%
#     summarise(
#         perc10 = quantile(abi[structure == 1], 0.10),
#         perc90 = quantile(abi[structure == 1], 0.90)
#     ) %>%
#     mutate(sp = factor(sp, labels=c(1, 2, 3)))

# abi_dist <- abi_dist %>% left_join(percentiles, by="sp")

# ggplot(abi_dist, aes(y=structure, x=abi))+
#     stat_halfeye(
#         aes(
#             fill_ramp = after_stat(level),
#             fill = after_stat(x < abi_80tiles[2] & x > abi_80tiles[1])
#         ), 
#         position="dodge", 
#         .width = c(.50, .80, 0.95, 1)
#     )+
#     # scale_fill_brewer(na.translate = FALSE) +
#     ggtitle("ABI Distribution")+
#     facet_wrap(~sp)+
#     theme_bw()+
#     theme(legend.position = "bottom")

# inspect_after_stat(p) %>% as_tibble() %>% pull(group) %>% unique

# abi_dist %>% 
#     ungroup() %>% 
#     group_by(structure, sp) %>%
#     reframe(
#         mu = mean(abi),
#         sigma = sd(abi),
#         q = quantile(abi, c(0.10, 0.50, 0.90))
#     ) %>%
#     mutate(
#         qname = rep(rep(c("ten", "fifty" ,"ninety"), 3), 8)
#     ) %>%
#     pivot_wider(names_from=qname, values_from=q, values_fn = mean) %>%
#     mutate(
#         p = pnorm(mu, mean=mu[which(structure == 1)], sd=sigma[which(structure == 1)])
#         # central_mass = pnorm(ten[which(structure == 1)], mean=mu, sd=sigma, lower.tail = FALSE) - pnorm(ninety[which(structure == 1)], mean=mu, sd=sigma, lower.tail = FALSE) 
#     ) %>%
#     select(structure, sp, p) %>%
#     pivot_wider(names_from=sp, values_from=p)



# # smallest <- which(ungrouped == min(abs(1-ungrouped)))[1]

# # ungrouped <- ungrouped[-smallest]
# # ungrouped

# sp_groups <- list()
# for(s in sp){
#     print(s)
#     ungrouped <- abis
#     groups <- list()
#     while(length(ungrouped) > 0){
#         smallest <- which(abs(1-ungrouped) == min(abs(1-ungrouped)))[1]
#         struct <- ungrouped[smallest]

#         group <- abi_dist %>% 
#             filter(sp == s, structure %in% ungrouped) %>%
#             ungroup() %>% 
#             group_by(structure, sp) %>%
#             reframe(
#                 mu = mean(abi),
#                 sigma = sd(abi)
#             ) %>%
#             mutate(
#                 p = pnorm(mu, mean=mu[which(structure == struct)], sd=sigma[which(structure == struct)])
#             ) %>%
#             filter(p < 0.905 & p > 0.105) %>%
#             pull(structure)

#         groups <- c(groups, list(group))
#         ungrouped <- ungrouped[! ungrouped %in% group]

#     }
#     sp_groups <- c(sp_groups, list(groups))
# }