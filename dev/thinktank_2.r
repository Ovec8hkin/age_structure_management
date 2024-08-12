rm(list=ls())

library(devtools)
library(patchwork)
library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(GGally)

afscOM_dir <- "~/Desktop/Projects/afscOM"
devtools::load_all(afscOM_dir)
source("R/reference_points.R")
source("R/harvest_control_rules.R")
source("R/simulate_TAC.R")
source("R/age_structure_stats.R")
source("R/data_utils.R")
source("R/data_processing.R")
source("R/run_mse.R")
source("R/run_mse_multiple.R")


tier3 <- function(ref_pts, naa, dem_params){
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    return(
        npfmc_tier3_F(ssb, ref_pts$B40, ref_pts$F40)
    )
}
naa_scalar <- function(ref_pts, naa, dem_params, as_func, ref_naa, ...){
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    return(
        as_scalar_threshold_f(
            ssb/ref_pts$B40, 
            naa=naa[,,1,], 
            ref_naa=ref_naa,
            as_func = as_func,
            ...,
            f_min=0,
            f_max=ref_pts$F40,
            lrp=0.05,
            urp=1.0
        )
    )
}

ref_as_00 <- readRDS("data/agestruct_f00.RDS")
ref_as_40 <- readRDS("data/agestruct_f40.RDS")

sable_om <- readRDS("data/sablefish_om.RDS")
sable_om$model_options$simulate_observations <- TRUE # Enable simulating observations
# Generated 50 age composition samples from the LL survey and fixed gear fishery
sable_om$model_options$obs_pars$surv_ll$ac_samps <- 100
sable_om$model_options$obs_pars$surv_tw$ac_samps <- 100
sable_om$model_options$obs_pars$fish_fx$ac_samps <- 100
# Generate age comp samples as integers rather than proportions
# (this is necesarry for the SpatialSablefishAssessment TMB model)
sable_om$model_options$obs_pars$surv_ll$as_integers = TRUE
sable_om$model_options$obs_pars$surv_tw$as_integers = TRUE
sable_om$model_options$obs_pars$fish_fx$as_integers = TRUE
# Decrease observation error
sable_om$model_options$obs_pars$surv_ll$rpn_cv <- 0.20
sable_om$model_options$obs_pars$surv_ll$rpw_cv <- 0.10
sable_om$model_options$obs_pars$surv_tw$rpw_cv <- 0.10
# Turn on estimation model
sable_om$model_options$run_estimation = FALSE

nsims=20
seeds <- sample(1:1e7, nsims)

om_tier3 <- run_mse_multiple(nsims=20, seeds=seeds, nyears=160, om=sable_om, hcr=tier3)
om_avgage <- run_mse_multiple(nsims=20, seeds=seeds, nyears=160, om=sable_om, hcr=naa_scalar, as_func=average_age, ages=2:31)

models <- afscOM::listN(om_tier3, om_avgage)

#' Plot regressions between age structure metrics
#' and SSB, Recruitment, and F
#' ----------------------------------------------

as_metrics <- bind_mse_outputs(models, "naa", extra_columns=list(hcr=c("tier3", "avgage"))) %>%
    as_tibble() %>%
    drop_na() %>%
    left_join(
        bind_mse_outputs(models, "out_f", extra_columns=list(hcr=c("tier3", "avgage"))) %>%
            rename("summary_F"=value), 
        by=c("time", "hcr", "sim")) %>%
    filter(sex == "F", time < 65, sim==min(sim), hcr=="tier3") %>%
    group_by(time, sim) %>%
    summarise(
        diversity = shannon_diversity(value),
        avgage    = average_age(value, ages=2:31),
        pmat      = prop_mature(value, mat=sable_om$dem_params$mat[1,,1,]),
        pmat_full = prop_fully_mature(value, mat=sable_om$dem_params$mat[1,,1,]),
        abi0      = abi(value, ref=ref_as_00$naa),
        ssb       = sum(value*sable_om$dem_params$waa[1,,1,]*sable_om$dem_params$mat[1,,1,]),
        rec       = value[age == 2],
        f         = mean(summary_F)
    ) %>%
    pivot_longer(-c(time, sim, ssb, rec, f), names_to="metric", values_to="value") %>%
    mutate(
        metric = factor(metric, levels=c("diversity", "avgage", "pmat", "pmat_full", "abi0"), labels=c("Shannon Diversity", "Average Age", "Proportion Mature", "Prop Fully Mature", "ABI"))
    )

rsq <- as_metrics %>% 
    split(~ metric) %>% 
    map(\(df) lm(value ~ ssb, data=df)) %>% 
    map(summary) %>%
    map(\(x) x$r.sq) %>%
    as.data.frame %>%
    pivot_longer(everything(), names_to="metric", values_to="rsq") %>%
    mutate(
        y = c(1.75, 3.875, 0.125, 0.03125, 0.125),
        metric = factor(metric, levels=c("Shannon.Diversity", "Average.Age", "Proportion.Mature", "Prop.Fully.Mature", "ABI"), labels=c("Shannon Diversity", "Average Age", "Proportion Mature", "Prop Fully Mature", "ABI"))
    ) %>% 
    left_join(
        as_metrics %>% 
        split(~ metric) %>% 
        map(\(df) lm(value ~ rec, data=df)) %>% 
        map(summary) %>%
        # map(\(x) coef(x)["ssb", "Estimate"])%>%
        # map(\(x) coef(x)["(Intercept)", "Estimate"])%>%
        map(\(x) x$r.sq) %>%
        as.data.frame %>%
        pivot_longer(everything(), names_to="metric", values_to="rsq") %>%
        mutate(
            y = c(3.25, 9, 0.875, 0.225, 0.875),
            metric = factor(metric, levels=c("Shannon.Diversity", "Average.Age", "Proportion.Mature", "Prop.Fully.Mature", "ABI"), labels=c("Shannon Diversity", "Average Age", "Proportion Mature", "Prop Fully Mature", "ABI"))
        ),
        by=c("metric")
    ) %>% 
    left_join(
        as_metrics %>% 
        split(~ metric) %>% 
        map(\(df) lm(f ~ value, data=df)) %>% 
        map(summary) %>%
        # map(\(x) coef(x)["ssb", "Estimate"])%>%
        # map(\(x) coef(x)["(Intercept)", "Estimate"])%>%
        map(\(x) x$r.sq) %>%
        as.data.frame %>%
        pivot_longer(everything(), names_to="metric", values_to="rsq") %>%
        mutate(
            y = c(1.75, 3.875, 0.125, 0.03125, 0.125),
            metric = factor(metric, levels=c("Shannon.Diversity", "Average.Age", "Proportion.Mature", "Prop.Fully.Mature", "ABI"), labels=c("Shannon Diversity", "Average Age", "Proportion Mature", "Prop Fully Mature", "ABI"))
        ),
        by=c("metric")
    )

#' Correlation scatterplots between AS metrics and SSB
ggplot(as_metrics, aes(x=ssb, y=value, color=metric))+
    geom_point()+
    geom_smooth(aes(x=ssb, y=value), color="black", method='lm', formula= y~x, fullrange=TRUE)+
    geom_text(data=rsq, aes(x=200, y=y.x, label=paste0(expression(R^2), "=", round(rsq.x, 3))), vjust=0, hjust=0)+
    facet_wrap(~metric, scales="free_y")+
    facetted_pos_scales(
        y=list(
            scale_y_continuous(limits=c(1.5, 3.5)),
            scale_y_continuous(limits=c(3, 10)),
            scale_y_continuous(limits=c(0, 1)),
            scale_y_continuous(limits=c(0, 0.25)),
            scale_y_continuous(limits=c(0, 1))
        )
    )+
    scale_x_continuous(limits=c(0, 300))+
    labs(x="SSB", y="Metric Value")+
    coord_cartesian(expand=0)+
    guides(color="none")+
    theme_bw()
ggsave("~/Desktop/think_tank/asm_ssb_scatter.png", bg="transparent", width=10, unit="in")

ggplot(as_metrics, aes(x=2*rec, y=value, color=metric))+
    geom_point()+
    geom_smooth(aes(x=2*rec, y=value), color="black", method='lm', formula= y~x, fullrange=TRUE)+
    geom_text(data=rsq, aes(x=60, y=y.y, label=paste0(expression(R^2), "=", round(rsq.y, 3))), vjust=0, hjust=0)+
    facet_wrap(~metric, scales="free_y")+
    facetted_pos_scales(
        y=list(
            scale_y_continuous(limits=c(1.5, 3.5)),
            scale_y_continuous(limits=c(3, 10)),
            scale_y_continuous(limits=c(0, 1)),
            scale_y_continuous(limits=c(0, 0.25)),
            scale_y_continuous(limits=c(0, 1))
        )
    )+
    scale_x_continuous(limits=c(0, 100))+
    labs(x="Age-2 Recruitment", y="Metric Value")+
    coord_cartesian(expand=0)+
    guides(color="none")+
    theme_bw()
ggsave("~/Desktop/think_tank/asm_rec_scatter.png", bg="transparent", width=10, unit="in")

ggplot(as_metrics, aes(x=value, y=f, color=metric))+
    geom_point()+
    geom_smooth(aes(x=value, y=f), color="black", method='lm', formula= y~x, fullrange=TRUE)+
    geom_text(data=rsq, aes(x=y, y=y, label=paste0(expression(R^2), "=", round(rsq, 3))), vjust=0, hjust=0)+
    facet_wrap(~metric, scales="free_x")+
    facetted_pos_scales(
        x=list(
            scale_x_continuous(limits=c(1.5, 3.5)),
            scale_x_continuous(limits=c(3, 10)),
            scale_x_continuous(limits=c(0, 1)),
            scale_x_continuous(limits=c(0, 0.25)),
            scale_x_continuous(limits=c(0, 1))
        )
    )+
    scale_y_continuous(limits=c(0, 0.15))+
    labs(x="F", y="Metric Value")+
    coord_cartesian(expand=0)+
    guides(color="none")+
    theme_bw()
ggsave("~/Desktop/think_tank/asm_f_scatter.png", bg="transparent", width=10, unit="in")

#' Correlation scatterplots between age structure metrics
as_metrics %>% 
    ungroup() %>%
    pivot_wider(names_from="metric", values_from="value") %>%
    select(`Shannon Diversity`, `Average Age`, `Proportion Mature`, `Prop Fully Mature`, `ABI`) %>%
    ggpairs(diag = list("continuous"="blankDiag"))+
    theme_bw()
ggsave("~/Desktop/think_tank/as_corr_scatter.png", width=10, unit="in")

#'
#' 

bind_mse_outputs(models, "naa", extra_columns=list(hcr=c("tier3", "avgage"))) %>%
    as_tibble() %>%
    left_join(melt(sable_om$dem_params$mat, value.name="mat"), by=c("time", "age", "sex")) %>%
    left_join(melt(sable_om$dem_params$waa, value.name="waa"), by=c("time", "age", "sex")) %>%
    mutate(
        bio = value*waa,
        ssb = value*waa*mat
    ) %>%
    select(-c("region.y", "region.x", "mat", "waa")) %>%
    filter(sex == "F", time < 65, sim==min(sim), hcr=="tier3") %>%
    group_by(time, sim) %>%
    summarise(
        diversity_naa = shannon_diversity(value),
        avgage_naa    = average_age(value, ages=2:31),
        pmat_naa      = prop_mature(value, mat=sable_om$dem_params$mat[1,,1,]),
        pmatfull_naa = prop_fully_mature(value, mat=sable_om$dem_params$mat[1,,1,]),
        abi0_naa      = abi(value, ref=ref_as_00$naa),
        diversity_bio = shannon_diversity(bio),
        avgage_bio    = average_age(bio, ages=2:31),
        pmat_bio      = prop_mature(bio, mat=sable_om$dem_params$mat[1,,1,]),
        pmatfull_bio = prop_fully_mature(bio, mat=sable_om$dem_params$mat[1,,1,]),
        abi0_bio      = abi(bio, ref=ref_as_00$bio),
        diversity_ssb = shannon_diversity(ssb),
        avgage_ssb    = average_age(ssb, ages=2:31),
        pmat_ssb      = prop_mature(ssb, mat=sable_om$dem_params$mat[1,,1,]),
        pmatfull_ssb = prop_fully_mature(ssb, mat=sable_om$dem_params$mat[1,,1,]),
        abi0_ssb      = abi(ssb, ref=ref_as_00$ssb),
        rec           = value[age == 2]
    ) %>%
    pivot_longer(-c(time, sim, rec), names_to="metric", values_to="value") %>%
    separate(metric, into=c("metric", "quantity"), sep="_") %>%
    mutate(
        metric = factor(metric, levels=c("diversity", "avgage", "pmat", "pmatfull", "abi0"), labels=c("Shannon Diversity", "Average Age", "Proportion Mature", "Prop Fully Mature", "ABI")),
        quantity = factor(quantity, levels=c("naa", "bio", "ssb"), labels=c("Numbers -at-age", "Biomass -at-age", "SSB -at-age"))
    ) %>%

    # ggplot()+
    #     geom_line(aes(x=time, y=value, color=metric, group=quantity, linetype=quantity))+
    #     facet_wrap(~metric, scales="free_y")+
    #     facetted_pos_scales(
    #         y=list(
    #             scale_y_continuous(limits=c(0, 3.5)),
    #             scale_y_continuous(limits=c(0, 20)),
    #             scale_y_continuous(limits=c(0, 1)),
    #             scale_y_continuous(limits=c(0, 0.6)),
    #             scale_y_continuous(limits=c(0, 1))
    #         )
    #     )+
    #     theme_bw()

    ggplot()+
        geom_line(aes(x=time, y=value, color=metric))+
        geom_line(aes(x=time, y=rec/100))+
        facet_grid(vars(metric), vars(quantity), scales="free_y")+
        facetted_pos_scales(
            y=list(
                scale_y_continuous(limits=c(0, 3.5), sec.axis=sec_axis(~ .*10)),
                scale_y_continuous(limits=c(0, 20)),
                scale_y_continuous(limits=c(0, 1)),
                scale_y_continuous(limits=c(0, 0.6)),
                scale_y_continuous(limits=c(0, 1))
            )
        )+
        theme_bw()











#' Plot SSB across model runs
#' --------------------------
bind_mse_outputs(models, c("naa"), extra_columns = list(hcr=c("tier3", "avgage"))) %>% as_tibble() %>%
    drop_na() %>%
    left_join(
        melt(sable_om$dem_params$waa, value.name="weight"), 
        by=c("time", "age", "sex")
    ) %>%
    left_join(
        melt(sable_om$dem_params$mat, value.name="maturity"), 
        by=c("time", "age", "sex")
    ) %>%
    drop_na() %>%
    mutate(
        biomass = value*weight,
        spbio = value*weight*maturity
    ) %>%
    filter(sex == "F") %>%
    group_by(time, hcr, sim) %>%
    summarise(spbio=sum(spbio)) %>%
    group_by(time, hcr) %>%
    median_qi(spbio, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
    reformat_ggdist_long(n=2) %>%

    ggplot() + 
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper)) +
        scale_fill_grey(start=0.9, end=0.60)+
        facet_wrap(~hcr)+
        scale_y_continuous(limits=c(0, 300))+
        theme_bw()




selectivity <- data.frame(
    age=2:31,
    fix_sel_f = as.vector(sable_om$dem_params$sel[64,,1,1,1,drop=FALSE]),
    fix_sel_m = as.vector(sable_om$dem_params$sel[64,,2,1,1,drop=FALSE]), 
    twl_sel_f = as.vector(sable_om$dem_params$sel[64,,1,1,2,drop=FALSE]),
    twl_sel_m = as.vector(sable_om$dem_params$sel[64,,2,1,2,drop=FALSE]), 
    joint_sel_f=as.vector(apply(sable_om$dem_params$sel[64,,1,,,drop=FALSE], c(1, 2), sum)/max(apply(sable_om$dem_params$sel[64,,1,,,drop=FALSE], c(1, 2), sum))),
    joint_sel_m=as.vector(apply(sable_om$dem_params$sel[64,,2,,,drop=FALSE], c(1, 2), sum)/max(apply(sable_om$dem_params$sel[64,,2,,,drop=FALSE], c(1, 2), sum)))
) %>% as_tibble() %>%
    pivot_longer(-c(age), names_to="q", values_to = "sel") %>%
    separate(q, into=c("fishery", "meh", "sex")) %>%
    mutate(
        fishery = factor(fishery, levels=c("joint", "fix", "twl"), labels=c("Joint", "Longline", "Trawl")),
        sex = factor(sex, labels=c("Female", "Male"))
    )

ggplot(selectivity)+
    geom_line(aes(x=age, y=sel, color=sex, group=sex), size=1)+
    facet_wrap(~fishery, nrow=1)+
    scale_y_continuous(limits=c(0, 1.1), expand=c(0, 0))+
    scale_x_continuous(breaks=seq(0, 30, 5))+
    theme_bw()
ggsave("~/Desktop/think_tank/selectivity.png")


bind_mse_outputs(models, "naa", extra_columns=list(hcr=c("tier3", "avgage"))) %>%
    as_tibble() %>%
    drop_na() %>%
    filter(sex == "F", time < 65, sim==min(sim), hcr=="tier3", age < 12 | age > 21) %>%

    ggplot()+
        geom_line(aes(x=time, y=value))+
        scale_x_continuous(breaks=seq(0, 60, 20), labels=seq(1960, 2020, 20))+
        labs(x="Year", y="Number of individuals (millions)")+
        facet_wrap(~age, scales="free_y")+
        theme_bw()


bind_mse_outputs(models, "naa", extra_columns = list(hcr=c("tier3", "avgage"))) %>%
    as_tibble() %>%
    filter(sex == "F" & time < 65) %>%
    drop_na() %>%
    mutate(
        # class = factor(
        #     case_when(age < 3 ~ "1/2", age < 5 ~ "2/3", age < 7 ~ "3/4", age < 9 ~ "4/5", age < 15 ~ "5/7", age > 14 ~ "7+"), 
        #     levels=c("1/2", "2/3", "3/4", "4/5", "5/7", "7+"), 
        #     labels=c("Grade 1/2 (1-2yo)", "Grade 2/3 (3-4yo)", "Grade 3/4 (5-6yo)", "Grade 4/5 (7-8yo)", "Grade 5/7 (9-14yo)", "Grade 7+ (15+yo)")
        # )
        class = factor(
            case_when(age < 7 ~ "Young", age < 16 ~ "Immature", age > 15 ~ "Mature"), 
            levels=c("Young", "Immature", "Mature"), 
            labels=c("Young (<7 yo)", "Immature (7-15 yo)", "Mature (16+ yo)")
        )
    ) %>%
    group_by(time, class, sim, hcr) %>%
    summarise(
        total_num = sum(value)
    ) %>%
    group_by(time, class, hcr) %>%
    summarise(
        total_num = mean(total_num)
    ) %>%
    #filter(ref == "00", hcr == "tier3") %>%

    ggplot()+
        geom_col(aes(x=time, y=total_num, fill=class))+
        #scale_fill_viridis(direction=-1, discrete=TRUE, option="magma")+
        geom_hline(yintercept = 45)+
        #facet_wrap(~hcr)+
        #facet_grid(rows=vars(ref), cols=vars(hcr))+
        scale_x_continuous(breaks=seq(1, nyears+1, 20), labels=seq(1960, 1960+nyears+1, 20))+
        coord_cartesian(expand=0)+
        labs(x="Year", title="Numbers-at-Age", fill="Age Group")+
        guides(fill=guide_legend(reverse=TRUE))+
        theme_bw()+
        theme(
            axis.text = element_text(size=12),
            axis.title.y=element_blank(), 
            legend.position = "bottom"
        )
ggsave("~/Desktop/think_tank/sablefish_boffff.png")
