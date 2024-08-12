library(devtools)
library(patchwork)
library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)

afscOM_dir <- "~/Desktop/Projects/afscOM"
devtools::load_all(afscOM_dir)
source("R/reference_points.R")
source("R/harvest_control_rules.R")
source("R/simulate_TAC.R")
source("R/age_structure_stats.R")
source("R/data_utils.R")
source("R/data_processing.R")


tier3 <- function(ref_pts, naa, dem_params){
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    return(
        npfmc_tier3_F(ssb, ref_pts$B40, ref_pts$F40)
    )
}

shannon_naa_scalar <- function(ref_pts, naa, dem_params, ref_naa, ...){
    sbaa <- naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE]
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    return(
        as_scalar_threshold_f(
            ssb/ref_pts$B40, 
            naa=naa[,,1,], 
            ref_naa=ref_naa,
            as_func = shannon_diversity,
            f_min=0,
            f_max=ref_pts$F40,
            lrp=0.05,
            urp=1.0
        )
    )
}

avgage_naa_scalar <- function(ref_pts, naa, dem_params, ref_naa, ages){
    sbaa <- naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE]
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    return(
        as_scalar_threshold_f(
            ssb/ref_pts$B40, 
            naa=naa[,,1,], 
            ref_naa=ref_naa,
            as_func = average_age,
            ages=ages,
            f_min=0,
            f_max=ref_pts$F40,
            lrp=0.05,
            urp=1.0
        )
    )
}

abi_naa_scalar <- function(ref_pts, naa, dem_params, ref_naa, start_age){
    sbaa <- naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE]
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    return(
        as_scalar_threshold_f(
            ssb/ref_pts$B40, 
            naa=naa[,,1,], 
            ref_naa=ref_naa,
            as_func = abi,
            ref=ref_naa,
            start_age=start_age,
            f_min=0,
            f_max=ref_pts$F40,
            lrp=0.05,
            urp=1.0
        )
    )
}

ref_as_00 <- readRDS("data/agestruct_f00.RDS")
ref_as_40 <- readRDS("data/agestruct_f40.RDS")

nsims <- 10

sable_om <- readRDS("data/sablefish_om.RDS")

om_tier3_00     <- run_mse(sable_om, tier3, nsims=nsims)
om_shannon_00   <- run_mse(sable_om, shannon_naa_scalar, ref_naa = ref_as_00$naa, nsims=nsims)
om_avgage_00    <- run_mse(sable_om, avgage_naa_scalar, ref_naa = ref_as_00$naa, ages=2:31, nsims=nsims)
om_abi_00       <- run_mse(sable_om, abi_naa_scalar, ref_naa = ref_as_00$naa, start_age=3, nsims=nsims)
om_abi2_00      <- run_mse(sable_om, abi_naa_scalar, ref_naa = ref_as_00$naa, start_age=7, nsims=nsims)

om_tier3_40     <- run_mse(sable_om, tier3, nsims=nsims)
om_shannon_40   <- run_mse(sable_om, shannon_naa_scalar, ref_naa = ref_as_40$naa, nsims=nsims)
om_avgage_40    <- run_mse(sable_om, avgage_naa_scalar, ref_naa = ref_as_40$naa, ages=2:31, nsims=nsims)
om_abi_40       <- run_mse(sable_om, abi_naa_scalar, ref_naa = ref_as_40$naa, start_age=3, nsims=nsims)
om_abi2_40      <- run_mse(sable_om, abi_naa_scalar, ref_naa = ref_as_40$naa, start_age=7, nsims=nsims)

model_runs <- list(
    om_tier3_00, 
    om_shannon_00, 
    om_avgage_00, 
    om_abi_00, 
    om_abi2_00, 
    om_tier3_40, 
    om_shannon_40,
    om_avgage_40,
    om_abi_40,
    om_abi2_40
)
extra_columns <- list(
    hcr = rep(c("tier3", "shannon", "avgage", "abi", "abi2"), 2),
    ref = rep(c("00", "40"), each=5)
)


shannon_biomass_df <- create_summary_biomass_df(om_shannon$naa, dem_params$waa, dem_params$mat)
shannon_catch_df <- create_summary_catch_df(om_shannon$caa, om_shannon$out_f)
shannon_d <- create_biomass_catch_summary_df(shannon_biomass_df,shannon_catch_df)

avgage_biomass_df <- create_summary_biomass_df(om_avgage$naa, dem_params$waa, dem_params$mat)
avgage_catch_df <- create_summary_catch_df(om_avgage$caa, om_avgage$out_f)
avgage_d <- create_biomass_catch_summary_df(avgage_biomass_df,avgage_catch_df)

abi_biomass_df <- create_summary_biomass_df(om_abi$naa, dem_params$waa, dem_params$mat)
abi_catch_df <- create_summary_catch_df(om_abi$caa, om_abi$out_f)
abi_d <- create_biomass_catch_summary_df(biomass_df,catch_df)

d <- shannon_d %>% mutate(hcr = "shannon") %>%
    bind_rows(avgage_d %>% mutate(hcr="avgage")) %>%
    bind_rows(abi_d %>% mutate(hcr="abi"))

# Summary Plot
ggplot(d)+
    geom_line(aes(x=time, y=median, color=hcr), size=0.5)+
    geom_vline(xintercept=2023-1960+1, linetype="dashed")+
    geom_hline(data=ref_points_df, aes(yintercept=rp), linetype="dashed", color="red")+
    scale_fill_brewer(palette = "Blues")+
    scale_y_continuous(limits=c(0, 100))+
    facet_wrap(~name, scales="free_y")+
    facetted_pos_scales(y=pos_scales)+
    scale_x_continuous(breaks=seq(1, nyears+1, 20), labels=seq(1960, 1960+nyears+1, 20))+
    coord_cartesian(expand=0)+
    theme_bw()+
    theme(
        strip.background = element_blank(),
        strip.text = element_text(size=14, hjust=0),
        panel.spacing.x = unit(0.75, "cm"),
        axis.title = element_blank(),
        legend.position = "bottom"
    )

# Spawning Biomass Trajectories
bind_mse_outputs(model_runs, "naa", extra_columns) %>% as_tibble() %>%
    drop_na() %>%
    left_join(
        melt(dem_params$waa, value.name="weight"), 
        by=c("time", "age", "sex")
    ) %>%
    left_join(
        melt(dem_params$mat, value.name="maturity"), 
        by=c("time", "age", "sex")
    ) %>%
    drop_na() %>%
    mutate(
        biomass = value*weight,
        spbio = value*weight*maturity
    ) %>%
    filter(sex == "F") %>%
    group_by(time, sim, hcr, ref) %>%
    summarise(spbio=sum(spbio)) %>%
    group_by(time, hcr, ref) %>%
    median_qi(spbio, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
    reformat_ggdist_long(n=3) %>%

    ggplot() + 
        geom_line(aes(x=time, y=median, color=hcr)) +
        facet_wrap(~ref)+
        scale_y_continuous(limits=c(0, 300))+
        theme_bw()

# Total Catch Trajectories
bind_mse_outputs(model_runs, "caa", extra_columns) %>% as_tibble() %>%
    drop_na() %>%
    group_by(time, sims, hcr, ref) %>%
    summarise(catch=sum(value)) %>%
    group_by(time, hcr, ref) %>%
    median_qi(catch, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
    reformat_ggdist_long(n=3) %>%

    ggplot() + 
        geom_line(aes(x=time, y=median, color=hcr)) +
        facet_wrap(~ref)+
        scale_y_continuous(limits=c(0, 50))+
        theme_bw()

# Regressions
metrics <- bind_mse_outputs(model_runs, "naa", extra_columns) %>%
    as_tibble() %>%
    drop_na() %>%
    filter(sex == "F", time < 65, sim == 1, hcr=="tier3", ref=="00") %>%
    group_by(time) %>%
    summarise(
        diversity = shannon_diversity(value),
        avgage    = average_age(value, ages=2:31),
        pmat      = prop_mature(value, mat=dem_params$mat[1,,1,]),
        pmat_full = prop_fully_mature(value, mat=dem_params$mat[1,,1,]),
        abi0      = abi(value, ref=ref_as_00$naa),
        ssb       = sum(value*dem_params$waa[1,,1,]*dem_params$mat[1,,1,]),
        rec       = value[age == 2]
    ) %>%
    pivot_longer(-c(time, ssb, rec), names_to="metric", values_to="value") %>%
    mutate(
        metric = factor(metric, levels=c("diversity", "avgage", "pmat", "pmat_full", "abi0"), labels=c("Shannon Diversity", "Average Age", "Proportion Mature", "Prop Fully Mature", "ABI"))
    ) 

rsq <- metrics %>% 
    split(~ metric) %>% 
    map(\(df) lm(value ~ ssb, data=df)) %>% 
    map(summary) %>%
    # map(\(x) coef(x)["ssb", "Estimate"])%>%
    # map(\(x) coef(x)["(Intercept)", "Estimate"])%>%
    map(\(x) x$r.sq) %>%
    as.data.frame %>%
    pivot_longer(everything(), names_to="metric", values_to="rsq") %>%
    mutate(
        y = c(1.75, 3.875, 0.125, 0.03125, 0.125),
        metric = factor(metric, levels=c("Shannon.Diversity", "Average.Age", "Proportion.Mature", "Prop.Fully.Mature", "ABI"), labels=c("Shannon Diversity", "Average Age", "Proportion Mature", "Prop Fully Mature", "ABI"))
    ) %>% 
    left_join(
        metrics %>% 
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
    )

ggplot(metrics, aes(x=2*rec, y=value, color=metric))+
    geom_point()+
    geom_smooth(aes(x=2*rec, y=value), method='lm', formula= y~x)+
    geom_text(data=rsq, aes(x=70, y=y.y, label=paste0(expression(R^2), "=", round(rsq.y, 3))), vjust=0, hjust=0)+
    scale_x_continuous(limits=c(0, 100))+
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
    labs(x="Age-2 Recruitment", y="Metric Value")+
    coord_cartesian(expand=0)+
    guides(color="none")+
    theme_bw()

pdf(file="~/Desktop/as_ssb_acf.pdf")
par(mfrow=c(2, 3))
ms <- metrics %>% pull(metric) %>% unique
for(m in ms){
    met <- metrics %>% filter(metric == m)
    ccf(met$ssb, met$value, main=paste(m, "vs SSB"))
}
dev.off()

pdf(file="~/Desktop/as_rec_acf.pdf")
par(mfrow=c(2, 3))
ms <- metrics %>% pull(metric) %>% unique
for(m in ms){
    met <- metrics %>% filter(metric == m)
    ccf(met$rec, met$value, main=paste(m, "vs Recruitment"))
}
dev.off()

ggsave("~/Desktop/as_regressions_recruitment.pdf")
    # group_by(time, sim, hcr, ref) %>%
    # summarise(
    #     diversity = shannon_diversity(value),
    #     avgage    = average_age(value, ages=2:31),
    #     pmat      = prop_mature(value, mat=dem_params$mat[1,,1,]),
    #     pmat_full = prop_fully_mature(value, mat=dem_params$mat[1,,1,]),
    #     abi0      = abi(value, ref=ref_as_00$naa),
    #     avg_ssb   = mean(sum(value*dem_params$waa[1,,1,]*dem_params$mat[1,,1,]))
    # ) %>%
    # group_by(hcr, ref) %>%
    # median_qi(diversity, avgage, pmat, pmat_full, abi0, avg_ssb, .width=c(0.50, 0.85)) %>%
    # reformat_ggdist_long(n=2) %>%
    # mutate(
    #     name = recode_factor(name, !!!c("abi0"="ABI", "avgage"="Average Age", "diversity"="Shannon Diversity", "pmat"="Proportion Mature", "pmat_full"="Prop Fully Mature", "avg_ssb"="Average SSB")),
    #     hcr = recode_factor(hcr, !!!c("abi"="ABI", "abi2"="ABI2", "avgage"="Average Age", "shannon"="Shannon Diversity", "tier3"="Tier 3a"))
    # )

# Age Structure Performance Metrics
bind_mse_outputs(model_runs, "naa", extra_columns) %>%
    left_join(
        bind_mse_outputs(model_runs, "caa", extra_columns) %>% as_tibble() %>% filter(fleet == "Fixed"),
        by=c("time", "age" , "sex", "sim"="sims", "hcr", "ref")
    ) %>%
    as_tibble() %>%
    drop_na() %>%
    select(-c(L1.x, L1.y, starts_with("region"))) %>%
    rename("naa"="value.x", "caa"="value.y") %>%
    filter(time > 64) %>%
    group_by(time, sex, sim, hcr, ref) %>%
    summarise(
        diversity = shannon_diversity(naa),
        avgage    = average_age(naa, ages=2:31),
        # pmat      = prop_mature(naa[sex=="F"], mat=dem_params$mat[1,,1,]),
        # pmat_full = prop_fully_mature(naa[sex=="F"], mat=dem_params$mat[1,,1,]),
        abi0      = abi(naa, ref=ref_as_00$naa),
        avg_ssb   = mean(sum(naa*dem_params$waa[1,,1,]*dem_params$mat[1,,1,])),
        avg_catch = mean(sum(caa)),
        prop_big_catch = sum(caa[age > 7])/sum(caa)
    ) %>%
    group_by(hcr, ref) %>%
    median_qi(diversity, avgage, abi0, avg_ssb, avg_catch, prop_big_catch, .width=c(0.50, 0.85)) %>%
    reformat_ggdist_long(n=2) %>%
    mutate(
        name = recode_factor(name, !!!c("abi0"="ABI", "avgage"="Average Age", "diversity"="Shannon Diversity", "pmat"="Proportion Mature", "pmat_full"="Prop Fully Mature", "avg_ssb"="Average SSB")),
        hcr = recode_factor(hcr, !!!c("abi"="ABI", "abi2"="ABI2", "avgage"="Average Age", "shannon"="Shannon Diversity", "tier3"="Tier 3a"))
    ) %>%

    ggplot()+
        geom_pointinterval(aes(x=median, xmin=lower, xmax=upper, y=hcr, color=hcr))+
        geom_vline(
            data = . %>% filter(hcr == "Tier 3a" & .width == 0.5),
            aes(xintercept=median),
            linetype="dashed"
        )+
        #facet_wrap(~name, scales="free_x")+
        facet_grid(rows=vars(ref), cols=vars(name), scales="free_x")+
        facetted_pos_scales(
            list(
                scale_x_continuous(limits=c(0, 1.5)),
                scale_x_continuous(limits=c(5, 10)),
                scale_x_continuous(limits=c(2, 3)),
                scale_x_continuous(limits=c(0, 300)),
                scale_x_continuous(limits=c(0, 20)),
                scale_x_continuous(limits=c(0, 1))
            )
        )+
        theme_bw()+
        theme(legend.position = "none")

# Population NAA timeseries
p1 <- bind_mse_outputs(model_runs, "naa", extra_columns) %>%
    as_tibble() %>%
    filter(sex == "F" & time > 65) %>%
    drop_na() %>%
    mutate(
        class = factor(
            case_when(age < 3 ~ "1/2", age < 5 ~ "2/3", age < 7 ~ "3/4", age < 9 ~ "4/5", age < 15 ~ "5/7", age > 14 ~ "7+"), 
            levels=c("1/2", "2/3", "3/4", "4/5", "5/7", "7+"), 
            labels=c("Grade 1/2 (1-2yo)", "Grade 2/3 (3-4yo)", "Grade 3/4 (5-6yo)", "Grade 4/5 (7-8yo)", "Grade 5/7 (9-14yo)", "Grade 7+ (15+yo)")
        )
    ) %>%
    group_by(time, class, sim, hcr, ref) %>%
    summarise(
        total_num = sum(value)
    ) %>%
    group_by(time, class, hcr, ref) %>%
    summarise(
        total_num = mean(total_num)
    ) %>%
    #filter(ref == "00", hcr == "tier3") %>%

    ggplot()+
        geom_bar(aes(x=time, y=total_num, fill=forcats::fct_rev(class)), position="fill", stat="identity")+
        scale_fill_viridis(direction=-1, discrete=TRUE, option="magma")+
        # geom_hline(yintercept = 0.1)+
        geom_hline(yintercept = 0.5)+
        #facet_wrap(~hcr)+
        facet_grid(rows=vars(ref), cols=vars(hcr))+
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
p1
ggsave("~/Desktop/agestruct_timeseries_1960_2023.jpeg")


bind_mse_outputs(model_runs, "naa", extra_columns) %>% 
    as_tibble() %>%
    filter(sex == "F" & time > 64) %>%
    drop_na() %>%
    mutate(
        class = factor(
            case_when(age < 7 ~ "Young", age < 16 ~ "Immature", age > 15 ~ "Mature"), 
            levels=c("Young", "Immature", "Mature"), 
            labels=c("Young (<7yo)", "Immature (7-15yo)", "Mature (16+yo)")
        )
    ) %>%
    group_by(time, class, sim, hcr, ref) %>%
    summarise(
        total_num = sum(value)
    ) %>%
    group_by(time, class, hcr, ref) %>%
    summarise(
        total_num = mean(total_num)
    ) %>%
    group_by(time, hcr, ref) %>%
    summarise(
        prop = total_num/sum(total_num),
        class = class
    ) %>%

    ggplot()+
        geom_line(aes(x=time, y=prop, color=class))+
        #geom_hline(data=ref_as_df, aes(yintercept=prop, color=class), linetype="dashed")+
        facet_grid(rows=vars(ref), cols=vars(hcr))+
        scale_x_continuous(breaks=seq(1, nyears+1, 20), labels=seq(1960, 1960+nyears+1, 20))+
        scale_y_continuous(limits=c(0, 1))+
        coord_cartesian(expand=0)+
        theme_bw()




# Catch-at-age timeseries
p2 <- bind_mse_outputs(model_runs, "caa", extra_columns) %>% 
    as_tibble() %>%
    filter(sex == "F", fleet=="Fixed", time > 65) %>%
    drop_na() %>%
    mutate(
        class = factor(
            case_when(age < 3 ~ "1/2", age < 5 ~ "2/3", age < 7 ~ "3/4", age < 9 ~ "4/5", age < 15 ~ "5/7", age > 14 ~ "7+"), 
            levels=c("1/2", "2/3", "3/4", "4/5", "5/7", "7+"), 
            labels=c("Grade 1/2 (1-2yo)", "Grade 2/3 (3-4yo)", "Grade 3/4 (5-6yo)", "Grade 4/5 (7-8yo)", "Grade 5/7 (9-14yo)", "Grade 7+ (15+yo)")
        )
    ) %>%
    group_by(time, class, sims, hcr, ref) %>%
    summarise(
        total_num = sum(value)
    ) %>%
    group_by(time, class, hcr, ref) %>%
    summarise(
        total_num = mean(total_num)
    ) %>%
    #filter(ref == "00", hcr=="tier3") %>%

    ggplot()+
        geom_bar(aes(x=time, y=total_num, fill=forcats::fct_rev(class)), position="fill", stat="identity")+
        scale_fill_viridis(direction=-1, discrete=TRUE, option="magma")+
        # geom_hline(yintercept = 0.1)+
        geom_hline(yintercept = 0.5)+
        geom_hline(yintercept = 0.75)+
        #facet_wrap(~hcr)+
        facet_grid(rows=vars(ref), cols=vars(hcr))+
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
p2
ggsave("~/Desktop/caa_timeseries_1960_2023.jpeg")

library(patchwork)

(p1/p2)+plot_layout(guides="collect")
ggsave("~/Desktop/aa_timeseries.pdf", width=8.5)

bind_mse_outputs(model_runs, "caa", extra_columns) %>% 
    as_tibble() %>%
    filter(sex == "F", time > 64) %>%
    drop_na() %>%
    mutate(
        class = factor(
            case_when(age < 7 ~ "Young", age < 16 ~ "Immature", age > 15 ~ "Mature"), 
            levels=c("Young", "Immature", "Mature"), 
            labels=c("Young (<7yo)", "Immature (7-15yo)", "Mature (16+yo)")
        )
    ) %>%
    group_by(time, class, sims, hcr, ref) %>%
    summarise(
        total_num = sum(value)
    ) %>%
    group_by(time, class, hcr, ref) %>%
    summarise(
        total_num = mean(total_num)
    ) %>%
    group_by(time, hcr, ref) %>%
    summarise(
        prop = total_num/sum(total_num),
        class = class
    ) %>%

    ggplot()+
        geom_line(aes(x=time, y=prop, color=class))+
        facet_grid(rows=vars(ref), cols=vars(hcr))+
        scale_x_continuous(breaks=seq(1, nyears+1, 20), labels=seq(1960, 1960+nyears+1, 20))+
        scale_y_continuous(limits=c(0, 1))+
        coord_cartesian(expand=0)+
        theme_bw()


ref_as_df <- data.frame(
    ref = rep(c("00", "40"), each=3),
    class = rep(c("Young (<7yo)", "Immature (7-15yo)", "Mature (16+yo)"), 2),
    prop = c(
        sum(ref_as_00$naa[1:5])/sum(ref_as_00$naa),
        sum(ref_as_00$naa[6:15])/sum(ref_as_00$naa),
        sum(ref_as_00$naa[16:30])/sum(ref_as_00$naa),
        sum(ref_as_40$naa[1:5])/sum(ref_as_40$naa),
        sum(ref_as_40$naa[6:15])/sum(ref_as_40$naa),
        sum(ref_as_40$naa[16:30])/sum(ref_as_40$naa)
    )
)