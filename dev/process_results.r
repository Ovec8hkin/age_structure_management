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

custom_theme <- theme_bw()+
    theme(
        panel.spacing.y = unit(0.5, "cm"),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=14),
        axis.text = element_text(size=14),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.position = "bottom"
    )

##########################################################################
##########################################################################
##############            Population Resiliency             ##############
##########################################################################
##########################################################################

productivity <- read_csv("data/productivity_data_full2.csv", col_types="fffffd") %>% filter(steepness == 0.75)

ggplot(productivity, aes(x=recruitent, y=structure, color=structure, group=sp))+
    stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
    # scale_fill_brewer(na.translate = FALSE) +
    facet_grid(sp ~ sigma_r, scales="free_x")+
    ggtitle("Productivity")+
    theme_bw()

##########################################################################
##########################################################################
##############               Spawning Biomass               ##############
##########################################################################
##########################################################################
ssb_data <- read_csv("data/ssb_data_full2.csv", col_types="iffffffd") %>% filter(steepness == 0.75)

# SSB Distribution: distribution of SSB when at equilibrium
ssb_dist <- ssb_data %>%
    group_by(steepness, sigma_r, sp) %>%
    mutate(
        Bref = median(ssb[which(time > 50, time < 100, structure == 1)]),
        depletion = ssb/Bref
    ) %>%
    filter(time > 50, time < 174, om == "Equilibrium") %>%
    group_by(sigma_r, steepness, sp, structure) %>%
    median_qi(ssb, .width=c(0.50, 0.80, 0.95, 1))

ggplot(
    ssb_dist, aes(y=structure, x=ssb, xmin=.lower, xmax=.upper, color=structure, shape=sp, group=sp))+
    geom_pointinterval(position="dodge")+
    facet_grid(sp ~ sigma_r)+
    ggtitle("SSB Distribution at Equilibrium")+
    theme_bw()

# SSB Timeseries: timeseries plots of spawning biomass
ssb_tseries <- ssb_data %>% as_tibble() %>%
    group_by(steepness, sigma_r, sp) %>%
    mutate(Bref = median(ssb[which(time > 50 & time < 100 & structure == 1)])) %>%
    group_by(time, steepness, structure, sp, sigma_r, om) %>%
    median_qi(ssb, Bref, .width=c(0.50))

ggplot(ssb_tseries)+
    geom_line(aes(x=time, y=ssb, color=structure, group=structure), linewidth=0.9)+
    geom_hline(aes(yintercept=Bref), linetype="dashed")+
    facet_grid(sp ~ sigma_r)+
    coord_cartesian(ylim=c(0, 200), expand=0)+
    ggtitle("Spawning Biomass Timeseries")+
    guides(color=guide_legend(nrow = 2))+
    labs(x="Year", y="SSB", color="ABI")+
    theme_bw()+custom_theme
ggsave("~/Desktop/ssbtraj.png", width=14, height=12.5, units="in")


ssb_tseries1 <- ssb_data %>% as_tibble() %>%
    filter(sim == 140474, sigma_r == 0.4, sp == 0.5, structure == 1) %>% 
    mutate(Bref = median(ssb[which(time > 50 & time < 100 & structure == 1)])) %>%
    group_by(time, steepness, structure, sp, sigma_r, om) %>%
    median_qi(ssb, Bref, .width=c(0.50))

ggplot(ssb_tseries1)+
    geom_line(aes(x=time, y=ssb), color="black", linewidth=0.9)+
    geom_hline(aes(yintercept=Bref), linetype="dashed")+
    coord_cartesian(ylim=c(0, 75), xlim=c(0, 350), expand=0)+
    scale_x_continuous("Simulation Year", breaks=seq(0, 350, 50), labels=seq(0, 350, 50))+
    ggtitle("Example Simulation Trajectory")+
    custom_theme+
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm")
    )
ggsave("~/Desktop/simtraj.png", width=8.5, height=7, units="in")



# SSB Decline: depletion relative to pre-crash conditions 
rebound_timerate <- ssb_data %>% as_tibble() %>%
    group_by(steepness, sigma_r, sp, structure) %>%
    mutate(Bref = median(ssb[which(time > 50 & time < 100)])) %>%
    filter(time > 250) %>%
    mutate(
        min_ssb = min(ssb),
        depletion = ssb/Bref
    ) %>%
    group_by(sim, structure, sigma_r, steepness, sp, om) %>%
    mutate(
        # rate = lead(depletion) - depletion,
        rebuilt = depletion > 0.95,
        rebuild_time = min(time[which(rebuilt == TRUE)])-250,
        rebuild_rate = ((0.95*Bref)-min_ssb)/rebuild_time
    ) %>%
    select(time, sim, steepness, structure, sp, sigma_r, om, rebuild_time, rebuild_rate)
    
    
rebound_timerate_summ <- rebound_timerate %>%
    group_by(sigma_r, steepness, sp, structure) %>%
    median_qi(rebuild_rate, rebuild_time, .width=c(0.50, 0.80, 0.95))

ggplot(rebound_timerate_summ, aes(y=structure, x=rebuild_rate, xmin=rebuild_rate.lower, xmax=rebuild_rate.upper, color=structure, group=sp))+
    geom_pointinterval(position="dodge")+
    # scale_fill_brewer(na.translate = FALSE) +
    facet_grid(sp ~ sigma_r, scales="free_x")+
    ggtitle("SBS Rebuilt Rate")+
    theme_bw()

ggplot(rebound_timerate_summ, aes(y=structure, x=rebuild_time, xmin=rebuild_time.lower, xmax=rebuild_time.upper, color=structure, group=sp))+
    geom_pointinterval(position="dodge")+
    # scale_fill_brewer(na.translate = FALSE) +
    facet_grid(sp ~ sigma_r, scales="free_x")+
    ggtitle("SBS Rebuilt Time")+
    theme_bw()


# Spawning Biomass Stability: average annual variation in SSB
ssb_stability <- ssb_data %>%
    group_by(sim, steepness, structure, sp, sigma_r, om) %>%
    summarise(ssb_aav = aav(ssb))

ggplot(ssb_stability, aes(y=structure, x=ssb_aav, color=structure, group=sp))+
    stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
    # scale_fill_brewer(na.translate = FALSE) +
    facet_grid(sp ~ sigma_r)+
    ggtitle("SSB AAV")+
    theme_bw()



##########################################################################
##########################################################################
##############                Fishery Catch                 ##############
##########################################################################
##########################################################################
catch_data <- read_csv("data/catch_data_full2.csv", col_types="iffffffd") %>% filter(steepness == 0.75)
catch_abi_data <- read_csv("data/catch_abi_data_full2.csv", col_types="iffffdff") %>% filter(steepness == 0.75)
catch_value_data <- read_csv("data/catch_value_data_full2.csv", col_types="iffffffd") %>% filter(steepness == 0.75)

# Catch Distribution: distribution of catch when at equilibrium
catch_dist <- catch_data %>%
    group_by(steepness, sigma_r, sp) %>%
    mutate(
        Cref = median(catch[which(time > 50 & time < 100 & structure == 1)]),
        depletion = catch/Cref
    ) %>%
    filter(time > 50, time < 100, om == "Equilibrium") %>%
    group_by(sigma_r, steepness, sp, structure) %>%
    median_qi(catch, Cref, .width=c(0.50, 0.80, 0.95, 1))

ggplot(catch_dist, aes(y=structure, x=catch, xmin=catch.lower, xmax=catch.upper, color=structure))+
    geom_pointinterval()+
    geom_vline(aes(xintercept=Cref), linetype="dashed")+
    geom_vline(aes(xintercept=0.80*Cref), linetype="dotted")+
    facet_grid(sp ~ sigma_r)+
    ggtitle("Catch Distribution at Equilibrium")+
    theme_bw()

# Catch Timeseries: timeseries plots of spawning biomass
catch_tseries <- catch_data %>% as_tibble() %>%
    group_by(steepness, sigma_r, sp) %>%
    mutate(Cref = median(catch[which(time > 50 & time < 100 & structure == 1)])) %>%
    group_by(time, steepness, structure, sp, sigma_r, om) %>%
    median_qi(catch, Cref, .width=c(0.50))

ggplot(catch_tseries)+
    geom_line(aes(x=time, y=catch, color=structure, group=structure), size=0.9)+
    geom_hline(aes(yintercept=Cref), linetype="dashed")+
    geom_hline(aes(yintercept=0.80*Cref), linetype="dotted")+
    facet_grid(sp ~ sigma_r)+
    coord_cartesian(ylim=c(0, 25))+
    ggtitle("Catch Timeseries")+
    theme_bw()


# Catch Decline: depletion relative to pre-crash conditions 
catch_data %>% as_tibble() %>%
    group_by(steepness, sigma_r, sp, structure) %>%
    mutate(Cref = median(catch[which(time > 50 & time < 100)])) %>%
    filter(time > 110, time < 200) %>%
    mutate(
        depletion = catch/Cref
    ) %>%
    group_by(sim, structure, sigma_r, steepness, sp, om) %>%
    mutate(
        rate = lead(depletion) - depletion,
        rebuilt = depletion > 0.80
    ) %>%
    summarise(rate = mean(rate[which(rebuilt == FALSE)], na.rm=TRUE)) %>%

    ggplot(aes(y=structure, x=rate, color=structure))+
        stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
        facet_grid(sp ~ sigma_r)+
        ggtitle("Rate of Catch Increase")+
        theme_bw()

# Catch Stability: average annual variation in SSB
catch_stability <- catch_data %>%
    group_by(sim, steepness, structure, sp, om, sigma_r) %>%
    summarise(catch_aav = aav(catch))

ggplot(catch_stability, aes(y=structure, x=catch_aav, color=structure, group=sp))+
    stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
    # scale_fill_brewer(na.translate = FALSE) +
    facet_grid(sp ~ sigma_r)+
    ggtitle("Catch AAV")+
    theme_bw()

# Catch ABI: age structure of catch
catch_abi_stability <- catch_abi_data %>%
    group_by(sim, steepness, structure, sp, sigma_r, om) %>%
    summarise(abi_aav = aav(abi))

ggplot(catch_abi_stability, aes(y=structure, x=abi_aav, color=structure, group=sp))+
    stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
    # scale_fill_brewer(na.translate = FALSE) +
    facet_grid(sp ~ sigma_r, scales="free_x")+
    ggtitle("Catch ABI AAV")+
    theme_bw()

# Catch ABI Distribution: distribution of ABI levels 
# when stock is near equilibrium biomass levels
catch_abi_dist <- catch_abi_data %>%
    group_by(steepness, sigma_r, sp) %>%
    mutate(
        ASref = median(abi[which(time > 50 & time < 100 & structure == 1)]),
        depletion = abi/ASref
    ) %>%
    filter(time > 50, time < 100, om == "Equilibrium") %>%
    group_by(sigma_r, steepness, sp, structure) %>%
    median_qi(abi, .width=c(0.50, 0.80, 0.95, 1))

ggplot(catch_abi_dist, aes(y=structure, x=abi, xmin=.lower, xmax=.upper, color=structure, group=sp))+
    geom_pointinterval(position="dodge")+
    # scale_fill_brewer(na.translate = FALSE) +
    facet_grid(sp ~ sigma_r, scales="free_x")+
    ggtitle("Catch ABI Distribution")+
    theme_bw()

# Catch Value Distribution: distribution of catch value 
# when stock is near equilibrium biomass levels
catch_value_dist <- catch_value_data %>%
    group_by(steepness, sigma_r, sp) %>%
    mutate(
        Valref = median(rel_value[which(time > 50 & time < 100 & structure == 1)])
    ) %>%
    filter(time > 50, time < 100, om == "Equilibrium") %>%
    group_by(sigma_r, steepness, sp, structure) %>%
    median_qi(rel_value, .width=c(0.50, 0.80, 0.95, 1))

ggplot(catch_value_dist, aes(y=structure, x=rel_value, xmin=.lower, xmax=.upper, color=structure, group=sp))+
    geom_pointinterval(position="dodge")+
    # scale_fill_brewer(na.translate = FALSE) +
    facet_grid(sp ~ sigma_r, scales="free_x")+
    ggtitle("Catch Value Distribution")+
    theme_bw()


##########################################################################
##########################################################################
##############                Age Structure                 ##############
##########################################################################
##########################################################################
abi_data <- read_csv("data/abi_data_full2.csv", col_types="iffffdff") %>% filter(steepness == 0.75)

# ABI Distribution: distribution of ABI levels when
# stock is near equilibrium biomass levels
abi_dist <- abi_data %>%
    group_by(steepness, sigma_r, sp) %>%
    mutate(
        ASref = median(abi[which(time > 50 & time < 174 & structure == 1)]),
        depletion = abi/ASref
    ) %>%
    filter(time > 50, time < 174, om == "Equilibrium") %>%
    group_by(sigma_r, steepness, sp, structure) %>%
    mean_qi(abi, .width=c(0.50, 0.80, 0.95, 1))

ggplot(abi_dist, aes(y=structure, x=abi, xmin=.lower, xmax=.upper, color=structure, group=sp))+
    geom_pointinterval(position="dodge")+
    # scale_fill_brewer(na.translate = FALSE) +
    facet_grid(sp ~ sigma_r, scales="free_x")+
    ggtitle("ABI Distribution")+
    theme_bw()

# ABI Timeseries: timeseries plots of ABI
abi_tseries <- abi_data %>% as_tibble() %>%
    group_by(steepness, sigma_r, sp) %>%
    mutate(ASref = median(abi[which(time > 50 & time < 100 & structure == 1)])) %>%
    group_by(time, steepness, structure, sp, sigma_r, om) %>%
    median_qi(abi, ASref, .width=c(0.50))

ggplot(abi_tseries)+
    geom_line(aes(x=time, y=abi, color=structure, group=structure), size=0.9)+
    geom_hline(aes(yintercept=ASref), linetype="dashed")+
    facet_grid(sp ~ sigma_r)+
    coord_cartesian(ylim=c(0, 3.5))+
    ggtitle("ABI Timeseries")+
    theme_bw()

# ABI Stability: average annual variation in ABI
abi_stability <- abi_data %>%
    group_by(sim, steepness, structure, sp, sigma_r, om) %>%
    summarise(abi_aav = aav(abi))

ggplot(abi_stability, aes(y=structure, x=abi_aav, color=structure, group=sp))+
    stat_pointinterval(aes(color_ramp = after_stat(level)), position="dodge", .width = c(.50, .80, 0.99, 1))+
    # scale_fill_brewer(na.translate = FALSE) +
    facet_grid(sp ~ sigma_r, scales="free_x")+
    ggtitle("ABI AAV")+
    theme_bw()

# Diversity Distribution: distribution of ABI levels when
# stock is near equilibrium biomass levels
popdiversity_data <- read_csv("data/popdiv_data_full2.csv", col_types="iffffdff") %>% filter(steepness == 0.75)

diversity_dist <- popdiversity_data %>%
    group_by(steepness, sigma_r, sp) %>%
    mutate(
        Dref = median(div[which(time > 50 & time < 100 & structure == 1)])
    ) %>%
    filter(time > 50, time < 100, om == "Equilibrium") %>%
    group_by(sigma_r, steepness, sp, structure) %>%
    median_qi(div, .width=c(0.50, 0.80, 0.95, 1))

ggplot(diversity_dist, aes(y=structure, x=div, xmin=.lower, xmax=.upper, color=structure, group=sp))+
    geom_pointinterval(position="dodge")+
    # scale_fill_brewer(na.translate = FALSE) +
    facet_grid(sp ~ sigma_r, scales="free_x")+
    ggtitle("Shannon Diversity Distribution")+
    theme_bw()

##########################################################################
##########################################################################
##############             Performance/Utility              ##############
##########################################################################
##########################################################################
relativize_performance <- function(df){
    df %>%
        mutate(rebuild_time = ifelse(is.infinite(rebuild_time), -100, rebuild_time)) %>%
        group_by(sigma_r, steepness, sp, sim) %>%
        mutate(
            across(catch:rebuild_rate, ~ (. - .[which(structure == 1)])/ .[which(structure == 1)])
        ) %>%
        mutate_at(vars(ends_with("aav")), ~ -1*.) %>%
        mutate_at(vars(matches("rebuild_time")), ~ -1*.)
}

performance_metrics <- catch_data %>%
    group_by(sim, sigma_r, steepness, sp, structure) %>%
    summarise(catch = mean(catch)) %>%
    left_join(
        ssb_stability %>%
            group_by(sim, sigma_r, steepness, sp, structure) %>%
            summarise(ssb_aav = mean(ssb_aav)),
        by = c("sim", "sp", "structure", "sigma_r", "steepness")
    ) %>%
    left_join(
        abi_stability %>%
            group_by(sim, sigma_r, steepness, sp, structure) %>%
            summarise(abi_aav = mean(abi_aav)),
        by = c("sim", "sp", "structure", "sigma_r", "steepness")
    ) %>%
    left_join(
        catch_stability %>%
            group_by(sim, sigma_r, steepness, sp, structure) %>%
            summarise(catch_aav = mean(catch_aav)),
        by = c("sim", "sp", "structure", "sigma_r", "steepness")
    ) %>%
    left_join(
        catch_value_data %>%
            group_by(sim, sigma_r, steepness, sp, structure) %>%
            summarise(catch_value = mean(rel_value)),
        by = c("sim", "sp", "structure", "sigma_r", "steepness")
    ) %>%
    left_join(
        catch_abi_data %>%
            group_by(sim, sigma_r, steepness, sp, structure) %>%
            summarise(catch_abi = mean(abi, na.rm=TRUE)),
        by = c("sim", "sp", "structure", "sigma_r", "steepness")
    ) %>%
    # left_join(
    #     productivity,
    #     by = c("sim", "sp", "structure", "sigma_r", "steepness")
    # ) %>%
    left_join(
        ssb_data %>%
            group_by(sim, sigma_r, steepness, sp, structure) %>%
            summarise(ssb = mean(ssb)),
        by = c("sim", "sp", "structure", "sigma_r", "steepness")
    ) %>%
    left_join(
        popdiversity_data %>%
            group_by(sim, sigma_r, steepness, sp, structure) %>%
            summarise(div = mean(div)),
        by = c("sim", "sp", "structure", "sigma_r", "steepness")
    ) %>%
    left_join(
        rebound_timerate %>%
            group_by(sim, sigma_r, steepness, sp, structure) %>%
            summarise(
                rebuild_time = mean(rebuild_time),
                rebuild_rate = mean(rebuild_rate)
            ),
        by = c("sim", "sp", "structure", "sigma_r", "steepness")
    )
 
p <- performance_metrics %>%
    relativize_performance() %>%
    group_by(sigma_r, steepness, sp, structure) %>%
    median_qi(ssb, ssb_aav, catch, catch_aav, catch_abi, catch_value, div, rebuild_time, .width=c(0.50, 0.80, 0.95)) %>%  
    reformat_ggdist_long(n=4) %>%
    mutate(type = ifelse(name %in% c("catch", "catch_aav", "catch_abi", "catch_value"), "fishery", "population")) %>%
    mutate(
        name = factor(name, levels = c("catch", "catch_aav", "catch_abi", "catch_value", "ssb", "ssb_aav", "div", "recruitent", "rebuild_time"), labels=c("Catch", "Catch AAV", "Catch ABI", "Catch Value", "SSB", "SSB AAV", "Shannon Diversity", "Average Recruitment", "Rebuilding Time"))
    )

# Population Dynamics Performance: performance plots for 
# population related performance metrics
ggplot(p %>% filter(type == "population"), aes(x=median, xmin=lower, xmax=upper, y=structure, color=sp, group=steepness, shape=steepness))+
    geom_pointrange(position=position_dodge(width = 1))+
    geom_vline(xintercept=0, linetype="dashed")+
    facet_grid(sigma_r ~ name, scales="free_x")+
    theme_bw()

# Fishery Performance: performance plots for fishery related 
# performance metrics
ggplot(p %>% filter(type == "fishery"), aes(x=median, xmin=lower, xmax=upper, y=structure, color=sp, group=steepness, shape=steepness))+
    geom_pointrange(position=position_dodge(width = 1))+
    geom_vline(xintercept=0, linetype="dashed")+
    facet_grid(sigma_r ~ name, scales="free_x")+
    theme_bw()

# Total Performance: performance plots for all performance metrics
ggplot(p, aes(x=median, xmin=lower, xmax=upper, y=structure, color=sp, shape=steepness))+
    geom_pointrange(position=position_dodge(width = 1))+
    # geom_vline(xintercept=0, linetype="dashed")+
    facet_grid(sigma_r ~ name, scales="free_x")+
    guides(shape="none")+
    custom_theme+
    labs(x="Value", y="ABI40%", color="Spawning Potential", title="Performance Metrics")
ggsave('~/Desktop/performance.jpg', width=18, height=12, units="in")

##########################################################################
##########################################################################
##############             K-Means Clustering               ##############
##########################################################################
##########################################################################

fishery_variables <- c("catch3", "catch_aav", "catch_abi", "catch_value")
pop_variables <- c("ssb", "ssb_aav", "div", "rebuild_time")

compute_kmean_data <- function(variables){

    match_indices <- function(x){
        uni = unique(x)
        return(match(x, uni))
    }

    performance_metrics %>%
        relativize_performance() %>%
        mutate(
            catch1 = 2*(catch),
            catch2 = -2*(catch^2),
            catch3 = (catch-0.1)^3,
        ) %>% 
        ungroup() %>% 
        select(c("sim", "sigma_r", "steepness", "sp", variables)) %>%
        nest(variables) %>%
        mutate(
            model = map(data, kmeans, centers=5),    # iterate model over list column of input data
            assignments = map(model, \(x) data.frame(x["cluster"]) %>% mutate(groupid = match_indices(cluster)) %>% pull(groupid)),
        ) %>% unnest(assignments) %>% unnest(assignments) %>% mutate(structure = factor(rep(unique(performance_metrics$structure), length.out=3420)))

}

fishery_kmean = compute_kmean_data(fishery_variables)
pop_kmean = compute_kmean_data(pop_variables)
tot_kmean = compute_kmean_data(c(fishery_variables, pop_variables))

# Age Cluster Performance Metric Boxplots: boxplots of
# distribution of performance metrics within each identified
# age structure cluster
performance_metrics %>% ungroup() %>% 
    arrange(sim, sigma_r, steepness, sp, structure) %>% 
    left_join(
        tot_kmean %>% select(-c(data, model)),
        by = c("sim", "sigma_r", "steepness", "sp", "structure")
    ) %>%
    pivot_longer(cols=c("ssb", "ssb_aav", "div", "rebuild_time", "recruitent", "catch", "catch_aav", "catch_abi"), names_to="name", values_to="value") %>%
    mutate(
        assignments = factor(assignments, labels=c("Very Young", "Young", "MSY", "Old", "Very Old")),
        name = factor(name, levels = c("catch", "catch_aav", "catch_abi", "ssb", "ssb_aav", "div", "recruitent", "rebuild_time"), labels=c("Catch", "Catch AAV", "Catch ABI", "SSB", "SSB AAV", "Shannon Diversity", "Average Recruitment", "Rebuilding Time"))
    ) %>%

    ggplot(aes(x=assignments, y=value, fill=factor(assignments), group=interaction(assignments, sp)))+
        geom_boxplot()+
        facet_grid(name ~ sigma_r, scales='free_y')+
        theme_bw()+
        custom_theme+
        labs(x="Age Group Cluster", y="Value", color="Age Group Cluster", title="In Cluster Performance Metric Distributions")

ggsave("~/Desktop/boxplots.jpg", width=14, height=18, units="in")

# Age Structure Cluster Assignments: K-Means cluster assignments
# for each age structure policy
performance_metrics %>% ungroup() %>% 
    arrange(sim, sigma_r, steepness, sp, structure) %>% 
    left_join(
        tot_kmean %>% select(-c(data, model)),
        by = c("sim", "sigma_r", "steepness", "sp", "structure")
    ) %>%
    group_by(sigma_r, sp, structure) %>%
    mode_qi(assignments) %>%
    mutate(assignments = factor(assignments, labels=c("Very Young", "Young", "MSY", "Old", "Very Old"))) %>%

    ggplot(aes(x=assignments, xmin=.lower, xmax=.upper, y=structure, color=structure)) +
        geom_point(size=3)+
        facet_grid(sigma_r ~ sp)+
        custom_theme+
        guides(color="none")+
        labs(x="Age Group Cluster", y="ABI40%", color="ABI40%", title="Age Structure Clusters")

ggsave("~/Desktop/age_group_clusters.jpg", width=14, height=12.5, units="in")

# Utility Plots: population, fishery, and total utility for each
# age structure policy colored by age structure cluster 
p1 <- performance_metrics %>%
    relativize_performance() %>%
    mutate(
        catch1 = 0, #catch,
        catch2 = -2*(catch^2),
        catch3 = 0, #3*(catch-0.1)^3,
        fish_util = catch_aav+catch1+catch2+catch3+catch_abi+catch_value,
        pop_util = ssb+ssb_aav+div+rebuild_time,
        util = (fish_util+pop_util)/2,
    ) %>%
    group_by(sigma_r, steepness, sp, structure) %>%
    median_qi(fish_util, pop_util, util, .width=c(0.50, 0.80, 0.95)) %>%
    reformat_ggdist_long(n=4) %>%
    arrange(sigma_r, steepness, sp, structure) %>% 
    # left_join(
    #     tot_kmean %>% select(-c(data, model)) %>% group_by(sigma_r, steepness, sp, structure) %>% mode_qi(assignments) %>% select(sigma_r, steepness, sp, structure, assignments),
    #     by = c("sigma_r", "steepness", "sp", "structure")
    # ) %>%
    mutate(
        # assignments = factor(assignments, labels=c("Very Young", "Young", "MSY", "Old", "Very Old")),
        name = factor(name, labels = c("Fishery Utility", "Population Utility", "Total Utility"))
    )


pmax <- p1 %>% filter(.width == 0.5, name == "Total Utility") %>% group_by(sigma_r, sp) %>% summarise(max_util = max(median))

ggplot(p1, aes(x=median, xmin=lower, xmax=upper, y=structure, shape=name, color=structure, group=structure))+
    geom_pointrange(position=position_dodge(width = 1))+
    geom_vline(xintercept=0, linetype="dashed")+
    geom_vline(data=pmax, aes(xintercept=max_util), linetype="dotted")+
    facet_grid(sp ~ sigma_r)+
    coord_cartesian(xlim=c(-1.5, 1.5))+
    scale_shape_manual(values=c(16, 1, 4))+
    custom_theme+
    labs(x="Utility", y="ABI40%", color="Age Structure", title="ABI Harvest Strategy Utility", shape="Utility Function")

ggsave("~/Desktop/utility.jpg", width=14, height=12.5, units="in")

# Aggregated Utility Plots: population, fishery, and total utility
# aggregated into age structure clusters
performance_metrics %>%
    relativize_performance() %>%
    mutate(
        catch1 = 2*catch,
        catch2 = 0, #-2*(catch^2),
        catch3 = 0, #3*(catch-0.1)^3,
        fish_util = catch_aav+catch1+catch2+catch3+catch_abi,
        pop_util = ssb+ssb_aav+div+rebuild_time,
        util = (fish_util+pop_util)/2,
    ) %>%
    arrange(sigma_r, steepness, sp, structure) %>% 
    left_join(
        tot_kmean %>% select(-c(data, model)) %>% group_by(sigma_r, steepness, sp, structure) %>% mode_qi(assignments) %>% select(sigma_r, steepness, sp, structure, assignments),
        by = c("sigma_r", "steepness", "sp", "structure")
    ) %>%
    mutate(assignments = factor(assignments, labels=c("Very Young", "Young", "MSY", "Old", "Very Old"))) %>%
    group_by(sigma_r, steepness, sp, assignments) %>%
    select(sim, sigma_r, steepness, sp, structure, assignments, fish_util, pop_util, util) %>%
    median_qi(fish_util, pop_util, util, .width=c(0.50, 0.80, 0.95)) %>%
    reformat_ggdist_long(n=4) %>%

    ggplot(aes(x=median, xmin=lower, xmax=upper, y=assignments, color=assignments, group=assignments, shape=name))+
        geom_pointrange(position=position_dodge(width = 1))+
        geom_vline(xintercept=0, linetype="dashed")+
        facet_grid(sp ~ sigma_r)+
        coord_cartesian(xlim=c(-1, 2))+
        scale_shape_manual(values=c(16, 1, 4))+
        theme_bw()
