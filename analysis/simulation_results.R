###############################################################
########## Figures and results for simulation study ###########
###############################################################

library(tidyverse)
library(latex2exp)
library(plotly)
library(htmlwidgets)

## TODO: Update path later
trace <- readRDS("analysis/processed_simulation_posteriors/posteriors.RData")
ess <- readRDS("analysis/processed_simulation_posteriors/ess.RData")
tree_stats <- readRDS("analysis/processed_simulation_posteriors/tree_imbalance.RData")

# ess 50% gave the best results (ess[[9]])
ess_passed <- names(which(unlist(
    lapply(
        ess[[9]],
        function(x) {
            all(x > 200, na.rm = TRUE)
        }
    )
)))

### Wrangle posteriors
## Filter out non-converged chains
posterior <- trace %>%
    mutate(
        id = paste(
            organism, clockModel, treePrior,
            resolution, replicate,
            sep = "_"
        )
    ) %>%
    filter(Sample >= (0.5 * 5e8)) %>%
    filter(id %in% ess_passed) %>%
    ## REVERSING MISTAKES, I will be able to delete the next chunk when I rerun
    rename(
        "treePrior" = clockModel,
        "clockModel" = treePrior
    ) %>%
    mutate(
        resolution = case_when(
            resolution == "Day-Fixed" ~ "Year",
            resolution == "Month-Fixed" ~ "Month",
            resolution == "Year-Fixed" ~ "Day"
        )
    ) %>%
    ## R0 for CE and combine reproductive number values
    mutate(reproductiveNumber = case_when(
        treePrior == "CE" & organism == "h1n1" ~ (growthRate / 91.3125) + 1,
        treePrior == "CE" & organism == "sars-cov-2" ~ (growthRate / 36.525) + 1,
        .default = reproductiveNumber
    ))

## Making species a factor to plot with
posterior$factor <- factor(
    levels = c("h1n1", "sars-cov-2", "saureus", "tb"),
    posterior$organism,
    labels = c(
        "H1N1", "SARS-CoV-2",
        "italic(S.~aureus)", "italic(M.~tuberculosis)"
    )
)

### Handle Tree Stats
tree_stats <- tree_stats %>%
    mutate(
        id = paste(
            organism, treePrior, clockodel,
            resolution, replicate,
            sep = "_"
        )
    ) %>%
    filter(id %in% ess_passed) %>%
    ## REVERSING MISTAKES, I will be able to delete the next chunk when I rerun
    rename(
        "treePrior" = treePrior,
        "clockModel" = clockodel
    ) %>%
    mutate(
        resolution = case_when(
            resolution == "Day-Fixed" ~ "Year",
            resolution == "Month-Fixed" ~ "Month",
            resolution == "Year-Fixed" ~ "Day"
        )
    )

## Making species a factor to plot with
tree_stats$factor <- factor(
    levels = c("h1n1", "sars-cov-2", "saureus", "tb"),
    tree_stats$organism,
    labels = c(
        "H1N1", "SARS-CoV-2",
        "italic(S.~aureus)", "italic(M.~tuberculosis)"
    )
)
## Filter cols for binding later
tree_stats <- tree_stats %>%
    select(!c(organism, id))

## Simulation values to compare against
rate_true <- data.frame(
    organism = c("h1n1", "sars-cov-2", "saureus", "tb"),
    factor = factor(
        levels = c("h1n1", "sars-cov-2", "saureus", "tb"),
        c("h1n1", "sars-cov-2", "saureus", "tb"),
        labels = c(
            "H1N1", "SARS-CoV-2",
            "italic(S.~aureus)", "italic(M.~tuberculosis)"
        )
    ),
    rate_true = c(4e-3, 1e-3, 1e-6, 1e-8)
)
origin_true <- data.frame(
    organism = c("h1n1", "sars-cov-2", "saureus", "tb"),
    factor = factor(
        levels = c("h1n1", "sars-cov-2", "saureus", "tb"),
        c("h1n1", "sars-cov-2", "saureus", "tb"),
        labels = c(
            "H1N1", "SARS-CoV-2",
            "italic(S.~aureus)", "italic(M.~tuberculosis)"
        )
    ),
    origin_true = c(0.25, 0.16, 25, 25)
)
reproductive_number_true <- data.frame(
    organism = c("h1n1", "sars-cov-2", "saureus", "saureus", "tb", "tb"),
    factor = factor(
        levels = c("h1n1", "sars-cov-2", "saureus", "tb"),
        c("h1n1", "sars-cov-2", "saureus", "saureus", "tb", "tb"),
        labels = c(
            "H1N1", "SARS-CoV-2",
            "italic(S.~aureus)", "italic(M.~tuberculosis)"
        )
    ),
    reproductive_number_true = c(1.3, 2.5, 2.0, 1.0, 2.5, 1.1),
    interval = c(
        "mean_R0", "mean_R0",
        "mean_Re1", "mean_Re2", "mean_Re1", "mean_Re2"
    )
)

rel_error <- posterior %>%
    ungroup() %>%
    select(!organism) %>%
    group_by(factor, treePrior, resolution, replicate) %>%
    summarise(
        mean_R0 = mean(reproductiveNumber),
        mean_Re1 = mean(reproductiveNumber.1),
        mean_Re2 = mean(reproductiveNumber.2),
        mean_age = mean(Tree.height),
        mean_clock_rate = mean(clockRate),
        .groups = "keep"
    ) %>%
    pivot_wider(
        names_from = resolution,
        names_glue = "{.value}_{resolution}",
        values_from = starts_with("mean_")
    ) %>%
    mutate(
        #R0
        error_R0_Day = 0,
        error_R0_Month = mean_R0_Month - mean_R0_Day,
        error_R0_Year = mean_R0_Year - mean_R0_Day,
        #Re1
        error_Re1_Day = 0,
        error_Re1_Month = mean_Re1_Month - mean_Re1_Day,
        error_Re1_Year = mean_Re1_Year - mean_Re1_Day,
        #Re2
        error_Re2_Day = 0,
        error_Re2_Month = mean_Re2_Month - mean_Re2_Day,
        error_Re2_Year = mean_Re2_Year - mean_Re2_Day,
        #age
        error_age_Day = 0,
        error_age_Month = mean_age_Month - mean_age_Day,
        error_age_Year = mean_age_Year - mean_age_Day,
        #rate
        error_clock_rate_Day = 0,
        error_clock_rate_Month = mean_clock_rate_Month - mean_clock_rate_Day,
        error_clock_rate_Year = mean_clock_rate_Year - mean_clock_rate_Day,
    ) %>%
    select(!starts_with("mean_")) %>%
    pivot_longer(
        starts_with("error_"),
        names_to = c(".value", "resolution"),
        names_pattern = "(error_.+)_(Day|Month|Year)"
    ) %>%
    pivot_longer(
            matches("error_R(0|e1|e2)"),
            names_to = "interval",
            values_to = "error_reproductive_number",
            values_drop_na = TRUE
    ) %>%
    mutate(interval = sub(interval, pattern = "error_", replacement = "")) %>%
    left_join(
        tree_stats,
        by = c("factor", "replicate", "resolution", "treePrior")
    )

####################################
############# Plotting #############
####################################

## Clock Parallel Coords
clock_plot <- (posterior %>%
    group_by(factor, treePrior, resolution, replicate) %>%
    summarise(mean_clock_rate = mean(clockRate)) %>%
    ggplot() +
    geom_segment(
        data = rate_true,
        aes(y = rate_true, yend = rate_true, x = -Inf, xend = Inf),
        col = "black", alpha = 0.6, lty = 3
    ) +
    geom_line(
        aes(
            x = resolution,
            y = mean_clock_rate,
            group = interaction(factor, treePrior, replicate),
            col = treePrior
        ),
        alpha = 0.5, linewidth = 0.1
    ) +
    geom_boxplot(
        aes(
            x = resolution, y = mean_clock_rate,
            fill = treePrior
        ),
        position = "dodge", alpha = 0.3, outlier.shape = NA
    ) +
    scale_y_continuous(
        trans = "log10",
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_discrete_manual(
        aesthetics = c("col", "fill"),
        labels = c("Birth Death", "Coalescent Exponential"),
        values = c("dodgerblue", "darkorange")
    ) +
    facet_wrap(
        ~factor, nrow = 1,
        scales = "free_y",
        labeller = label_parsed
    ) +
    ylab(TeX("Substitution rate \\textit{(subs/site/time)}")) +
    xlab("Date resolution") +
    guides(fill = "none", col = "none") +
    theme_bw() +
    theme(
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 14)
    ))
clock_plot
ggsave("clock_traj.pdf", dpi = 300, width = 10, height = 5, units = "in")

## Origin Parallel Coords

origin_plot <- (posterior %>%
    group_by(factor, treePrior, resolution, replicate) %>%
    summarise(mean_outbreak_age = mean(Tree.height)) %>%
    ggplot() +
    geom_segment(
        data = origin_true,
        aes(y = origin_true, yend = origin_true, x = -Inf, xend = Inf),
        col = "black", alpha = 0.6, lty = 3
    ) +
    geom_line(
        aes(
            x = resolution,
            y = mean_outbreak_age,
            group = interaction(factor, treePrior, replicate),
            col = treePrior
        ),
        alpha = 0.5, linewidth = 0.1
    ) +
    geom_boxplot(
        aes(
            x = resolution, y = mean_outbreak_age,
            fill = treePrior
        ),
        position = "dodge", alpha = 0.3, outlier.shape = NA
    ) +
    scale_y_continuous(
        breaks = c(0, (1 / 12), seq(0.25, 0.75, by = 0.25), seq(1, 35, by = 5)),
        labels = c(
            0, paste0(c(1, 3, 6, 9), "m"),
            paste0(seq(1, 35, by = 5), "y")
        )
    ) +
    scale_discrete_manual(
        aesthetics = c("col", "fill"),
        labels = c("Birth Death", "Coalescent Exponential"),
        values = c("dodgerblue", "darkorange")
    ) +
    facet_wrap(
        ~factor, nrow = 1,
        scales = "free_y",
        labeller = label_parsed
    ) +
    ylab(TeX("tMRCA (months or years)")) +
    xlab("Date resolution") +
    #guides(fill = "none", col = "none") +
    theme_bw() +
    theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        text = element_text(size = 14)
    ))
origin_plot
ggsave("origin_traj.pdf", dpi = 300, width = 10, height = 5, units = "in")




## Reproductive Number Parallel Coords
reproductive_number_true$factor <- factor(
    reproductive_number_true$organism,
    levels = c("h1n1", "sars-cov-2", "saureus", "tb"),
        labels = c(
            "H1N1", "SARS-CoV-2",
            "italic(S.~aureus)", "italic(M.~tuberculosis)"
        )
    )

reproductive_number_plot <- (posterior %>%
    filter(
        !(
            (organism == "h1n1" | organism == "sars-cov-2")
            & resolution == "Year"
        )
    ) %>%
    group_by(factor, treePrior, resolution, replicate) %>%
    summarise(
        mean_R0 = mean(reproductiveNumber),
        mean_Re1 = mean(reproductiveNumber.1),
        mean_Re2 = mean(reproductiveNumber.2),
    ) %>%
    pivot_longer(
            matches("mean_R0|mean_Re1|mean_Re2"),
            names_to = "interval",
            values_to = "mean_reproductive_number",
            values_drop_na = TRUE
    ) %>%
    group_by(factor, treePrior, resolution, replicate, interval) %>%
    ggplot() +
    geom_segment(
        data = reproductive_number_true,
        aes(
            y = reproductive_number_true, yend = reproductive_number_true,
            x = -Inf, xend = Inf
        ),
        col = "black", alpha = 0.6, lty = 3
    ) +
    geom_line(
        aes(
            x = resolution,
            y = mean_reproductive_number,
            group = interaction(factor, interval, replicate, treePrior),
            col = interaction(interval, treePrior)
        ),
        alpha = 1, linewidth = 0.1
    ) +
    geom_boxplot(
        aes(
            x = resolution, y = mean_reproductive_number,
            group = interaction(interval, resolution, treePrior),
            fill = interaction(interval, treePrior)
        ),
        position = "dodge", alpha = 0.3,
        outlier.shape = NA
    ) +
    scale_y_continuous(
        trans = "log10",
        breaks = c(seq(1, 2.5, by = 0.5), 3, seq(5, 15, by = 5))
    ) +
    scale_discrete_manual(
        aesthetics = c("colour", "fill"),
        labels = c(
            ~italic(R[0])~Birth~Death, ~italic(R[e[1]])~Birth~Death,
            ~italic(R[e[2]])~Birth~Death, ~italic(R[0])~Coalescent~Exponential
         ),
         values = c(
            "mean_R0.BD" = "dodgerblue", "mean_R0.CE" = "darkorange",
            "mean_Re1.BD" = "purple", "mean_Re2.BD" = "green"
        )
    ) +
    facet_wrap(
        ~factor, nrow = 1,
        scales = "free",
        labeller = label_parsed
    ) +
    ylab(TeX("\\textit{$R_{\\bullet}$}")) +
    xlab("Date resolution") +
    theme_bw() +
    theme(
        legend.title = element_blank(),
        text = element_text(size = 14),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    ))
reproductive_number_plot
ggsave(
    "reproductive_number_traj.pdf", dpi = 300,
    width = 10, height = 5, units = "in"
)

## Simulation Combined Parms plot
cowplot::plot_grid(
    labels = "AUTO",
    clock_plot, origin_plot, reproductive_number_plot,
    align = "v",
    nrow = 3
)
ggsave(
    "figures/simulation_parm_panel.pdf",
    dpi = 300
)

## Likelihood plots. TODO Make interactive later
sim_likielihood_plot <- posterior %>%
    pivot_longer(
        matches("BDSKY_Serial|CoalescentExponential"),
        names_to = "beast_phylodynamic_likelihood_name",
        values_to = "phylodynamicLikelihood",
        values_drop_na = TRUE
    ) %>%
    mutate(phylogeneticLikelihood = treeLikelihood) %>%
    ungroup() %>%
    select(!organism) %>%
    group_by(factor, treePrior, resolution, replicate) %>%
    summarise(
        mean_phylogenetic = mean(phylogeneticLikelihood),
        mean_phylodynamic = mean(phylodynamicLikelihood),
        .groups = "keep"
    ) %>%
    pivot_wider(
        names_from = resolution,
        names_glue = "{.value}_{resolution}",
        values_from = starts_with("mean_phylo")
    ) %>%
    mutate(
        # Adjusted Phylogenetic Likelihood
        phylogenetic_Day = 0,
        phylogenetic_Month = mean_phylogenetic_Month - mean_phylogenetic_Day,
        phylogenetic_Year = mean_phylogenetic_Year - mean_phylogenetic_Day,
        # Adjusted Phylodynamic Likelihood
        phylodynamic_Day = 0,
        phylodynamic_Month = mean_phylodynamic_Month - mean_phylodynamic_Day,
        phylodynamic_Year = mean_phylodynamic_Year - mean_phylodynamic_Day,
    ) %>%
    select(!starts_with("mean_phylo")) %>%
    pivot_longer(
        matches("(phylogenetic|phylodynamic)_(Day|Month|Year)"),
        names_to = c(".value", "resolution"),
        names_pattern = "(phylogenetic|phylodynamic)_(Day|Month|Year)"
    ) %>%
    filter(!(factor == "H1N1" & phylodynamic > 25000)) %>% # Just for scale
    ggplot(
        aes(
            x = phylogenetic, y = phylodynamic,
            fill = resolution, group = interaction(replicate, treePrior)
        )
    ) +
    xlab("Phylogenetic likelihood (difference from day)") +
    ylab("Phylodynamic Likelihood (difference from day)") +
    geom_point(shape = 21, alpha = 0.75, size = 3) +
    scale_discrete_manual(
        aesthetics = c("fill", "col"),
        values = c("red", "dodgerblue", "black")
    ) +
    facet_wrap(
        ~factor, nrow = 2,
        scales = "free",
        labeller = label_parsed
    ) +
    theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 24)
    )

ggsave(
    plot = sim_likielihood_plot,
    filename = "figures/simulation_likelihood.pdf",
    dpi = 300
)

### EI branch ratio vs resolution
ei_plot <- (rel_error %>%
    filter(
        !(
            (factor == "H1N1" | factor == "SARS-CoV-2")
            & resolution == "Year"
        )
    ) %>%
    group_by(factor, treePrior, resolution, replicate) %>%
    ggplot() +
    geom_line(
        aes(
            x = resolution,
            y = ei_ratio,
            group = interaction(factor, replicate, treePrior),
            col = interaction(treePrior)
        ),
        alpha = 1, linewidth = 0.5
    ) +
    facet_wrap(
        ~factor, nrow = 1,
        scales = "free",
        labeller = label_parsed
    ) +
    ylab("External to Internal Branch Length Ratio") +
    xlab("Date resolution") +
    theme_bw() +
    theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    ))
ei_plot
ggsave("ei_ratio.pdf", dpi = 300)

### Error plots

## Age error
rel_error %>%
    filter(
        resolution != "Day"
        &
        !(factor == "H1N1" & resolution == "Year")
        &
        !(factor == "SARS-CoV-2" & resolution == "Year")
    ) %>%
    pivot_longer(
            matches("n_tips|coless_imbalance|sampling_span"),
            names_to = "stat",
            values_to = "value",
            values_drop_na = TRUE
    ) %>%
    ggplot(
        aes(
            x = value, y = error_age,
            fill = interaction(interval, treePrior, resolution),
            col = interaction(interval, treePrior, resolution)
        )
    ) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point(shape = 21, alpha = 0.2, col = "black") +
    facet_wrap(
        ~ factor * stat, labeller = label_parsed,
        ncol = 3, scales = "free"
    ) +
    theme(
        legend.position = "bottom", legend.title = element_blank()
    )
ggsave("age_error.pdf", dpi = 300)

## clock error
rel_error %>%
    filter(
        resolution != "Day"
        &
        !(factor == "H1N1" & resolution == "Year")
        &
        !(factor == "SARS-CoV-2" & resolution == "Year")
    ) %>%
    pivot_longer(
            matches("n_tips|coless_imbalance|sampling_span"),
            names_to = "stat",
            values_to = "value",
            values_drop_na = TRUE
    ) %>%
    ggplot(
        aes(
            x = value, y = error_clock_rate,
            fill = interaction(interval, treePrior, resolution),
            col = interaction(interval, treePrior, resolution)
        )
    ) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point(shape = 21, alpha = 0.2, col = "black") +
    facet_wrap(
        ~ factor * stat, labeller = label_parsed,
        ncol = 3, scales = "free"
    ) +
    theme(
        legend.position = "bottom", legend.title = element_blank()
    )
ggsave("clock_rate_error.pdf", dpi = 300)

## reproductive number error
rel_error %>%
    filter(
        resolution != "Day"
        &
        !(factor == "H1N1" & resolution == "Year")
        &
        !(factor == "SARS-CoV-2" & resolution == "Year")
    ) %>%
    pivot_longer(
            matches("n_tips|coless_imbalance|sampling_span"),
            names_to = "stat",
            values_to = "value",
            values_drop_na = TRUE
    ) %>%
    ggplot(
        aes(
            x = value, y = error_reproductive_number,
            fill = interaction(interval, treePrior, resolution),
            col = interaction(interval, treePrior, resolution)
        )
    ) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point(shape = 21, alpha = 0.2, col = "black") +
    facet_wrap(
        ~ factor * stat, labeller = label_parsed,
        ncol = 3, scales = "free"
    ) +
    theme(
        legend.position = "bottom", legend.title = element_blank()
    )
ggsave("reproductie_number_error.pdf", dpi = 300)
