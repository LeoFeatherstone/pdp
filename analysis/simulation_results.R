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
#tree_stats <- readRDS("analysis/processed_simulation_posteriors/tree_imbalance.RData")

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
# tree_stats <- tree_stats %>%
#     mutate(
#         id = paste(
#             organism, treePrior, clockodel,
#             resolution, replicate,
#             sep = "_"
#         )
#     ) %>%
#     filter(id %in% ess_passed) %>%
#     ## REVERSING MISTAKES, I will be able to delete the next chunk when I rerun
#     rename(
#         "treePrior" = treePrior,
#         "clockModel" = clockodel
#     ) %>%
#     mutate(
#         resolution = case_when(
#             resolution == "Day-Fixed" ~ "Year",
#             resolution == "Month-Fixed" ~ "Month",
#             resolution == "Year-Fixed" ~ "Day"
#         )
#     )

# ## Making species a factor to plot with
# tree_stats$factor <- factor(
#     levels = c("h1n1", "sars-cov-2", "saureus", "tb"),
#     tree_stats$organism,
#     labels = c(
#         "H1N1", "SARS-CoV-2",
#         "italic(S.~aureus)", "italic(M.~tuberculosis)"
#     )
# )
# ## Filter cols for binding later
# tree_stats <- tree_stats %>%
#     select(!c(organism, id))

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

## Reproductive Number Parallel Coords
reproductive_number_plot <- (posterior %>%
    filter(
        !(
            (organism == "h1n1" | organism == "sars-cov-2")
            & resolution == "Year"
        )
    ) %>%
    group_by(factor, treePrior, resolution, replicate) %>%
    summarise(
        mean_R0 = mean(reproductiveNumber), # a for legend ordering
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
    arrange(interval) %>%
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
            ~italic(R[0])~Birth~Death, ~italic(R[0])~Coalescent~Exponential,
            ~italic(R[e[1]])~Birth~Death, ~italic(R[e[2]])~Birth~Death
         ),
         values = c(
            "mean_R0.BD" = "dodgerblue", "mean_Re1.BD" = "purple",
            "mean_Re2.BD" = "green", "mean_R0.CE" = "darkorange"
        ),
        limits = c(
            "mean_R0.BD", "mean_R0.CE",
            "mean_Re1.BD", "mean_Re2.BD"
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
    )
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

### Variance plots - repeat above with variance instead of mean
variance_data <- posterior %>%
    group_by(factor, treePrior, resolution, replicate) %>%
    summarise(
        var_R0 = var(reproductiveNumber),
        var_Re1 = var(reproductiveNumber.1),
        var_Re2 = var(reproductiveNumber.2),
        var_outbreak_age = var(Tree.height),
        var_clock_rate = var(clockRate)
    ) %>%
    pivot_longer(
            matches("var_R0|var_Re1|var_Re2"),
            names_to = "interval",
            values_to = "var_reproductive_number",
            values_drop_na = TRUE
    ) %>%
    filter(
        !(factor == "H1N1" & resolution == "Year")
        &
        !(factor == "SARS-CoV-2" & resolution == "Year")
    ) 

## Origin Error vs Variance
origin_variance <- variance_data %>%
    ggplot(aes(x = resolution, y = var_outbreak_age)) +
        geom_line(
            aes(
                group = interaction(factor, replicate, treePrior),
                col = treePrior
            ),
            alpha = 0.5
        ) +
        geom_boxplot(
        aes(
            fill = treePrior
        ),
        position = "dodge", alpha = 0.3, outlier.shape = NA
    ) +
        facet_wrap(~factor, scales = "free", labeller = label_parsed, ncol = 4) +
        xlab("Date Resolution") +
        ylab("Variance in tMRCA") +
            scale_discrete_manual(
        aesthetics = c("col", "fill"),
        labels = c("Birth Death", "Coalescent Exponential"),
        values = c("dodgerblue", "darkorange")
    ) +
        theme_bw() +
        theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            text = element_text(size = 14),
            panel.grid.minor = element_blank()
        )

## Clock Rate Error vs Variance
clock_variance <- variance_data %>%
    ggplot(aes(x = resolution, y = var_clock_rate, col = treePrior)) +
        geom_line(
            aes(
                group = interaction(factor, replicate, treePrior),
                col = treePrior
            ),
            alpha = 0.5
        ) +
            geom_boxplot(
        aes(
            fill = treePrior
        ),
        position = "dodge", alpha = 0.3, outlier.shape = NA
    ) +
        facet_wrap(~factor, scales = "free", labeller = label_parsed, ncol = 4) +
        xlab("Date Resolution") +
        ylab("Variance in substitution rate") +
            scale_discrete_manual(
        aesthetics = c("col", "fill"),
        labels = c("Birth Death", "Coalescent Exponential"),
        values = c("dodgerblue", "darkorange")
    ) +
        theme_bw() +
        theme(
            legend.title = element_blank(),
            legend.position = "none",
            text = element_text(size = 14),
            panel.grid.minor = element_blank()
        )

## Reproductive Number Error vs Variance
reproductive_number_variance <- variance_data %>%
    ggplot(aes(x = resolution, y = var_reproductive_number, col = interaction(interval, treePrior))) +
    geom_line(
            aes(
                group = interaction(factor, interval, replicate, treePrior),
                col = interaction(interval, treePrior)
            ),
            alpha = 0.5
        ) +
        geom_boxplot(
        aes(
            x = resolution, y = var_reproductive_number,
            group = interaction(interval, resolution, treePrior),
            fill = interaction(interval, treePrior)
        ),
        position = "dodge", alpha = 0.3,
        outlier.shape = NA
    ) +
        facet_wrap(~factor, scales = "free", labeller = label_parsed, ncol = 4) +
        xlab("Date Resolution") +
        ylab(latex2exp::TeX("Variance in $\\textit{R_{\\bullet}}$")) +
        scale_discrete_manual(
        aesthetics = c("col", "fill"),
        labels = c(
            ~italic(R[0])~Birth~Death, ~italic(R[0])~Coalescent~Exponential,
            ~italic(R[e[1]])~Birth~Death, ~italic(R[e[2]])~Birth~Death
         ),
         values = c(
            "var_R0.BD" = "dodgerblue", "var_Re1.BD" = "purple",
            "var_Re2.BD" = "green", "var_R0.CE" = "darkorange"
        ),
        limits = c(
            "var_R0.BD", "var_R0.CE",
            "var_Re1.BD", "var_Re2.BD"
        )
    ) +
        theme_bw() +
        theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            text = element_text(size = 14),
            panel.grid.minor = element_blank()
        )


## Simulation Combined Parms plot
var_plot <- cowplot::plot_grid(
    labels = "AUTO",
    clock_variance, origin_variance, reproductive_number_variance,
    align = "v",
    nrow = 3
)
ggsave(
    plot = var_plot,
    file = "figures/simulation_variance_panel.pdf",
    dpi = 300
)

### Likelihood plots. TODO Make interactive later
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
## Commented out b/c result not informative 09-12-24
# ei_plot <- (rel_error %>%
#     filter(
#         !(
#             (factor == "H1N1" | factor == "SARS-CoV-2")
#             & resolution == "Year"
#         )
#     ) %>%
#     group_by(factor, treePrior, resolution, replicate) %>%
#     ggplot() +
#     geom_line(
#         aes(
#             x = resolution,
#             y = ei_ratio,
#             group = interaction(factor, replicate, treePrior),
#             col = interaction(treePrior)
#         ),
#         alpha = 1, linewidth = 0.5
#     ) +
#     facet_wrap(
#         ~factor, nrow = 1,
#         scales = "free",
#         labeller = label_parsed
#     ) +
#     ylab("External to Internal Branch Length Ratio") +
#     xlab("Date resolution") +
#     theme_bw() +
#     theme(
#         legend.title = element_blank(),
#         legend.position = "bottom",
#         panel.grid.minor = element_blank()
#     ))
# ei_plot
ggsave("figures/ei_ratio.pdf", dpi = 300)
