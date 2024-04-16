###############################################################
########## Figures and results for simulation study ###########
###############################################################

library(tidyverse)
library(latex2exp)
library(plotly)

## TODO: Update path later
trace <- readRDS("analysis/processed_simulation_posteriors/posteriors.RData")
ess <- readRDS("analysis/processed_simulation_posteriors/ess.RData")
#treeStats <- readRDS("analysis/processed_simulation_posteriors/ess.RData")

## Find which burnin percentage yielded best convergence
for (list in ess) {
    print(
        sum(
            sapply(list, function(x) {
                all(x > 200, na.rm = TRUE)
            })
        )
    )
} # ess 50% gave the best results (ess[[9]])

ess_passed <- names(which(unlist(
    lapply(
        ess[[9]],
        function(x) {
            all(x > 200, na.rm = TRUE)
        }
    )
)))

write.table(
    table(gsub(ess_passed, pattern = "-Fixed.+", replacement = "")),
    row.names = FALSE, col.names = FALSE, quote = FALSE,
    file = "/data/cephfs/punim0819/leo/pdp/analysis/simulation_results/ess_50pc.txt"
)

## Filter out non-converged chains
posterior <- trace %>%
    mutate(
        id = paste(
            organism, clockModel, treePrior,
            resolution, replicate,
            sep = "_"
        )
    ) %>%
    filter(id %in% ess_passed) %>%
    filter(Sample >= (0.5 * 5e8))

## REVERSING MISTAKES, I will be able to delete the next chunk when I rerun
posterior <- posterior %>%
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
    )

## Making species a factor to plot with
posterior$factor <- factor(
    levels = c("h1n1", "sars-cov-2", "saureus", "tb"),
    posterior$organism,
    labels = c(
        "H1N1", "SARS-CoV-2",
        "italic(S.~aureus)", "italic(M.~tuberculosis)"
    )
)

## Clock Parallel Coords
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
    rate = c(4e-3, 1e-3, 1e-6, 1e-8)
)
clock_plot <- (posterior %>%
    group_by(factor, treePrior, resolution, replicate) %>%
    summarise(mean_clock_rate = mean(clockRate)) %>%
    ggplot() +
    geom_segment(
        data = rate_true,
        aes(y = rate, yend = rate, x = -Inf, xend = Inf),
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
        panel.grid.minor = element_blank()
    ))
clock_plot
ggsave("clock_traj.pdf", dpi = 300, width = 10, height = 5, units = "in")

## Origin Parallel Coords
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
    origin = c(0.25, 0.16, 25, 25)
)
origin_plot <- (posterior %>%
    group_by(factor, treePrior, resolution, replicate) %>%
    summarise(mean_outbreak_age = mean(Tree.height)) %>%
    ggplot() +
    geom_segment(
        data = origin_true,
        aes(y = origin, yend = origin, x = -Inf, xend = Inf),
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
    ylab(TeX("Outbreak age (months or years)")) +
    xlab("Date resolution") +
    #guides(fill = "none", col = "none") +
    theme_bw() +
    theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    ))
origin_plot
ggsave("origin_traj.pdf", dpi = 300, width = 10, height = 5, units = "in")




## Reproductive Number Parallel Coords
reproductive_number_true <- data.frame(
    organism = c("h1n1", "sars-cov-2", "saureus", "saureus", "tb", "tb"),
    reproductive_number = c(1.3, 2.5, 2.0, 1.0, 2.5, 1.1)
)
reproductive_number_true$factor <- factor(
    reproductive_number_true$organism,
    levels = c("h1n1", "sars-cov-2", "saureus", "tb"),
        labels = c(
            "H1N1", "SARS-CoV-2",
            "italic(S.~aureus)", "italic(M.~tuberculosis)"
        )
    )
# R0 for CE and combine reproductive number values
posterior <- (posterior %>%
    mutate(reproductiveNumber = case_when(
        treePrior == "CE" & organism == "h1n1" ~ (growthRate / 91.3125) + 1,
        treePrior == "CE" & organism == "sars-cov-2" ~ (growthRate / 36.525) + 1,
        .default = reproductiveNumber
    )))

reproductive_number_plot <- (posterior %>%
    filter(
        !((organism == "h1n1" | organism == "sars-cov-2") & resolution == "Year")
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
            y = reproductive_number, yend = reproductive_number,
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
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    ))
reproductive_number_plot
ggsave(
    "reproductive_number_traj.pdf", dpi = 300,
    width = 10, height = 5, units = "in"
)

## Simulation Combines Parms plot
cowplot::plot_grid(
    labels = "AUTO",
    clock_plot, origin_plot, reproductive_number_plot,
    align = "v",
    nrow = 3
)
ggsave(
    "simulation_parm_panel.pdf",
    dpi = 300
)

## Likelihood plots - make interactive
likelihood_plot <- posterior %>%
    #filter(replicate == 1) %>%
    filter(factor == "SARS-CoV-2" & treePrior == "BD") %>%
    mutate(
        phylogeneticLikelihood = posterior - likelihood - prior,
        phylodynamicLikelihood = likelihood
    ) %>%
    ggplot(
        aes(
            x = phylogeneticLikelihood, y = phylodynamicLikelihood,
            group = interaction(resolution, treePrior, factor),
            col = interaction(resolution, treePrior, factor),
            frame = replicate
        )
    ) +
    geom_density2d() +
    facet_wrap(
        ~factor, nrow = 2,
        scales = "free",
        labeller = label_parsed
    ) +
    #ylab(TeX("\\textit{$R_{\\bullet}$}")) +
    #xlab("Date resolution") +
    coord_fixed() +
    theme_bw() +
    theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    )

likelihood_widget <- ggplotly(likelihood_plot)
#ggsave("simulation_likelihood.pdf", dpi = 300, plot = likelihood_plot)
htmlwidgets::saveWidget(likelihood_widget, "simulation_likelihood.html")


## Tree Statistic Plots
load("analysis/processed_simulation_posteriors/tree_imbalance.RData")
