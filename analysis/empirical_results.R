###################################################
############# Emprirical data plots ###############
###################################################

library(tracerer)
library(tidyverse)
library(latex2exp)
conflicted::conflicts_prefer(latex2exp::TeX)
library(cowplot)
conflicted::conflicts_prefer(dplyr::filter())
library(treedataverse)
library(NELSI)
library(phangorn)
library(phytools)
library(kableExtra)


# Get the list of .log files in the directory
log_files <- list.files("analysis/empirical_data", pattern = "\\.log$", full.names = TRUE)

# Read and combine the .log files into lists
traces <- lapply(log_files, tracerer::parse_beast_tracelog_file)
names(traces) <- basename(log_files)

# Merge lists into a single data frame
traces <- traces %>%
    bind_rows(.id = "file") %>%
    separate_wider_delim(
        file,
        delim = "_", names = c("organism", "treePrior", "clock", "resolution"),
    ) %>%
    select(!matches("rate.+|freq.+")) %>%
    filter(Sample > 25000000) %>%
    mutate(resolution = gsub(pattern = "-Fixed[.]log", replacement = "", resolution)) %>%
    mutate(reproductiveNumber = case_when(
        treePrior == "CE" & organism == "h1n1" ~ (growthRate / 91.3125) + 1,
        treePrior == "CE" & organism == "sars-cov-2" ~ (growthRate / 36.525) + 1,
        .default = reproductiveNumber
    )) %>%
    filter(!(organism %in% c("h1n1", "sars-cov-2") & resolution == "Year"))

# Convert organism to factor for plotting
traces$organism <- factor(
    levels = c("h1n1", "sars-cov-2", "saureus", "tb"),
    traces$organism,
    labels = c(
        "H1N1", "SARS-CoV-2",
        "italic(S.~aureus)", "italic(M.~tuberculosis)"
    )
)

clock_plot <- ggplot(traces, aes(x = resolution, y = clockRate, fill = treePrior)) +
    geom_violin(alpha = 0.5, scale = "width", draw_quantiles = c(0.025, 0.5, 0.975)) +
    scale_discrete_manual(
        aesthetics = c("colour", "fill"),
        labels = c(Birth ~ Death, Coalescent ~ Exponential),
        values = c("BD" = "dodgerblue", "CE" = "darkorange")
    ) +
    scale_y_continuous(
        trans = "log10",
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    facet_wrap(~organism, scales = "free", labeller = label_parsed, nrow = 1) +
    ylab(TeX("Substitution rate")) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        text = element_text(size = 14)
    )

# Plot posterior tree height. Facet by organism and colour by resolution
age_plot_covid_ce_unfiltered <- traces %>%
    #filter(!(organism == "SARS-CoV-2" & treePrior == "CE")) %>% # Posterior is flat
    ggplot(aes(x = resolution, y = Tree.height, fill = treePrior)) +
    geom_violin(alpha = 0.5, scale = "width", draw_quantiles = c(0.025, 0.5, 0.975)) +
    scale_discrete_manual(
        aesthetics = c("colour", "fill"),
        labels = c(Birth ~ Death, Coalescent ~ Exponential),
        values = c("BD" = "dodgerblue", "CE" = "darkorange")
    ) +
    scale_y_continuous(
        breaks = c(0, (1 / 12), (2 / 12), seq(0.25, 0.75, by = 0.25), seq(1, 30, by = 5), seq(30, 60, by = 10)),
        labels = c(
            0, paste0(c(1, 2, 3, 6, 9), "m"),
            paste0(seq(1, 30, by = 5), "y"),
            paste0(seq(30, 60, by = 10), "y")
        )
    ) +
    facet_wrap(~organism, scales = "free", labeller = label_parsed, nrow = 1) +
    ylab(TeX("tMRCA")) + # was "Outbreak Age"
    theme_bw() +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 14)
    )
# Plot posterior tree height. Facet by organism and colour by resolution
age_plot_covid_ce_filtered <- traces %>%
    filter(!(organism == "SARS-CoV-2" & treePrior == "CE")) %>% # Posterior is flat
    ggplot(aes(x = resolution, y = Tree.height, fill = treePrior)) +
    geom_violin(alpha = 0.5, scale = "width", draw_quantiles = c(0.025, 0.5, 0.975)) +
    scale_discrete_manual(
        aesthetics = c("colour", "fill"),
        labels = c(Birth ~ Death, Coalescent ~ Exponential),
        values = c("BD" = "dodgerblue", "CE" = "darkorange")
    ) +
    scale_y_continuous(
        breaks = c(0, (1 / 12), (2 / 12), seq(0.25, 0.75, by = 0.25), seq(1, 30, by = 5), seq(30, 60, by = 10)),
        labels = c(
            0, paste0(c(1, 2, 3, 6, 9), "m"),
            paste0(seq(1, 30, by = 5), "y"),
            paste0(seq(30, 60, by = 10), "y")
        )
    ) +
    facet_wrap(~organism, scales = "free", labeller = label_parsed, nrow = 1) +
    ylab(TeX("tMRCA")) + # was "Outbreak Age"
    theme_bw() +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 14)
    )
# Plot the posterior reproductive number. Facet by organism and colour by resolution
reproductive_plot <- traces %>%
    filter(!(organism %in% c("H1N1", "SARS-CoV-2") & resolution == "Year")) %>%
    pivot_longer(
        cols = starts_with("reproductiveNumber"),
        names_to = "interval", values_to = "reproductiveNumber"
    ) %>%
    ggplot(
        aes(
            x = resolution, y = reproductiveNumber,
            fill = interaction(interval, treePrior)
        )
    ) +
    geom_violin(alpha = 0.4, scale = "width", na.rm = TRUE, draw_quantiles = c(0.025, 0.5, 0.975)) +
    scale_discrete_manual(
        aesthetics = c("colour", "fill"),
        limits = c(
            "reproductiveNumber.BD",
            "reproductiveNumber.CE",
            "reproductiveNumber.1.BD",
            "reproductiveNumber.2.BD"
        ),
        labels = c(
            ~ italic(R[0]) ~ Birth ~ Death,
            ~ italic(R[0]) ~ Coalescent ~ Exponential,
            ~ italic(R[e[1]]) ~ Birth ~ Death,
            ~ italic(R[e[2]]) ~ Birth ~ Death
        ),
        values = c(
            "reproductiveNumber.BD" = "dodgerblue",
            "reproductiveNumber.CE" = "darkorange",
            "reproductiveNumber.1.BD" = "purple",
            "reproductiveNumber.2.BD" = "green"
        )
        
    ) +
    scale_y_continuous(
        trans = "log",
        breaks = c(0.5, 1.1, 2, 3, 5, 10, 15)
    ) +
    facet_wrap(~organism, scales = "free", labeller = label_parsed, nrow = 1) +
    labs(y = TeX("\\textit{$R_{\\bullet}$}"), x = "Date resolution") +
    theme_bw() +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 14)
    )

# Combine the plots into a single panel
panel <- plot_grid(
    clock_plot, age_plot_covid_ce_filtered, reproductive_plot,
    nrow = 3, labels = "AUTO"
)

ggsave(
    plot = panel, 
    "figures/empirical_parms.pdf",
    width = 8, height = 8,
    units = "in", dpi = 300
)

panel <- plot_grid(
    clock_plot, age_plot_covid_ce_unfiltered, reproductive_plot,
    nrow = 3, labels = "AUTO"
)

ggsave(
    plot = panel, 
    "figures/empirical_parms_covid_ce_unfiltered.pdf",
    width = 8, height = 8,
    units = "in", dpi = 300
)


## Table of mean posterior estimates and HPDs
mean_HPD <- function(vec) {
    mean <- format(mean(vec), digits = 3, scientific = TRUE)
    hpd <- format(
        quantile(vec, probs = c(0.025, 0.975), na.rm = TRUE),
        digits = 2,
        scientific = TRUE
    )

    str <- paste0(mean, " (", paste(hpd, collapse = ", "), ")")

    return(str)
}


tab <- traces %>%
    select(
        organism, treePrior, resolution,
        reproductiveNumber, reproductiveNumber.1, reproductiveNumber.2,
        Tree.height, clockRate
    ) %>%
    group_by(organism, treePrior, resolution) %>%
    summarise(
        clock = mean_HPD(clockRate),
        reproductiveNumber = mean_HPD(reproductiveNumber),
        reproductiveNumber.1 = mean_HPD(reproductiveNumber.1),
        reproductiveNumber.2 = mean_HPD(reproductiveNumber.2),
        age = mean_HPD(Tree.height)
    ) %>%
    mutate(across(everything(), ~ {
        gsub(.x, pattern = "NA [(]NA, NA[)]", replacement = "-") %>%
        gsub(pattern = "[e][+]00", replacement = "") %>%
        gsub(pattern = "[e][-]0", replacement = "e-")
    }))

clock_origin_tab <- tab %>% select(resolution, clock, age)
re_tab <- tab %>% select(resolution, starts_with("reproductiveNumber"))

latex_tab1 <- kable(clock_origin_tab, format = "latex", booktabs = TRUE)
latex_tab2 <- kable(re_tab, format = "latex", booktabs = TRUE)

writeLines(latex_tab1, "empirical_clock_age.tex")
writeLines(latex_tab2, "empirical_re.tex")
  
## Emprirical densitree plots

# List .trees files in the directory
tree_files <- list.files("empirical_data", pattern = "\\.trees$", full.names = TRUE)

# Read and combine the .trees files into lists
trees <- lapply(tree_files, read.nexus)

# Remove 10% burn-in and sample 100 trees from each posterior
trees <- lapply(trees, function(x) x[sample((0.1 * length(x)):length(x), 100)])

# Name each posterior by file name without path and extension
names(trees) <- basename(tree_files) %>%
    gsub(pattern = "\\.trees", replacement = "", x = .)

# Create a vector of modified names from trees
tree_names <- names(trees) %>%
    str_remove(".trees") %>%
    str_remove("-Fixed") %>%
    str_remove("SC_") %>%
    str_replace_all("_", "~") %>%
    str_replace("saureus", "italic(S.~aureus)") %>%
    str_replace("tb", "italic(M.~tuberculosis)") %>%
    str_replace("h1n1", "H1N1") %>%
    str_replace("sars-cov-2", "SARS-CoV-2")

# colour vector matching order of colouring for resolution and tree prior above
tree_colour <- c(
    rep("dodgerblue", 3), rep("darkorange", 3),
    rep("dodgerblue", 3), rep("darkorange", 3),
    rep("dodgerblue", 3), rep("dodgerblue", 3)
)

# lower x limits for each organism
x_limits <- c(rep(-0.5, 3), rep(-0.6, 3), rep(-0.2, 3), rep(-0.65, 3), rep(-35, 3), rep(-35, 3))

tree_plot <- lapply(seq_along(trees), function(i) {
    ggdensitree(
        trees[[i]],
        colour = tree_colour[i],
        alpha = 0.1
    ) +
    ggtitle(parse(text = tree_names[i])) +
    xlim(x_limits[i], 0) +
    theme_tree2() +
    ggtree::theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
    )
})

# Combine the tree plots into a single panel
tree_panel <- cowplot::plot_grid(
    plotlist = tree_plot, ncol = 6,
    labels = "AUTO", byrow = FALSE
)

# Save the final plot as a PDF
ggsave(
    plot = tree_panel, filename = "../figures/empirical_densitrees.pdf",
    width = 17, height = 10, units = "in", dpi = 300
)

## Empirical likelihood plots
emp_likelihood_plot <- traces %>% 
    group_by(organism, resolution, treePrior) %>%
    pivot_longer(
        matches("BDSKY_Serial|CoalescentExponential"),
        names_to = "beast_phylodynamic_likelihood_name",
        values_to = "phylodynamicLikelihood",
        values_drop_na = TRUE
    ) %>%
    mutate(phylogeneticLikelihood = treeLikelihood) %>%
    mutate(treePrior = case_when(
        treePrior == "CE" ~ "Coalescent~Exponential",
        treePrior == "BD" ~ "Birth~Death"
    )) %>%
    ggplot(
        aes(
            x = phylogeneticLikelihood,
            y = phylodynamicLikelihood,
            fill = resolution,
            col = resolution
        )) +
    stat_ellipse() +
    geom_point(shape = 21, size = 3, alpha = 0.1) +
    facet_wrap(
        ~ organism * treePrior,
        scales = "free",
        labeller = label_parsed,
        dir = "v", ncol = 3
    ) +
    labs(
        x = "Phylogenetic Likelihood",
        y = "Phylodynamic Likelihood",
    ) +
    scale_discrete_manual(
        aesthetics = c("fill", "col"),
        values = c("red", "dodgerblue", "black")
    ) +
    theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 24)
    )

ggsave(
    plot = emp_likelihood_plot,
    filename = "figures/empirical_likelihood.pdf",
    width = 17, height = 10, units = "in", dpi = 300
)

## Checking posterior population size under CE
epop_plot <- traces %>% 
    group_by(organism, resolution, treePrior) %>%
    filter(treePrior == "CE") %>%
    ggplot(
        aes(
            x = ePopSize,
            y = growthRate,
            fill = resolution,
            col = resolution
        )) +
    stat_ellipse() +
    geom_point(shape = 21, size = 3, alpha = 0.1) +
    facet_wrap(
        ~ organism, scales = "free_y"
    ) +
    scale_x_log10() +
    labs(
        x = "Scaled Effective Population Size",
        y = "Growth Rate",
    ) +
    scale_discrete_manual(
        aesthetics = c("fill", "col"),
        values = c("red", "dodgerblue", "black")
    ) +
    theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 12)
    )

ggsave(
    plot = epop_plot,
    filename = "figures/empirical_effpop.pdf",
    dpi = 300, width = 6, height = 3
)
