###################################################
############# Emprirical data plots ###############
###################################################

library(tracerer)
library(tidyverse)
library(latex2exp)
library(cowplot)
conflicted::conflicts_prefer(dplyr::filter())

# Get the list of .log files in the directory
log_files <- list.files("empirical_data", pattern = "\\.log$", full.names = TRUE)

# Read and combine the .log files into lists
traces <- lapply(log_files, tracerer::parse_beast_tracelog_file)
names(traces) <- basename(log_files)

# Merge lists into a single data frame
traces <- traces %>%
  bind_rows(.id  = "file") %>%
  separate_wider_delim(
    file, delim = "_", names = c("organism", "treePrior", "clock", "resolution"),
  ) %>%
  select(!matches("rate.+|freq.+")) %>%
  filter(Sample > 25000000) %>%
  mutate(resolution = gsub(pattern = "-Fixed[.]log", replacement = "", resolution)) %>%
  mutate(reproductiveNumber = case_when(
        treePrior == "CE" & organism == "h1n1" ~ (growthRate  / 91.3125) + 1,
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

# Plot posterior clock rate. Facet by organism and colour by resolution
clock_plot <- traces %>% 
    ggplot(aes(x = resolution, y = clockRate, fill = treePrior)) +
    geom_violin(alpha = 0.5) +
    scale_discrete_manual(
        aesthetics = c("colour", "fill"),
        labels = c(Birth~Death, Coalescent~Exponential),
        values = c("BD" = "dodgerblue", "CE" = "darkorange")
    ) +
    scale_y_continuous(
        trans = "log10",
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    facet_wrap(~ organism, scales = "free", labeller = label_parsed, nrow = 1) +
    ylab(TeX("Substitution rate \\textit{(subs/site/time)}")) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank()
    )

# Plot posterior tree height. Facet by organism and colour by resolution
age_plot <- traces %>% 
    filter(!(organism == "SARS-CoV-2" & treePrior == "CE")) %>% # Posterior is flat
    ggplot(aes(x = resolution, y = Tree.height, fill = treePrior)) +
    geom_violin(alpha = 0.5) +
    scale_discrete_manual(
        aesthetics = c("colour", "fill"),
        labels = c(Birth~Death, Coalescent~Exponential),
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
    facet_wrap(~ organism, scales = "free", labeller = label_parsed, nrow = 1) +
    ylab(TeX("Outbreak age (months or years)")) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank()
    )

# Plot the posterior reproductive number. Facet by organism and colour by resolution
reproductive_plot <- traces %>% 
    filter(!(organism %in% c("H1N1", "SARS-CoV-2") & resolution == "Year")) %>%
    pivot_longer(cols = starts_with("reproductiveNumber"), names_to = "interval", values_to = "reproductiveNumber") %>%
    ggplot(
        aes(
            x = resolution, y = reproductiveNumber,
            fill = interaction(interval, treePrior)
        )
    ) +
    geom_violin(alpha = 0.4, scale = "width", na.rm = TRUE) +
    scale_discrete_manual(
        aesthetics = c("colour", "fill"),
        labels = c(
            ~italic(R[0])~Birth~Death, ~italic(R[e[1]])~Birth~Death,
            ~italic(R[e[2]])~Birth~Death, ~italic(R[0])~Coalescent~Exponential
         ),
         values = c(
            "reproductiveNumber.BD" = "dodgerblue", "reproductiveNumber.CE" = "darkorange",
            "reproductiveNumber.1.BD" = "purple", "reproductiveNumber.2.BD" = "green"
        )
    ) +
    scale_y_continuous(
        trans = "log",
        breaks = c(0.5, 1, 1.1, 1.2, 2, 3, 5, 10, 15)
    ) +
    facet_wrap(~ organism, scales = "free", labeller = label_parsed, nrow = 1) +
    labs(y = TeX("\\textit{$R_{\\bullet}$}"), x = "Date resolution") +
    theme_bw() +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank()
    )

# Combine the plots into a single panel
panel <- plot_grid(
    clock_plot, age_plot, reproductive_plot,
    nrow = 3, labels = "AUTO"
)

ggsave(plot = panel, "empirical_parms.pdf", width = 10, height = 10, units = "in", dpi = 300)

