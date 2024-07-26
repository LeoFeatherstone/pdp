###################################################################
###### Central idea, figure on plane as for Biek et al 2015 #######
###################################################################

library(tidyverse)
library(latex2exp)
library(scales)

# define log7 transform
log_7 <- function(x) log(x) / log(7)
pow_7 <- function(x) exp(7)

# write data manually for now
data <- data.frame(
    "rate" = c(0.004, 0.001, 0.000001, 0.00000001),
    "g_length" = c(13158, 29903, 2900000, 4300000),
    "name" = c(
        "H1N1",
        "SARS-CoV-2",
        "italic(S.~aureus)",
        "italic(M.~tuberculosis)"
    )
)
data <- data %>%
    mutate(m_t = 365.25 / (rate * g_length)) %>%
    slice(rep(seq_len(n()), 3))

data <- cbind.data.frame(data, "resolution" = rep(c(7, 30, 365.25)))

data <- data %>%
    mutate(m_t = log_7(m_t)) %>%
    mutate(resolution = log_7(resolution)) %>%
    mutate(col = case_when(
        !(name %in% c(
            "H1N1", "SARS-CoV-2",
            "italic(S.~aureus)", "italic(M.~tuberculosis)"
        )) ~ "dodgerblue",
        .default = "coral"
    ))


t_breaks <- sapply(c(1, 30, 365.25), function(x) log_7(x))

data %>%
    ggplot(
        aes(x = resolution, y = m_t, group = name, label = name, fill = col)
    ) +
    annotate("polygon", x = c(0, 0, 5.1), y = c(0, 5.1, 5.1), fill = "grey") +
    geom_line() +
    geom_point(pch = 21, size = 2) +
    geom_text(
        data = subset(data, resolution == 1),
        aes(x = resolution, y = m_t, label = name),
        parse = TRUE,
        nudge_y = c(-0.25, 0.25, 0.25, 0.25),
        nudge_x = c(-0.05, -0.05, -0.05, 0.75),
        size = 3,
        hjust = 0,
        check_overlap = TRUE
    ) +
    geom_abline(intercept = 1e-100, slope = 1) +
    coord_cartesian() +
    scale_x_continuous(
        breaks = c(0, 1, log_7(30), log_7(365.25), log_7(3652.5)),
        labels = c(
            "Days", "Weeks", "Months", "Years", "Decades"
        ),
        limits = c(0, 5.1),
        name = "Sampling date resolution",
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        breaks = c(0, 1, log_7(30), log_7(365.25), log_7(3652.5)),
        labels = c(
            "", "Weeks", "Months", "Years", "Decades"
        ),
        limits = c(0, 5.1),
        name = "Average substitution time",
        expand = c(0, 0)
    ) +
    annotate(
        "label",
        label = "Robust%<->%Biased",
        x = 3.2, y = 3.3,
        size = 3,
        parse = TRUE
    ) +
    coord_fixed(ratio = 1) +
    theme_classic() +
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        text = element_text(size = 12),
        axis.text = element_text(angle = 45, hjust = 1)
    )

ggsave(
    "figures/plane.pdf",
    dpi = 300,
    units = "in",
    width = 3,
    height = 3
)


