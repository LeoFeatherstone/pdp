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
    "rate" = c(0.004, 0.001, 0.0000006, 0.00000001),
    "g_length" = c(13158, 29903, 4825265, 4300000),
    "name" = c(
        "H1N1",
        "SARS-CoV-2",
        "italic(S.~sonnei)",
        "italic(M.~tuberculosis)"
    )
)
data <- data %>%
    mutate(m_t = 365.25 / (rate * g_length)) %>%
    slice(rep(seq_len(n()), 3))

data <- cbind.data.frame(data, "resolution" = rep(c(7, 30, 365.25)))

data <- data %>%
    mutate(m_t = log_7(m_t)) %>%
    mutate(resolution = log_7(resolution))


t_breaks <- sapply(c(1, 30, 365.25), function(x) log_7(x))

data %>%
    ggplot(aes(x = resolution, y = m_t, group = name, label = name)) +
    annotate("polygon", x = c(0, 0, 4.8), y = c(0, 4.8, 4.8), fill = "grey") +
    geom_line() +
    geom_point(pch = 21, size = 7, fill = alpha("red", 0.7)) +
    geom_text(
        data = subset(data, resolution == 1),
        aes(x = resolution, y = m_t, label = name),
        parse = TRUE,
        nudge_y = 0.1,
        nudge_x = -0.05,
        size = 6,
        hjust = 0
    ) +
    geom_abline(intercept = 1e-100, slope = 1) +
    coord_cartesian(xlim = c(0, 4.7), ylim = c(0, 4.7)) +
    scale_x_continuous(
        breaks = c(0, log_7(30), log_7(365.25)),
        labels = c(
            "Days", "Months", "Years"
        ),
        limits = c(0, 4.8),
        name = "Resolution",
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        breaks = c(0, log_7(30), log_7(365.25)),
        labels = c(
            "", "Months", "Years"
        ),
        limits = c(0, 4.8),
        name = "Average mutation time",
        expand = c(0, 0)
    ) +
    annotate(
        "label",
        label = "More Robust <---> Less Robust",
        x = 3, y = 3,
        size = 6
    ) +
    coord_fixed(ratio = 1) +
    theme_linedraw() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)
    )
ggsave(
    "plane.pdf",
    dpi = 300,
    units = "in"
)
