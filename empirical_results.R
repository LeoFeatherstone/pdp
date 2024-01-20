#########################################
### Wrangle and compare poterior data ###
#########################################

library(tidyverse)
library(coda)
library(xtable)
library(latex2exp)
library(scales)
library(cowplot)
library(ggExtra)
library(ggridges)

## function applies burnin
# df is data frame
# pc is percent burnin
apply_burnin <- function(df, states) {
    return(df[-which(df$Sample <= states), ])
}

##### Empirical Data #####
load("./empirical_data/posteriors.RData")

emp_data <- lapply(
    data,
    function(x) apply_burnin(x, 1e7)
)

# burnin
emp_ess <- lapply(
    emp_data,
    function(x) {
        effectiveSize(
            as.mcmc(x[, -grep(pattern = "Sample", colnames(x))])
        )
    }
)

# Find which treatments had ess < 200
low_ess <- lapply(
    seq_along(emp_ess),
    function(x) {
        if (any(emp_ess[[x]] < 200)) {
            return(
                paste(
                    names(which(emp_ess[[x]] < 200)),
                    emp_ess[[x]][which(emp_ess[[x]] < 200)]
                )
            )
        }
    }
)
names(low_ess) <- names(emp_data)
# ESS < 200 only for Tb with relaxed clock and h1n1 with year resolution.
# We don't report on both, so no issue.

# format variables names
rename_cols <- function(vec) {
    vec <- gsub(
        vec,
        pattern = "reproductiveNumber_BDSKY_Serial[.]|reproductiveNumber[.]",
        replacement = "Re"
    )
    vec <- gsub(
        vec,
        pattern = "reproductiveNumber_BDSKY_Serial$|reproductiveNumber$",
        replacement = "R0"
    )
    vec <- gsub(
        vec,
        pattern = "becomeUninfectiousRate_BDSKY_Serial|becomeUninfectiousRate",
        replacement = "delta"
    )
    vec <- gsub(
        vec,
        pattern = "samplingProportion_BDSKY_Serial",
        replacement = "p"
    )
    vec <- gsub(
        vec,
        pattern = "clockRate[.]c",
        replacement = "clockRate"
    )
    vec <- gsub(
        vec,
        pattern = "origin.+",
        replacement = "origin"
    )
    vec <- gsub(
        vec,
        pattern = "TreeHeight|Tree.height",
        replacement = "age"
    )

    return(vec)
}

emp_data <- lapply(
    emp_data,
    function(x) {
        colnames(x) <- rename_cols(colnames(x))
        return(x)
    }
)

# get parms of interest
emp_data <- lapply(
    emp_data,
    function(x) {
        return(
            x[
                ,
                grep(
                    colnames(x),
                    pattern = "R.|delta|origin|clock|p|age|growthRate|BDSKY_Serial|CoalescentExponential|likelihood"
                )
            ]
        )
    }
)

# bind data
emp_data <-
    bind_rows(emp_data, .id = "id") %>%
    separate_wider_delim(
        id,
        delim = "_",
        names = c("organism", "clock", "treePrior", "resolution")
    )

## Plotting
# Likelihood plane
p <- emp_data %>%
    subset(treePrior == "BD") %>%
    subset(!(organism == "h1n1" & resolution == "Year") & !(organism == "shigella")) %>%
    select(resolution, organism, clockRate, BDSKY_Serial, likelihood) %>%
    group_by(resolution, organism) %>%
    ggplot(aes(x = likelihood, y = BDSKY_Serial, col = resolution)) +
    geom_point(shape = 16, alpha = 0.2) +
    facet_wrap(~organism, scales = "free", nrow = 3) +
    scale_color_manual(
        values = ggsci::pal_cosmic("signature_substitutions", alpha = 0.2)(4),
        labels = parse_format(),
        name = "",
    ) +
    ylab("Phylodynamic Likelihood") +
    xlab("Phylogenetic Likelihood") +
    theme_minimal()
ggsave("likelihood_plane.pdf", dpi = 300)

# Rate and origin
emp_data %>%
    subset(treePrior == "BD") %>%
    subset(!(organism == "h1n1" & resolution == "Year") & !(organism == "shigella")) %>%
    select(resolution, organism, clockRate, origin) %>%
    group_by(resolution, organism) %>%
    ggplot(aes(x = origin, y = clockRate, col = resolution)) +
    geom_point(shape = 16, alpha = 0.2) +
    facet_wrap(~organism, scales = "free", nrow = 3) +
    scale_color_manual(
        values = ggsci::pal_cosmic("signature_substitutions", alpha = 0.2)(4),
        labels = parse_format(),
        name = "",
    ) +
    #ylab("Phylodynamic Likelihood") +
    #xlab("Phylogenetic Likelihood") +
    theme_minimal()
ggsave("empirical_rate_origin.pdf", dpi = 300)

# Transmission and origin
emp_data %>%
    pivot_longer(
        cols = c("R0", "Re1", "Re2"),
        names_to = "reproductiveNumber",
        values_to = "R",
    ) %>%
    subset(treePrior == "BD") %>%
    subset(!(organism == "h1n1" & resolution == "Year") & !(organism == "shigella")) %>%
    select(resolution, organism, R, reproductiveNumber, origin, clockRate) %>%
    group_by(resolution, organism) %>%
    ggplot(aes(x = clockRate, y = R, col = resolution, shape = reproductiveNumber)) +
    geom_point() +
    facet_wrap(~organism, scales = "free", nrow = 3) +
    scale_color_manual(
        values = ggsci::pal_cosmic("signature_substitutions", alpha = 0.2)(4),
        labels = parse_format(),
        name = "",
    ) +
    #ylab("Phylodynamic Likelihood") +
    #xlab("Phylogenetic Likelihood") +
    theme_minimal()
ggsave("empirical_transmission_origin.pdf", dpi = 300)

## Attempt panel
# h1n1
    # likelihood
    p <- emp_data %>%
    subset(treePrior == "BD" & organism == "h1n1" & resolution == "Year") %>%
    select(resolution, organism, clockRate, BDSKY_Serial, likelihood) %>%
    group_by(resolution) %>%
    ggplot(aes(x = likelihood, y = BDSKY_Serial, col = resolution)) +
    geom_point(shape = 16, alpha = 0.2) +
    scale_color_manual(
        values = ggsci::pal_cosmic("signature_substitutions", alpha = 0.2)(4),
        labels = parse_format(),
        name = "",
    ) +
    ylab("Phylodynamic Likelihood") +
    xlab("Phylogenetic Likelihood") +
    theme_minimal()
    p <- ggMarginal(p, type = "histogram", groupColour = TRUE, groupFill = TRUE)

pdf("h1n1_test_marginal.pdf", useDingbats = FALSE)
    ggMarginal(p, type = "histogram", groupColour = TRUE, groupFill = TRUE)
dev.off()

    # rate and origin
    h1n1_rate_marginal <- emp_data %>%
    subset(treePrior == "BD" & organism == "h1n1" & !(resolution == "Year")) %>%
    select(resolution, organism, clockRate, origin, likelihood) %>%
    group_by(resolution) %>%
    ggplot(aes(x = origin, y = clockRate, col = resolution)) +
    geom_point(shape = 16, alpha = 0.2) +
    scale_color_manual(
        values = ggsci::pal_cosmic("signature_substitutions", alpha = 0.2)(4),
        labels = parse_format(),
        name = "",
    ) +
    #ylab("Phylodynamic Likelihood") +
    #xlab("Phylogenetic Likelihood") +
    theme_minimal()
    h1n1_rate_marginal <- ggMarginal(h1n1_rate_marginal, type = "histogram", groupColour = TRUE, groupFill = TRUE)


cowplot::plot_grid(p, h1n1_rate_marginal, nrow = 1)
ggsave("test_panel.pdf", dpi = 300, width = 6, height = 4, units = "in")
    # transmission and origin


# sampling and origin plot
emp_aln <- paste0(
    "./empirical_data/",
    c("h1n1", "sars-cov-2", "tb"),
    ".fasta"
)
emp_aln <- lapply(
    emp_aln,
    function(x) ape::read.dna(x, format = "fasta")
)
names(emp_aln) <- c("h1n1", "sars-cov-2", "tb")

samp_times <- lapply(
    emp_aln,
    function(x) {
        rownames(x)
    }
)
samp_times <- lapply(
    seq_along(samp_times),
    function(i) {
        as.data.frame(cbind(
            samp_times[[i]]
        ))
    }
)
names(samp_times) <- names(emp_aln)

samp_times <- bind_rows(samp_times, .id = "id")
colnames(samp_times) <- c("organism", "label")
samp_times <- samp_times %>%
    separate(
        label,
        sep = "_",
        into = c("name", "Year", "Month", "Day")
    ) %>%
    mutate(
        Year = (as.Date(Year)),
        Month = (as.Date(Month)),
        Day = (as.Date(Day))
    ) %>%
    pivot_longer(
        c("Year", "Month", "Day"),
        names_to = "resolution",
        values_to = "date"
    ) %>%
    group_by(organism, resolution) %>%
    mutate(
        mrsd = max(date)
    )
# associate origin data
max_dates <- samp_times %>%
    select(organism, resolution, mrsd) %>%
    distinct()

origin_times <- emp_data %>%
    select(origin, organism, resolution) %>%
    left_join(
        max_dates,
        by = c("organism", "resolution")
    ) %>%
    group_by(organism, resolution) %>%
    mutate(
        hpd_upper = quantile(origin, 0.975, na.rm = TRUE),
        hpd_lower = quantile(origin, 0.025, na.rm = TRUE)
    ) %>%
    filter(origin >= hpd_lower & origin <= hpd_upper) %>%
    ungroup() %>%
    group_by(organism) %>%
    mutate(
        adjusted_origin = as.Date(date_decimal((decimal_date(max(mrsd)) - origin)))
    ) %>%
    ungroup() %>%
    group_by(organism, resolution) %>%
    mutate(
        abs_origin = as.Date(date_decimal((decimal_date(mrsd) - origin)))
    )

# filter
samp_time_filtered <- samp_times %>%
    subset(
        !(organism == "h1n1" & resolution == "Year")
        &
        !(organism == "sars-cov-2" & resolution == "Year")
    )
origin_times_filtered <- origin_times %>%
    subset(
        !(organism == "h1n1" & resolution == "Year")
        &
        !(organism == "sars-cov-2" & resolution == "Year")
    )

ggplot() +
    geom_density_ridges(
        stat = "binline",
        bins = 40,
        draw_baseline = FALSE,
        data = samp_time_filtered,
        scale = 1,
        jittered_points = TRUE,
        position = position_points_jitter(width = 0.05, height = 0),
        point_shape = "|", point_size = 3, point_alpha = 1,
        aes(x = date, y = resolution, fill = resolution, group = resolution)
    ) +
    geom_density_ridges(
        stat = "binline",
        data = origin_times_filtered,
        draw_baseline = FALSE,
        scale = 1,
        bins = 50,
        fill = "grey",
        aes(x = adjusted_origin, y = resolution, group = resolution)
    ) +
    scale_fill_manual(
        values = ggsci::pal_cosmic("signature_substitutions", alpha = 0.6)(4),
        labels = parse_format(),
        name = ""
    ) +
    scale_x_date(date_labels = "%b-%y") +
    facet_wrap(~organism, scales = "free", ncol = 1) +
    theme_minimal()
ggsave("adjusted_origin_sampling.pdf", dpi = 300, width = 4, height = 6, units = "in")
########## Older Plots ##########



## plots
emp_data$organism <- factor(emp_data$organism, labels = c(
    h1n1 = "H1N1",
    `sars-cov-2` = "SARS-CoV-2",
    shigella = TeX("\\textit{S. sonnei}"),
    tb = TeX("\\textit{M. tuberculosis}")
))

emp_data <- emp_data %>% mutate(organism = case_when(
    organism == "h1n1" ~ "H1N1",
    organism == "sars-cov-2" ~ "SARS-CoV-2",
    organism == "shigella" ~ as.character(TeX("\\textit{S. sonnei}")),
    organism == "tb" ~ as.character(TeX("\\textit{M. tuberculosis}"))
))

# clock rate plot
pdf(file = "empirical_clock.pdf", useDingbats = TRUE)
emp_data %>%
    subset(treePrior == "BD") %>%
    subset(!(organism == "h1n1" & resolution == "Year")) %>%
    select(resolution, organism, clockRate) %>%
    group_by(resolution, organism) %>%
    ggplot(aes(x = resolution, y = clockRate, fill = organism)) +
    geom_violin(
        alpha = 0.5,
        draw_quantiles = c(0.025, 0.5, 0.975),
    ) +
    scale_y_continuous(
        trans = "log10",
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_fill_manual(
        values = ggsci::pal_lancet("lanonc", alpha = 0.6)(4),
        labels = parse_format(),
        name = "",
    ) +
    facet_wrap(
        ~organism,
        scales = "free_y",
        labeller = label_parsed
    ) +
    ylab("Posterior substitution rate") +
    xlab("Date resolution") +
    theme_minimal() +
    theme(
        legend.text.align = 0,
        legend.position = "botton"
    )
dev.off()

# age of outbreak
pdf("empirical_age.pdf", useDingbats = FALSE)
emp_data %>%
    subset(treePrior == "BD") %>%
    # subset(!(organism == "H1N1" & resolution == "Year")) %>%
    # subset(!(organism == "SARS-CoV-2" & resolution == "Year")) %>%
    select(resolution, organism, age) %>%
    group_by(resolution, organism) %>%
    ggplot(aes(x = 365.25 * age, fill = resolution)) +
    geom_histogram(
        bins = 100,
        alpha = 0.5,
        position = "identity"
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        oob = scales::censor
    ) +
    scale_fill_manual(
        values = ggsci::pal_lancet("lanonc", alpha = 0.6)(3),
        labels = parse_format(),
        name = ""
    ) +
    scale_x_continuous(
        breaks = c(1, 7, 30, 60, 90, 180, 365.25, 913.125, 1095.75, 3652.5, 9131.25, 18262.5),
        labels = c("Day", "Week", "Month", "2 Months", "3 Months", "6 Months", "Year", "2.5 Years", "3 Years", "10 Years", "25 years", "50 Years")
    ) +
    facet_wrap(
        ~organism,
        scales = "free",
        labeller = label_parsed
    ) +
    ylab("Posterior Frequency outbreak age") +
    xlab("Time ") +
    theme_minimal() +
    theme(
        legend.text.align = 0,
        legend.position = "right",
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.minor = element_blank()
    )

dev.off()

## Reproductive numbers
# sampling times for run
# samp dates for rug
emp_aln <- paste0(
    "./empirical_data/",
    c("h1n1", "sars-cov-2", "shigella", "tb"),
    ".fasta"
)
emp_aln <- lapply(
    emp_aln,
    function(x) ape::read.dna(x, format = "fasta")
)
names(emp_aln) <- c("h1n1", "sars-cov-2", "shigella", "tb")

samp_times <- lapply(
    emp_aln,
    function(x) {
        rownames(x)
    }
)
samp_times <- lapply(
    seq_along(samp_times),
    function(i) {
        as.data.frame(cbind(
            samp_times[[i]]
        ))
    }
)
names(samp_times) <- names(emp_aln)

samp_times <- bind_rows(samp_times, .id = "id")
colnames(samp_times) <- c("organism", "label")
samp_times <- samp_times %>%
    separate(
        label,
        sep = "_",
        into = c("name", "Year", "Month", "Day")
    ) %>%
    mutate(
        Year = decimal_date(as.Date(Year)),
        Month = decimal_date(as.Date(Month)),
        Day = decimal_date(as.Date(Day))
    ) %>%
    pivot_longer(
        c("Year", "Month", "Day"),
        names_to = "resolution",
        values_to = "date"
    )

max_dates <- samp_times %>%
    group_by(organism, resolution) %>%
    summarise(max_date = max(date))

samp_times <- samp_times %>%
    mutate(date = as.Date(date_decimal(date))) %>%
    subset(
        !(organism == "h1n1" & resolution == "Year") &
        !(organism == "sars-cov-2" & resolution == "Year")
    )

emp_data_r <- left_join(
        emp_data,
        max_dates,
        by = c("organism", "resolution")
    ) %>%
    pivot_longer(
        cols = matches("^R0$|^Re1$|^Re2$"),
        names_to = "R",
        values_to = "value",
        values_drop_na = TRUE
    ) %>%
    group_by(organism, resolution) %>%
    mutate(mean_age = mean(age)) %>%
    mutate(date = case_when(
        R == "R0" ~ as.Date(date_decimal(max_date - (0.5 * mean_age))),
        R == "Re1" ~ as.Date(date_decimal(max_date - (0.75 * mean_age))),
        R == "Re2" ~ as.Date(date_decimal(max_date - (0.25 * mean_age))),
    )) %>%
    mutate(
        origin_time = as.Date(date_decimal((max_date - (age)))),
        change_time = as.Date(date_decimal((max_date - (age/2))))
    )


plot_r <- function(emp_data_r, condition, width, ylab, legend_pos) {
    library(tidyverse)
    library(ggplot2)

    emp_data_r %>%
    subset(
        !(organism == "h1n1" & resolution == "Year") &
        !(organism == "sars-cov-2" & resolution == "Year") &
        treePrior == "BD"
    ) %>%
    subset(organism == condition) %>%
    ggplot(aes(x = date, fill = resolution)) +
    geom_violin(
        aes(
            y = value,
            group = fct_cross(resolution, R),
        ),
        scale = "width",
        width = width,
        position = "dodge",
    ) +
    geom_histogram(
        aes(x = origin_time, y = 15 * ..count.. / sum(..count..)),
        position = "identity",
        bins = 50
    ) +
    geom_rug(
        data = subset(samp_times, organism == condition),
        aes(x = date, col = resolution)
    ) +
    #geom_vline(aes(xintercept = mean(change_time), col = resolution)) +
    #facet_wrap(~organism, scales = "free") +
    scale_discrete_manual(
        aesthetics = c("fill", "col"),
        values = ggsci::pal_lancet("lanonc", alpha = 0.4)(3),
        labels = parse_format(),
        name = ""
    ) +
    ylab(TeX(paste0("Posterior $", ylab, "$"))) +
    ggtitle(condition) +
    xlab("Date") +
    scale_x_date(date_labels = "%b-%Y", minor_breaks = "1 month") +
    theme_minimal() +
    theme(
        legend.text.align = 0,
        legend.position = legend_pos,
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)
    )
}

main <- plot_grid(
    plot_r(emp_data_r, "h1n1", 80, ylab = "R_{0}", "none"),
    plot_r(emp_data_r, "sars-cov-2", 30, ylab = "R_{0}", "none"),
    plot_r(emp_data_r, "shigella", 600, ylab = "R_{e_{\\bullet}}", "none"),
    plot_r(emp_data_r, "tb", 9000, ylab = "R_{e_{\\bullet}}", "none"),
    nrow = 2, labels = "AUTO"
)


pdf("empirical_R.pdf", useDingbats = FALSE)
    plot_grid(
        main,
        get_legend(plot_r(emp_data_r, "tb", 0, ylab = "", "bottom")),
        rel_heights = c(10, 1),
        ncol = 1
    )
dev.off()

# summary table with mean and 95% HPD for each parameter
get_hpd <- function(parm) {
    return(
        paste0(
            round(mean(parm), 3),
            ", ",
            "[",
            paste(
                round(
                    quantile(parm, na.rm = TRUE, probs = c(0.025, 0.975)),
                    3
                ),
                sep = "", collapse = ", "
            ),
            "]"
        )
    )
}
emp_table <- emp_data %>%
    subset(treePrior == "BD") %>%
    # mutate_if(is.numeric, round) %>%
    select(
        resolution,
        organism,
        clockRate,
        p,
        delta,
        R0,
        Re1,
        Re2,
        age
    ) %>%
    group_by(organism, resolution) %>%
    summarise(
        R0 = get_hpd(R0),
        Re1 = get_hpd(Re1),
        Re2 = get_hpd(Re2),
        p = get_hpd(p),
        delta = get_hpd(delta),
        Age = get_hpd(age),
        .groups = "drop"
    )

print(
    xtable(emp_table, type = "latex", digits = 3),
    file = "empirical_data_table.tex",
    include.rownames = FALSE
)

#############################################################################
################### Supplementary Coal' Exp' Stuff ##########################
#############################################################################

pdf("empirical_ce.pdf", useDingbats = FALSE)
    emp_data %>%
        subset(
            treePrior == "CE" & organism == "h1n1" & !(resolution == "Year")
        ) %>%
        select(growthRate, age, clockRate, resolution) %>%
        mutate(
            R0 = 1 + ((4 / 365.25) * growthRate)
        ) %>%
        pivot_longer(
            cols = c("age", "clockRate", "R0", "growthRate"),
            names_to = "parm",
            values_to = "posterior"
        ) %>%
        ggplot() +
        geom_violin(aes(y = posterior, fill = resolution, x = NA)) +
        facet_wrap(~parm, scales = "free")
dev.off()
