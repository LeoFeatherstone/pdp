#########################################
### Wrangle and compare poterior date ###
#########################################

library(tidyverse)
library(coda)
library(xtable)
library(latex2exp)
library(scales)

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
        pattern = "reproductiveNumber_BDSKY_Serial$",
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
            x[,
            grep(
                colnames(x),
                pattern = "R.|delta|origin|clock|p|age"
            )]
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

## plots
emp_data$organism <- factor(emp_data$organism, labels = c(
    h1n1 = "H1N1",
    `sars-cov-2` = "SARS-CoV-2",
    shigella = TeX("\\textit{S. sonnei}"),
    tb = TeX("\\textit{M. tuberculosis}"))
)


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
    #subset(!(organism == "H1N1" & resolution == "Year")) %>%
    #subset(!(organism == "SARS-CoV-2" & resolution == "Year")) %>%
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
    ylab("Posterior Frequency outbreak duration") +
    xlab("Time ") +
    theme_minimal() +
    theme(
        legend.text.align = 0,
        legend.position = "right",
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.minor = element_blank()
    )

dev.off()

# summary table with mean and 95% HPD for each parameter
emp_table <- emp_data %>%
    subset(treePrior == "BD") %>%
    #mutate_if(is.numeric, round) %>%
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
        meanR0 = mean(R0),
        R0HPD = paste0("[", paste(signif(quantile(R0, na.rm  = TRUE, probs = c(0.025, 0.975)),3), sep="", collapse=", "), "]"),
        meanRe1 = mean(Re1),
        Re1HPD = paste0("[", paste(signif(quantile(Re1, na.rm  = TRUE, probs = c(0.025, 0.975)),3), sep="", collapse=", "), "]"),
        meanRe2 = mean(Re2),
        Re2HPD = paste0("[", paste(signif(quantile(Re2, na.rm  = TRUE, probs = c(0.025, 0.975)),3), sep="", collapse=", "), "]"),
        meanP = mean(p),
        pHPD = paste0("[", paste(signif(quantile(p, na.rm  = TRUE, probs = c(0.025, 0.975)),3), sep="", collapse=", "), "]"),
        meanDelta = mean(delta),
        deltaHPD = paste0("[", paste(signif(quantile(delta, na.rm  = TRUE, probs = c(0.025, 0.975)),3), sep="", collapse=", "), "]"),
        meanOrigin = mean(origin),
        ageHPD = paste0("[", paste(signif(quantile(age, na.rm  = TRUE, probs = c(0.025, 0.975)),3), sep="", collapse=", "), "]"),
        .groups = "drop"
    )
print(
    xtable(emp_table, type = "latex", digits = 3),
    file = "empirical_data_table.tex",
    include.rownames = FALSE
)



# emp plots - NB below was done before formatting emp_data col names. Temporary!
cols <- c("coralred", "dodgerblue")
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
names(emp_aln) <- c("h1n1", "sars_cov_2", "shigella", "tb")
samp_times <- lapply(
    emp_aln,
    function(x) {
        decimal_date(
        as.Date(
            gsub(rownames(x), pattern = ".*_", replacement = "")
        ))
    }
)
max_samp_times <- unlist(sapply(samp_times, function(x) max(x)))

# make df
max(unlist(lapply(samp_times, function(y) length(y)))) # = 161

samp_times <- lapply(
    samp_times,
    function(x) { c(x, rep(NA, times = (161 - length(x)))) }
)
samp_times <- bind_rows(samp_times)

#h1n1
h1n1 <- emp_data %>%
    filter(
        organism == "H1N1"
        &
        resolution != "Year"
    ) %>%
    ggplot() +
    geom_violin(
        aes(
            x = (mean(as.Date(date_decimal(max_samp_times["h1n1"] - (origin * 0.5))))),
            y = R0,
            fill = resolution,
            group = resolution
        ),
        scale = "width",
        position = "identity",
        draw_quantiles = c(0.5),
        alpha = 0.5,
        width = 0.25
    ) +
    # geom_histogram(
    #     aes(
    #         x = max_samp_times["h1n1"] - origin,
    #         y = (..count.. / 5000) + 1,
    #         fill = resolution),
    #     bins = 100,
    #     alpha = 0.5,
    #     position = "identity"
    # ) +
    geom_rug(data = samp_times, aes(x = as.Date(date_decimal(h1n1)))) + 
    coord_cartesian(clip = "off") +
    coord_cartesian(ylim = c(1, 1.21)) +
    scale_y_continuous(
       expand = c(0, 0),
       oob = scales::censor
    ) +
    xlab("Date") + ylab("Reproductive Number") +
    theme_minimal()


# sars-cov-2 data
covid <- emp_data %>%
    filter(
        organism == "SARS-CoV-2"
    ) %>%
    ggplot() +
    geom_violin(
        aes(
            x = mean(max_samp_times["sars_cov_2"] - (origin * 0.5)),
            y = R0,
            fill = resolution
        ),
        scale = "width",
        position = "identity",
        draw_quantiles = c(0.5),
        alpha = 0.5,
        width = 0.05
    ) +
    # geom_histogram(
    #     aes(
    #         x = max_samp_times["sars_cov_2"] - origin,
    #         y = (..count.. / 100),
    #         fill = resolution),
    #     bins = 100,
    #     alpha = 0.5,
    #     position = "identity"
    # ) +
    geom_rug(data = samp_times, aes(x = sars_cov_2)) + 
    coord_cartesian(clip = "off") +
    xlab("Date") + ylab("Reproductive Number") +
    theme_minimal()

# tb data
tb <- emp_data %>%
    filter(
        organism == "TB"
    ) %>%
    ggplot() +
    geom_violin(
        aes(
            x = mean(max_samp_times["tb"] - (origin * 0.75)),
            y = Re1,
            fill = resolution
        ),
        scale = "area",
        position = "dodge",
        draw_quantiles = c(0.5),
        alpha = 0.5,
        width = 8
    ) +
    geom_violin(
        aes(
            x = mean(max_samp_times["tb"] - (origin * 0.25)),
            y = Re2,
            fill = resolution
        ),
        scale = "area",
        position = "dodge",
        draw_quantiles = c(0.5),
        alpha = 0.5,
        width = 8
     ) +
    # geom_histogram(
    #     aes(
    #         x = max_samp_times["tb"] - origin,
    #         y = (..count.. / 200),
    #         fill = resolution),
    #     bins = 100,
    #     alpha = 0.5,
    #     position = "identity"
    # ) +
    geom_rug(data = samp_times, aes(x = tb)) + 
    coord_cartesian(clip = "off") +
    xlim(1980, 2010) + ylim(0, 5) +
    xlab("Date") + ylab("Reproductive Number") + 
    theme_minimal()

# shigella
shigella <- emp_data %>%
    filter(
        organism == "Shigella"
    ) %>%
    ggplot() +
    geom_violin(
        aes(
            x = mean(max_samp_times["shigella"] - (origin * 0.75)),
            y = Re1,
            fill = resolution
        ),
        scale = "width",
        position = "dodge",
        draw_quantiles = c(0.5),
        alpha = 0.5,
        width = 2
    ) +
    geom_violin(
        aes(
            x = mean(max_samp_times["shigella"] - (origin * 0.25)),
            y = Re2,
            fill = resolution
        ),
        scale = "width",
        position = "dodge",
        draw_quantiles = c(0.5),
        alpha = 0.5,
        width = 2
     ) +
    # geom_histogram(
    #     aes(
    #         x = max_samp_times["shigella"] - origin,
    #         y = (..count.. / sum(..count..)) + 0.8,
    #         fill = resolution),
    #     bins = 100,
    #     alpha = 0.5,
    #     position = "identity"
    # ) +
    geom_rug(data = samp_times, aes(x = shigella)) +
    coord_cartesian(clip = "off") +
    coord_cartesian(ylim = c(0.8, 1.25)) +
    scale_y_continuous(
        expand = c(0, 0),
        oob = scales::censor
    ) +
    xlab("Date") + ylab("Reproductive Number") +
    theme_minimal()

empirical_plot <- cowplot::plot_grid(
    labels = "AUTO",
    covid, h1n1, shigella, tb,
    nrow = 2, ncol = 2
)
pdf("empirical_plot.pdf", useDingbats = FALSE)
    empirical_plot
dev.off()



## Tabular results
# mean
emp_data %>%
    group_by(organism, resolution) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE)
# precision
emp_data %>%
    group_by(organism, resolution) %>%
    summarise_if(is.numeric, var, na.rm = TRUE)

## parallel plot
norm <- function(x) {
    return((x - mean(x)) / sd(x))
}

pdf("empirical_normalised_parallel.pdf", useDingbats = FALSE)
    emp_data %>%
        filter(clock == "SC") %>%
        filter(resolution == "Day" | resolution == "Month") %>%
        group_by(organism) %>%
        mutate(
            R0 = norm(R0),
            Re1 = norm(Re1),
            Re2 = norm(Re2),
            origin = norm(origin),
            clockRate = norm(clockRate),
        ) %>%
        group_by(organism, resolution) %>%
        pivot_longer(
            matches("R0|Re1|Re2|origin|clockRate"),
            names_to = "parameter",
            values_to = "value",
            values_drop_na = TRUE
        ) %>%
        ggplot() +
        geom_line(
            aes(x = parameter, y = value, col = resolution, group = Sample),
            alpha = 0.1
        ) +
        scale_x_discrete(limits = c("origin", "clockRate", "R0", "Re1", "Re2")) +
        facet_wrap(~organism * resolution, scales = "free_x", nrow = 2)
dev.off()

