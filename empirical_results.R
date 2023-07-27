#########################################
### Wrangle and compare poterior date ###
#########################################

library(tidyverse)
library(coda)
library(xtable)

## function applies burnin
# df is data frame
# pc is percent burnin
apply_burnin <- function(df, pc) {
    return(df[-c(0:ceiling((pc / 100) * dim(df)[1])), ])
}

##### Empirical Data #####
emp_log <- dir(path = "./empirical_data", pattern = ".+[.]log")

emp_data <- lapply(
    emp_log,
    function(x) read.table(paste0("./empirical_data/", x), header = TRUE)
)
emp_data <- lapply(
    emp_data,
    function(x) apply_burnin(x, 10)
)
names(emp_data) <- gsub(
    emp_log,
    pattern = "[.]log",
    replacement = ""
)

# patch names with clock until rerun
names(emp_data) <- sapply(
    names(emp_data),
    function(x) {
        if (!grepl("_SC_BD", x) && !grepl("_RC_BD", x)) {
            return(gsub(x, pattern = "_BD", replacement = "_SC_BD"))
        }
        return(x)
    }
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
        pattern = "becomeUninfectiousRate_BDSKY_Serial",
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
                pattern = "R.|delta|origin|clock|p"
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
# capitalise names
emp_data <- emp_data %>%
    mutate(organism = case_when(
        organism == "h1n1" ~ "H1N1",
        organism == "sars-cov-2" ~ "SARS-CoV-2",
        organism == "shigella" ~ "Shigella",
        organism == "tb" ~ "TB"
    ))

## plots
# clock rate plot
pdf(file = "empirical_clock_trajectory.pdf", useDingbats = TRUE)
    emp_data %>%
        select(resolution, organism, clockRate) %>%
        group_by(resolution, organism) %>%
        ggplot() +
        geom_violin(
            aes(x = resolution, y = clockRate, fill = organism),
            alpha = 0.5,
            draw_quantiles = c(0.025, 0.5, 0.975),
        ) +
        scale_y_continuous(
            trans = "log10",
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        facet_wrap(~organism, scales = "freH1N1e_y") +
        ylab("Posterior substitution rate") +
        xlab("Date resolution") +
        theme_minimal()
dev.off()

# summary table with mean and 95% HPD for each parameter
emp_table <- emp_data %>%
    select(
        resolution,
        organism,
        clockRate,
        p,
        delta,
        R0,
        Re1,
        Re2,
        origin
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
        originHPD = paste0("[", paste(signif(quantile(origin, na.rm  = TRUE, probs = c(0.025, 0.975)),3), sep="", collapse=", "), "]"),
        .groups = "drop"
    )
print(xtable(emp_table, type = "latex", digits = 3), file = "empirical_data_table.tex")



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

pdf("empirical_origin.pdf", useDingbats = FALSE)
    emp_data %>%
    ggplot() +
    geom_histogram(
        aes(
            x = origin,
            fill = resolution
    ),
        bins = 100,
        alpha = 0.5,
        position = "identity"
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        oob = scales::censor
    ) +
    facet_wrap(~organism, scales = "free") +
    xlab("Origin") + ylab("Posterior Frequency") +
    theme_minimal()
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

