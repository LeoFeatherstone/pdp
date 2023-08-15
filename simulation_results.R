#########################################
### Wrangle and compare poterior data ###
#########################################

library(tidyverse)
library(coda)
library(tracerer)

## function applies burnin
# df is data frame
# pc is percent burnin
apply_burnin <- function(df, pc) {
    return(df[-c(0:ceiling((pc / 100) * dim(df)[1])), ])
}

# Get data
load("sim_study/posteriors.RData")
data <- lapply(
    data,
    function(x) apply_burnin(x, 30)
)


# format variables names
rename_cols <- function(vec) {
    vec <- gsub(
        vec,
        pattern = "reproductiveNumber[.]",
        replacement = "Re"
    )
    vec <- gsub(
        vec,
        pattern = "^reproductiveNumber$",
        replacement = "R0"
    )
    vec <- gsub(
        vec,
        pattern = "becomeUninfectiousRate",
        replacement = "delta"
    )
    vec <- gsub(
        vec,
        pattern = "samplingProportion",
        replacement = "p"
    )
    return(vec)
}

data <- lapply(
    data,
    function(x) {
        colnames(x) <- rename_cols(colnames(x))
        return(x)
    }
)

# get parms of interest
data <- lapply(
    data,
    function(x) {
        return(
            x[,
            grep(
                colnames(x),
                pattern = "^R.+|delta|origin|clock|^p$|Sample"
            )]
        )
    }
)

ess <- lapply(
    data,
    function(x) {
        # effectiveSize(
        #     as.mcmc(x)
        # )
        calc_esses(x, 10000)
    }
)

ess <- bind_rows(ess, .id = "id") %>%
    separate_wider_delim(
        id,
        delim = "_",
        names = c("organism", "clock", "treePrior", "resolution", "id")
    ) %>%
    group_by(organism, resolution)

ess_summary <- ess %>%
    summarise(
        clockCount = sum(clockRate < 200),
        R0Count = sum(R0 < 200),
        Re1Count = sum(Re1 < 200),
        Re2Count = sum(Re2 < 200),
        originCount = sum(origin < 200),
        pCount = sum(p < 200),
        deltaCount = sum(delta < 200),
    ) %>%
    filter()

ess %>%
    filter(
        clockRate < 200
        | origin < 200
        | R0 < 200
        | Re1 < 200
        | Re2 < 200
    ) %>%
    select(organism, resolution, id, clockRate, origin, R0, Re1, Re2) %>%
    print(n = 1000)

ess %>% filter(organism == "shigella" & id == "1")



# visualise to see where work needs to be done
ggplot(ess) +
    geom_histogram(
        aes(x = clockRate, fill = resolution),
        alpha = 0.45,
        bins = 20,
        position = "identity"
    ) +
    annotate("segment", xend = 200, x = 200, y = 0, yend = Inf) +
    facet_wrap(
        ~ organism,
        scales = "free"
    )

data %>%
    filter(organism == "SARS-CoV-2" & resolution != "Year") %>%
    ggplot() +
    geom_histogram(
        aes(x = delta, fill = resolution),
        alpha = 0.45,
        bins = 200,
        position = "identity"
    ) +
    xlim(0, 8)

# Find which treatments had ess < 200
lapply(
    seq_along(sim_ess[201:400]),
    function(x) {
        if (any(sim_ess[[x]] < 200)) {
            tmp <- sim_ess[[x]]
            print(names(sim_ess)[x])
            print(names(tmp[which(tmp < 200)]))
            print(tmp[which(tmp < 200)])

        }
    }
)

# Find which parms of interest had ess < 200
lapply(
    seq_along(sim_ess),
    function(x) {
        if (any(sim_ess[[x]] < 200)) {
            print(names(sim_ess)[x])
            print(sim_ess[x])

        }
    }
)
# bind data
sim_data <-
    bind_rows(sim_data, .id = "id") %>%
    separate_wider_delim(
        id,
        delim = "_",
        names = c("organism", "treePrior", "resolution", "replicate")
    )

sim_data <- sim_data %>%
    mutate(organism = case_when(
        organism == "h1n1" ~ "H1N1",
        organism == "sars-cov-2" ~ "SARS-CoV-2",
        organism == "shigella" ~ "Shigella",
        organism == "tb" ~ "TB"
    ))

### sim plots

## clock rate
# define true value df
rate_true <- data.frame(
    organism = c("h1n1", "sars-cov-2", "shigella", "tb"),
    rate = c(0.004, 0.001, 0.0000006, 0.00000001)
)

pdf(file = "sim_clock_trajectory.pdf", useDingbats = TRUE)
    sim_data %>%
        select(resolution, replicate, organism, clockRate) %>%
        group_by(replicate, resolution, organism) %>%
        summarise(mean_clock_rate = ((mean(clockRate)))) %>%
        ggplot() +
        geom_boxplot(
            aes(x = resolution, y = mean_clock_rate, col = organism)
        ) +
        geom_line(
            aes(
                x = resolution,
                y = mean_clock_rate,
                group = replicate,
                col = organism
            ),
            alpha = 0.5
        ) +
        geom_segment(
            data = rate_true,
            aes(y = rate, yend = rate, x = -Inf, xend = Inf),
            col = "black"
        ) +
        scale_y_continuous(
            trans = "log10",
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        facet_wrap(~organism, scales = "free_y") +
        ylab("Posterior mean  substitution rate") +
        xlab("Date resolution") +
        theme_minimal()
dev.off()

# R0/e plots

## clock rate
# define true value df
re_true <- data.frame(
    organism = c("h1n1", "sars-cov-2", "shigella", "tb"),
    R0 = c(1.3, 2.5, NA, NA),
    Re1 = c(NA, NA, 1.5, 2),
    Re2 = c(NA, NA, 1, 1.1)
) %>%
pivot_longer(
            matches("R0|Re1|Re2"),
            names_to = "interval",
            values_to = "Reff",
            values_drop_na = TRUE
        )

pdf(file = "sim_Re_trajectory.pdf", useDingbats = TRUE)
    sim_data %>%
        select(resolution, replicate, organism, R0, Re1, Re2) %>%
        pivot_longer(
            matches("R0|Re1|Re2"),
            names_to = "interval",
            values_to = "Reff",
            values_drop_na = TRUE
        ) %>%
        group_by(replicate, resolution, organism, interval) %>%
        summarise(meanR = ((mean(Reff)))) %>%
        ggplot() +
        geom_boxplot(
            aes(
                x = resolution,
                y = meanR,
                col = organism
            )
        ) +
        geom_line(
            aes(
                x = resolution,
                y = meanR,
                group = replicate,
                col = organism
            ),
            alpha = 0.5
        ) +
        geom_segment(
            data = re_true,
            aes(y = Reff, yend = Reff, x = -Inf, xend = Inf),
            col = "black"
         ) +
        scale_y_continuous(
            trans = "log10",
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        facet_wrap(~organism * interval, scales = "free_y") +
        ylab("Posterior mean  Reff") +
        xlab("Date resolution") +
        theme_minimal()
dev.off()

## origin
# define true value df
origin_true <- data.frame(
    organism = c("H1N1", "SARS-CoV-2", "Shigella", "TB"),
    rate = c(0.25, 0.16, 0.5, 25)
)

pdf(file = "sim_origin_trajectory.pdf", useDingbats = TRUE)
    sim_data %>%
        select(resolution, replicate, organism, origin) %>%
        group_by(replicate, resolution, organism) %>%
        summarise(mean_origin = ((mean(origin)))) %>%
        ggplot() +
        geom_boxplot(
            aes(x = resolution, y = mean_origin, col = organism)
        ) +
        geom_line(
            aes(
                x = resolution,
                y = mean_origin,
                group = replicate,
                col = organism
            ),
            alpha = 0.5
        ) +
        geom_segment(
            data = origin_true,
            aes(y = rate, yend = rate, x = -Inf, xend = Inf),
            col = "black"
        ) +
        facet_wrap(~organism, scales = "free_y") +
        ylab("Posterior mean  origin") +
        xlab("Date resolution") +
        theme_minimal() +
        theme(
            legend.position = "bottom"
        )
dev.off()

# summary table
sim_table <- sim_data %>%
    select(
        resolution,
        replicate,
        organism,
        clockRate,
        p,
        delta,
        R0,
        Re1,
        Re2,
        origin
    ) %>%
    group_by(resolution, replicate, organism) %>%
    summarise(
        meanR0 = mean(R0),
        meanRe1 = mean(Re1),
        meanRe2 = mean(Re2),
        meanP = mean(p),
        meanDelta = mean(delta),
        meanOrigin = mean(origin),
        .groups = "drop"
    ) %>%
    group_by(organism, resolution) %>%
    summarise(
        meanR0Err = case_when(
            organism == "sars-cov-2" ~ mean(abs(meanR0 - 2.5)),
            organism == "h1n1" ~ mean(abs(meanR0 - 1.3)),
            .default = NA
        ),
        meanRe1Err = case_when(
            organism == "shigella" ~ mean(abs(meanRe1 - 1.5)),
            organism == "tb" ~ mean(abs(meanRe1 - 2)),
            .default = NA
        ),
        meanRe2Err = case_when(
            organism == "shigella" ~ mean(abs(meanRe2 - 1)),
            organism == "tb" ~ mean(abs(meanRe2 - 1.1)),
            .default = NA
        ),
        meanPErr = case_when(
            organism == "h1n1" ~ mean(abs(meanP - 0.015)),
            organism == "sars-cov-2" ~ mean(abs(meanP - 0.8)),
            organism == "shigella" ~ mean(abs(meanP - 0.4)),
            organism == "tb" ~ mean(abs(meanP - 0.08)),
            .default = NA
        ),
        meanDeltaErr = case_when(
            organism == "h1n1" ~ mean(abs(meanDelta - 91.3125)),
            organism == "sars-cov-2" ~ mean(abs(meanDelta - 36.525)),
            organism == "shigella" ~ mean(abs(meanDelta - 52.18)),
            organism == "tb" ~ mean(abs(meanDelta - 0.125)),
            .default = NA
        ),
        meanOriginErr = case_when(
            organism == "h1n1" ~ mean(abs(meanOrigin - 0.25)),
            organism == "sars-cov-2" ~ mean(abs(meanOrigin - 0.16)),
            organism == "shigella" ~ mean(abs(meanOrigin - 0.5)),
            organism == "tb" ~ mean(abs(meanOrigin - 25)),
            .default = NA
        )
    ) %>%
    distinct()

# brilliant table to tex writer!
print(xtable(sim_table, type = "latex"), file = "sim_data_table.tex")

# figure comparing error and percent of mutation time
err_vs_res <- sim_data %>%
    select(
        resolution,
        replicate,
        organism,
        clockRate,
        p,
        delta,
        R0,
        Re1,
        Re2,
        origin
    ) %>%
    group_by(resolution, replicate, organism) %>%
    summarise(
        meanR0 = mean(R0),
        meanRe1 = mean(Re1),
        meanRe2 = mean(Re2),
        meanP = mean(p),
        meanDelta = mean(delta),
        meanOrigin = mean(origin),
        .groups = "drop"
    ) %>%
    group_by(organism, resolution) %>%
    mutate(
        R0Err = case_when(
            organism == "SARS-CoV-2" ~ ((meanR0 - 2.5)),
            organism == "H1N1" ~ ((meanR0 - 1.3)),
            .default = NA
        ),
        Re1Err = case_when(
            organism == "Shigella" ~ ((meanRe1 - 1.5)),
            organism == "TB" ~ ((meanRe1 - 2)),
            .default = NA
        ),
        Re2Err = case_when(
            organism == "Shigella" ~ ((meanRe2 - 1)),
            organism == "TB" ~ ((meanRe2 - 1.1)),
            .default = NA
        ),
        PErr = case_when(
            organism == "H1N1" ~ ((meanP - 0.015)),
            organism == "SARS-CoV-2" ~ ((meanP - 0.8)),
            organism == "Shigella" ~ ((meanP - 0.4)),
            organism == "TB" ~ ((meanP - 0.08)),
            .default = NA
        ),
        DeltaErr = case_when(
            organism == "H1N1" ~ ((meanDelta - 91.3125)),
            organism == "SARS-CoV-2" ~ ((meanDelta - 36.525)),
            organism == "Shigella" ~ mean((meanDelta - 52.18)),
            organism == "TB" ~ ((meanDelta - 0.125)),
            .default = NA
        ),
        OriginErr = case_when(
            organism == "H1N1" ~ ((meanOrigin - 0.25)),
            organism == "SARS-CoV-2" ~ ((meanOrigin - 0.16)),
            organism == "Shigella" ~ ((meanOrigin - 0.5)),
            organism == "TB" ~ ((meanOrigin - 25)),
            .default = NA
        )
    ) %>%
    mutate(
        mutation_time = case_when(
            organism == "H1N1" ~ 0.004 * 13258,
            organism == "SARS-CoV-2" ~ 0.001 * 29903,
            organism == "Shigella" ~ 0.0000006 * 4825265,
            organism == "TB" ~ 0.00000001 * 4300000,
        )
    ) %>%
    mutate(
        res_lost = case_when(
            resolution == "Day" ~ 1 / 365.25,
            resolution == "Month" ~ 31 / 365.25,
            resolution == "Year" ~ 1,
        )
    ) %>%
    mutate(
        num_conflated_mutations = res_lost * mutation_time
    )

## map error time

err_vs_res %>%
    # filter(
    #     resolution != "Year"
    # ) %>%
    ggplot(
        aes(x = num_conflated_mutations, y = abs(OriginErr), fill = organism)
    ) +
    geom_point(pch = 21) #+
    # geom_smooth(
    #     se = FALSE
    # ) +
    #facet_wrap(~organism, scales = "free")
