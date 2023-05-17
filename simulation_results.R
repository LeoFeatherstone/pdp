#########################################
### Wrangle and compare poterior data ###
#########################################

library(tidyverse)
library(coda)
library(GGally)

## function applies burnin
# df is data frame
# pc is percent burnin
apply_burnin <- function(df, pc) {
    return(df[-c(0:ceiling((pc / 100) * dim(df)[1])), ])
}

##### Empirical Data #####

sim_log <- dir(path = "./sim_study", pattern = ".+[.]log")

sim_data <- lapply(
    sim_log,
    function(x) {
        print(paste0("Read ", x))
        return(read.table(paste0("./sim_study/", x), header = TRUE))
    }
)
names(sim_data) <- gsub(
    sim_log,
    pattern = "[.]log",
    replacement = ""
)
# save for speed later
#save(sim_data, file = "raw_sim_posteriors.RData")
load("raw_sim_posteriors.RData")
# burnin
sim_data <- lapply(
    sim_data,
    function(x) apply_burnin(x, 25)
)

# patch names with clock until rerun. # UPTO - remove _SC_ from tb names
for (i in seq_along(names(sim_data))){
    if (grepl(pattern = "tb_", names(sim_data)[i])) {
        names(sim_data)[i] <- gsub(
            names(sim_data)[i],
            pattern = "_SC_BD",
            replacement = "_BD"
        )
    }
}

# format variables names
rename_cols <- function(vec) {
    vec <- gsub(
        vec,
        pattern = "reproductiveNumber_BDSKY_Serial[.]",
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

sim_data <- lapply(
    sim_data,
    function(x) {
        colnames(x) <- rename_cols(colnames(x))
        return(x)
    }
)




sim_ess <- lapply(
    sim_data,
    function(x) {
        effectiveSize(
            as.mcmc(x[, -grep(pattern = "Sample", colnames(x))])
        )
    }
)

# Find which treatments had ess < 200
lapply(
    seq_along(sim_ess),
    function(x) {
        if (any(sim_ess[[x]] < 200)) {
            print(names(sim_ess)[x])
            print(sim_ess[x])

        }
    }
)

# get parms of interest
sim_data <- lapply(
    sim_data,
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
sim_data <-
    bind_rows(sim_data, .id = "id") %>%
    separate_wider_delim(
        id,
        delim = "_",
        names = c("organism", "treePrior", "resolution", "replicate")
    )

# sim plots

#tmp <- 
sim_data %>%
    filter(organism == "h1n1") %>%
    select(resolution, replicate, reproductiveNumber_BDSKY_Serial) %>%
    group_by(replicate, resolution) %>%
    summarise(meanR0 = mean(reproductiveNumber_BDSKY_Serial)) %>%
    pivot_wider(
        names_from = resolution,
        values_from = c(meanR0)) %>%
    ggparcoord(
        columns = 2:3,
        scale = "globalminmax"
    )

sim_data %>%
    filter(
        organism == "h1n1" |
        organism == "sars-cov-2" |
        organism == "shigella"
    ) %>%
    select(resolution, replicate, organism, reproductiveNumber_BDSKY_Serial) %>%
    group_by(replicate, resolution, organism) %>%
    summarise(meanR0 = mean(reproductiveNumber_BDSKY_Serial)) %>%
    pivot_wider(
        names_from = resolution,
        values_from = c(meanR0)) %>%
    ggparcoord(
        columns = 3:4,
        scale = "globalminmax",
        groupColumn = "organism"
    )


sim_data %>%
    filter(organism == "shigella" | organism == "tb") %>%
    select(resolution, replicate, organism, reproductiveNumber_BDSKY_Serial.1) %>%
    group_by(replicate, resolution, organism) %>%
    summarise(meanR0 = mean(reproductiveNumber_BDSKY_Serial.1)) %>%
    pivot_wider(
        names_from = resolution,
        values_from = c(meanR0)) %>%
    ggparcoord(
        columns = 3:5,
        scale = "globalminmax",
        groupColumn = "organism"
    )

sim_data %>%
    filter(organism == "shigella" | organism == "tb") %>%
    select(resolution, replicate, organism, reproductiveNumber_BDSKY_Serial.2) %>%
    group_by(replicate, resolution, organism) %>%
    summarise(logMeanR0 = log10(mean(reproductiveNumber_BDSKY_Serial.2))) %>%
    pivot_wider(
        names_from = resolution,
        values_from = c(logMeanR0)) %>%
    ggparcoord(
        columns = 3:5,
        scale = "globalminmax",
        groupColumn = "organism"
    )


# clock rate
sim_data %>%
    #filter(organism != "h1n1") %>%
    select(resolution, replicate, organism, clockRate) %>%
    group_by(replicate, resolution, organism) %>%
    summarise(mean_clock_rate = ((mean(clockRate)))) %>%
    pivot_wider(
        names_from = resolution,
        values_from = c(mean_clock_rate)) %>%
    ggparcoord(
        columns = 3:5,
        scale = "globalminmax",
        alphaLines = 0.4,
        groupColumn = "organism"
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
