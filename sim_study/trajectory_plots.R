############################################################
####### Figures for likelihood and rate trajectories #######
############################################################

library(tidyverse)
library(coda)
library(tracerer)
library(gganimate)

load("./sim_study/posteriors.RData")

## function applies burnin
# df is data frame
# pc is percent burnin
apply_burnin <- function(df, pc) {
    return(df[-c(0:ceiling((pc / 100) * dim(df)[1])), ])
}

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

# bind data and add phylodynamic likelihood
data_df <-
    bind_rows(data, .id = "id") %>%
    separate_wider_delim(
        id,
        delim = "_",
        names = c("organism", "clock", "treePrior", "resolution", "replicate")
    ) %>%
    rename("phylogeneticLikelihood" = likelihood) %>%
    mutate(
        #phylodynamicLikelihood = posterior - phylogeneticLikelihood - prior
        phylodynamicLikelihood = prior
    )

# work with randomly chosen subset for now
df <- data_df %>%
    filter(
        organism == "h1n1",
        resolution == "Day" | resolution == "Month",
        replicate == "1"
    )

# Likelihood phase
# df %>%
#     ggplot(aes(x = phylogeneticLikelihood, y = phylodynamicLikelihood)) +
#     geom_point() +
#     facet_wrap(~resolution, scales = "free") +
#     theme_minimal() +
#     theme(aspect.ratio = 0.9) +
#     transition_states(
#         Sample,
#         transition_length = 1,
#         state_length = 1
#   ) 

# anim_save("likelihood_traj.gif")

plot_a <- df %>%
    ggplot(aes(x = phylogeneticLikelihood, y = phylodynamicLikelihood)) +
    geom_point(pch = 16, alpha = 0.2) +
    facet_wrap(~resolution) +
    theme_minimal() +
    theme(aspect.ratio = 0.9)

ggsave("likelihood_plane.pdf", dpi = 300)

# rates
plot_b <- df %>%
    ggplot(aes(x = R0, y = clockRate)) +
    geom_point(pch = 16, alpha = 0.2, fill = "red") +
    geom_hline(yintercept = 0.004) + 
    geom_vline(xintercept = 1.3) +
    facet_wrap(~resolution) +
    theme_minimal() 


ggsave("rates_plane.pdf", dpi = 300)

cowplot::plot_grid(plot_a, plot_b, nrow = 2)
ggsave("plane.pdf", dpi = 300)
