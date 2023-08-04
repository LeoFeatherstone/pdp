#############################################
###### Looking at simulated sequence ########
#############################################

library(tidyverse)

load("sim_study/alignment_data.RData")

data <- data %>%
    rownames_to_column(var = "pathogen") %>%
    mutate(pathogen = gsub(pattern = "_.+", replacement = "", pathogen))

# Site patterns plot
pdf("sim_data.pdf", useDingbats = FALSE)
    ggplot(data) +
    geom_point(
        aes(
            x = (num_samples),
            y = (site_patterns),
            col = pathogen
        )) +
    # scale_x_continuous(trans = "log10") +
    # scale_y_continuous(trans = "log10") +
    facet_wrap(~pathogen,scales = "free") +
    # xlim(0, Inf) + ylim(0, Inf) +
    geom_abline(intercept = 0, slope = 1) #+
    #theme(aspect.ratio = 1)
dev.off()

# sampling hist
date_df <- data.frame()
for (i in seq_along(dates)) {
    df <- cbind(dates[[i]], rep(names(dates[i]), times = length(dates[[i]])))
    colnames(df) <- c("date", "dataset")
    date_df <- rbind(date_df, df)
}

date_df <- date_df %>%
    mutate(dataset = gsub(pattern = "[.]fasta", replacement = "", dataset)) %>%
    separate_wider_delim(dataset, "_", names = c("Pathogen", "Replicate")) %>%
    mutate(date = as.Date(as.numeric(date), origin = "1970-01-01"))

pdf("sim_samp_traj.pdf", useDingbats = FALSE)
    ggplot(date_df) +
    geom_freqpoly(
        aes(x = date, color = Pathogen, group = Replicate),
        alpha = .25,
        binwidth = 7
    ) +
    ylab("Weekly Samples") +
    facet_wrap(~Pathogen, scales = "free")
dev.off()
