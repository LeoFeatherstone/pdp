library(plotly)

r <- seq(from = 1e-8, to = 1e-2, length.out = 100)
D <- seq(4 / 365.25, 2, length.out = 100)
G <- seq(1e04, 1e09, length.out = 100)

rG <- seq(min(r %*% t(G)), max(r %*% t(G)), length.out = 100)

meff <- rG %*% t(D)
rownames(meff) <- as.character(rG)
colnames(meff) <- as.character(D)

plot_ly(z = meff) %>% add_surface()
