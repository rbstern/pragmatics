library(magrittr)
library(tidyverse)

#################################################
## X ~ N(mu, sigma^2), where sigma is unknown. ##
## Pragmatic hypothesis for H0: mu = 0 using   ##
## KL divergence, epsilon = 0.1 and M^2 = 2    ##
#################################################
eps = 0.1
n = 100
mu0 = 0
mus = seq(from = -1, to = 1, length.out = n)
sigmas = seq(from = 1/n, to = 2, length.out = n)

sigma_grid = list()
for(sigma0 in sigmas)
{
  new_grid <- tibble()
  for(sigma1 in sigmas)
  {
    radius = sigma1*(1+2*eps-sigma0/sigma1-log(sigma1/sigma0))
    inter_left = mu0 - radius
    inter_right = mu0 + radius
    color = (mus > inter_left) & (mus < inter_right)
    new_grid %<>% rbind(tibble(mu = mus, sigma = sigma1, color = color))
  }
  sigma_grid[[as.character(sigma0)]] = new_grid
}

#write_rds(sigma_grid, "./data/KL_norm_0.1_100_0.rds")
#sigma_grid = sigma_grid = read_rds("./data/KL_norm_0.1_100_0.rds")

small_grid = sigma_grid[[round(n/2)]]
join_grid = tibble(mu = small_grid$mu, sigma = small_grid$sigma)
max_col = sigma_grid[[1]]$color
for(sigma0 in as.character(sigmas))
{
  max_col = (max_col | sigma_grid[[sigma0]]$color)
}
join_grid$color = max_col

join_grid %>% 
  filter(color == TRUE) %>%
  ggplot(aes(x = mu, y = sigma)) +
  geom_point()
