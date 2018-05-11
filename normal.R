library(magrittr)
library(tidyverse)

#################################################
## X ~ N(mu, sigma^2), where sigma is unknown. ##
## Pragmatic hypothesis for H0: mu = 0 using   ##
##  epsilon = 0.1 and M^2 = 2, KL and CD.      ##
#################################################
eps = 0.1
n = 300
mu0 = 0
mus = seq(from = -1, to = 1, length.out = n)
sigmas = seq(from = 1/n, to = 2, length.out = n)

# KL grids

normal_kl_grid = function(sigma0)
{
  new_grid <- tibble(mu = rep(mus, n), 
                     sigma = rep(sigmas, each = n), 
                     color = rep(NA, n*n))
  for(ii in 1:length(sigmas))
  {
    color = rep(FALSE, n)
    radius = sigmas[ii]*(1+2*eps-sigma0/sigmas[ii]-log(sigmas[ii]/sigma0))
    if(radius > 0)
    {
      inter_left = mu0 - sqrt(radius)
      inter_right = mu0 + sqrt(radius)
      color = (mus > inter_left) & (mus < inter_right)
    }
    new_grid$color[(n*(ii-1)+1):(n*ii)] = color
  }
  return(new_grid)
}

unit_sigma_kl_grid = normal_kl_grid(1)

join_grid = normal_kl_grid(sigmas[1])
for(sigma0 in sigmas)
{
  new_grid = normal_kl_grid(sigma0)
  join_grid$color = (join_grid$color | new_grid$color)
}

#Save data
#write_rds(join_grid, "./data/KL_norm_0.1_300_0.rds")
#join_kl_grid = read_rds("./data/KL_norm_0.1_100_0.rds")

# Classification distance grids

normal_cd_grid = function(sigma0)
{
  new_grid <- tibble(mu = rep(mus, n),
                     sigma = rep(sigmas, each = n), 
                     color = rep(NA, n*n))
  for(ii in 1:length(sigmas))
  {
    for(jj in 1:length(mus))
    {
      l_1 = function(x) abs(dnorm(x, mus[jj], sigmas[ii])-dnorm(x, mu0, sigma0))
      cd = 0.25*integrate(l_1, lower = -Inf, upper = Inf)$value
      new_grid$color[n*(ii-1)+jj] = (cd <= eps)
    }
  }
  return(new_grid)
}

unit_sigma_cd_grid = normal_cd_grid(1)

join_grid = normal_cd_grid(sigmas[1])
for(sigma0 in sigmas)
{
  print(sigma0)
  new_grid = normal_cd_grid(sigma0)
  join_grid$color = (join_grid$color | new_grid$color)
}

#Save data
#write_rds(join_grid, "./data/CD_norm_0.1_300_0.rds")
#join_cd_grid = read_rds("./data/CD_norm_0.1_300_0.rds")

# Plot function

plot_norm_grid <- function(grid)
{
  grid %>% 
    ggplot(aes(x = mu, y = sigma)) +
    theme_minimal() +
    geom_tile(aes(alpha = 0.75, 
                  fill = as.numeric(color)), 
              show.legend = FALSE) +
    scale_fill_gradient(low = "white",
                        high = "steelblue") +
    geom_segment(aes(x = 0, y = 0, 
                     xend = 0, yend = 2, 
                     colour = "red"),
                 show.legend = FALSE) +
    xlab(expression(mu)) + 
    ylab(expression(sigma^2))
}

#KL plots

plot_norm_grid(unit_sigma_kl_grid)
#ggsave("./figures/norm_kl_0.1_300_0_1.pdf")
#ggsave("./figures/norm_kl_0.1_300_0_1.png")

plot_norm_grid(join_kl_grid)
#ggsave("./figures/norm_kl_0.1_300_0.pdf")
#ggsave("./figures/norm_kl_0.1_300_0.png")

#CD plots

plot_norm_grid(unit_sigma_cd_grid)
#ggsave("./figures/norm_cd_0.1_300_0_1.pdf")
#ggsave("./figures/norm_cd_0.1_300_0_1.png")

plot_norm_grid(join_cd_grid)
#ggsave("./figures/norm_cd_0.1_300_0.pdf")
#ggsave("./figures/norm_cd_0.1_300_0.png")
