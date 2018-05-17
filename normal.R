library(magrittr)
library(tidyverse)

#################################################
## X ~ N(mu, sigma^2), where sigma is unknown. ##
## Pragmatic hypothesis for H0: mu = 0 using   ##
##  epsilon = 0.1 and M^2 = 2, KL and CD.      ##
#################################################

###################################################
## Global parameters                             ##
mu0 = 0 # null hypothesis mean                   ##
sigma1 = 1 # sigma for simple hypothesis         ##
B = 301 # grid size                              ##
eps = 0.1                                        ##
###################################################

###################################################
## Useful variables                              ##
mus = seq(from = -1, to = 1, length.out = B)     ##
sigmas = seq(from = 1/B, to = 2, length.out = B) ##
###################################################

# Useful functions

## Kullback-Leibler grid for null with sigma = sigma0
normal_kl_grid = function(sigma0)
{
  new_grid <- tibble(mu = rep(mus, B), 
                     sigma = rep(sigmas, each = B), 
                     color = rep(NA, B*B))
  for(ii in 1:B)
  {
    color = rep(FALSE, B)
    radius = sigmas[ii]*(1+2*eps-sigma0/sigmas[ii]-log(sigmas[ii]/sigma0))
    if(radius > 0)
    {
      inter_left = mu0 - sqrt(radius)
      inter_right = mu0 + sqrt(radius)
      color = (mus > inter_left) & (mus < inter_right)
    }
    new_grid$color[(B*(ii-1)+1):(B*ii)] = color
  }
  return(new_grid)
}

## Classification distance grid for null with sigma = sigma0
normal_cd_grid = function(sigma0)
{
  new_grid <- tibble(mu = rep(mus, B),
                     sigma = rep(sigmas, each = B), 
                     color = rep(NA, B*B))
  for(ii in 1:B)
  {
    for(jj in 1:B)
    {
      l_1 = function(x) abs(dnorm(x, mus[jj], sigmas[ii]) - dnorm(x, mu0, sigma0))
      cd = 0.25*integrate(l_1, lower = -Inf, upper = Inf)$value
      new_grid$color[B*(ii-1)+jj] = (cd <= eps)
    }
  }
  return(new_grid)
}

## grid with unknown sigma
generate_join_grid = function(grid_method)
{
  join_grid = grid_method(sigmas[1])
  for(sigma0 in sigmas)
  {
    print(paste(sigma0, "/", max(sigmas)))
    new_grid = grid_method(sigma0)
    join_grid$color = (join_grid$color | new_grid$color)
  }
  join_grid
}

# Experiments

file_path = function(method) paste("./data/normal_", method, "_", 
                                   mu0, "_", round(eps, 2), "_", B, 
                                   ".rds", sep = "")

join_kl_grid = generate_join_grid(normal_kl_grid)
write_rds(join_kl_grid, file_path("KL"))

join_cd_grid = generate_join_grid(normal_cd_grid)
write_rds(join_cd_grid, file_path("C"))

# Plots

## Useful plot function

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

figure_path = function(name, ext, extra="") paste("./figures/normal_",
                                                  name, "_", mu0, "_", round(eps,2), "_", B,
                                                  extra, ext, sep = "")

simple_extra = paste("_", round(sigma1, 2), sep = "")

## KL plots

unit_sigma_kl_grid = normal_kl_grid(sigma1)
plot_norm_grid(unit_sigma_kl_grid)
ggsave(figure_path("KL", ".pdf", extra = simple_extra))
ggsave(figure_path("KL", ".png", extra = simple_extra))

join_kl_grid = read_rds(file_path("KL"))
plot_norm_grid(join_kl_grid)
ggsave(figure_path("KL", ".pdf"))
ggsave(figure_path("KL", ".png"))

## CD plots

unit_sigma_cd_grid = normal_cd_grid(sigma1)
plot_norm_grid(unit_sigma_cd_grid)
ggsave(figure_path("C", ".pdf", extra = simple_extra))
ggsave(figure_path("C", ".png", extra = simple_extra))

join_cd_grid = read_rds(file_path("C"))
plot_norm_grid(join_cd_grid)
ggsave(figure_path("C", ".pdf"))
ggsave(figure_path("C", ".png"))
