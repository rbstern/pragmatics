library(tidyverse)
library(magrittr)

#################################################
## X ~ N(mu, sigma^2), where sigma is unknown. ##
## Pragmatic hypothesis for H0: mu = 0 using   ##
##  epsilon = 0.1 and M^2 = 2, KL and CD.      ##
#################################################

#######################################################
## Global parameters                                 ##
mu0 = 0 # null hypothesis mean                       ##
sigma0 = 1                                           ##
B = 301 # grid size                                  ##
eps = 0.1                                            ##
simple_hypothesis = tibble(mu = mu0, sigma = sigma0) ##
#######################################################

###################################################
## Useful variables                              ##
mus = seq(from = -1.2, to = 1.2, length.out = B) ##
sigmas = seq(from = 1/B, to = 2, length.out = B) ##
###################################################

# Useful functions

## Kullback-Leibler grid for null
## mu0 = parameters$mu and
## sigma0 = parameters$sigma
normal_KL_grid = function(parameters)
{
  mu0 = parameters$mu
  sigma0 = parameters$sigma
  new_grid = tibble(mu = rep(mus, B), 
                    sigma = rep(sigmas, each = B), 
                    color = rep(NA, B*B))
  for(ii in 1:B)
  {
    color = rep(FALSE, B)
    radius = sigmas[ii]^2*(1 + 2*eps - 2*log(sigmas[ii]/sigma0)) -sigma0^2
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

## Classification distance grid for null 
## mu0 = parameters$mu and
## sigma0 = parameters$sigma
normal_C_grid = function(parameters)
{
  mu0 = parameters$mu
  sigma0 = parameters$sigma
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

## grid for composite hypothesis
generate_join_grid = function(grid_method, null_hypothesis)
{
  join_grid = grid_method(null_hypothesis[1,])
  for(ii in 2:nrow(null_hypothesis))
  {
    print(paste(ii, "/", nrow(null_hypothesis), sep = ""))
    new_grid = grid_method(null_hypothesis[ii,])
    join_grid$color = (join_grid$color | new_grid$color)
  }
  join_grid
}

## grid plot function
plot_norm_grid <- function(grid, null_hypothesis)
{
  grid %>% 
    ggplot(aes(x = mu, y = sigma)) +
    theme_minimal() +
    geom_tile(aes(alpha = 0.75, 
                  fill = as.numeric(color)), 
              show.legend = FALSE) +
    scale_fill_gradient(low = "white",
                        high = "steelblue") +
    geom_line(aes(x = mu, y = sigma,  
                   colour = "red"),
               show.legend = FALSE, 
               data = null_hypothesis) +
    xlab(expression(mu)) + 
    ylab(expression(sigma))
}

#################################
## Experiments for H_0: mu = 0 ##
#################################

# Parameter grid for H_0
null_hypothesis = tibble(mu = rep(mu0, B),
                         sigma = sigmas)

file_path = function(method) paste("./data/normal_", method, "_", 
                                   mu0, "_", round(eps, 2), "_", B, 
                                   ".rds", sep = "")

join_KL_grid = generate_join_grid(normal_KL_grid, null_hypothesis)
write_rds(join_KL_grid, file_path("KL"))

join_C_grid = generate_join_grid(normal_C_grid, null_hypothesis)
write_rds(join_C_grid, file_path("C"))

# Plots

figure_path = function(name, ext, extra="") paste("./figures/normal_",
                                                  name, "_", mu0, "_", round(eps,2), "_", B,
                                                  extra, ext, sep = "")

simple_extra = paste("_", round(sigma0, 2), sep = "")

## KL plots

unit_sigma_KL_grid = normal_KL_grid(simple_hypothesis)
plot_norm_grid(unit_sigma_KL_grid, null_hypothesis)
ggsave(figure_path("KL", ".pdf", extra = simple_extra))
ggsave(figure_path("KL", ".png", extra = simple_extra))

join_KL_grid = read_rds(file_path("KL"))
plot_norm_grid(join_KL_grid, null_hypothesis)
ggsave(figure_path("KL", ".pdf"))
ggsave(figure_path("KL", ".png"))

## C plots
unit_sigma_C_grid = normal_C_grid(simple_hypothesis)
plot_norm_grid(unit_sigma_C_grid, null_hypothesis)
ggsave(figure_path("C", ".pdf", extra = simple_extra))
ggsave(figure_path("C", ".png", extra = simple_extra))

join_C_grid = read_rds(file_path("C"))
plot_norm_grid(join_C_grid, null_hypothesis)
ggsave(figure_path("C", ".pdf"))
ggsave(figure_path("C", ".png"))

#######################################
## Experiments for H_0: sigma/mu = 1 ##
#######################################

###################################################
## Update for positive mu                        ##
mus = seq(from = 1/B, to = 3, length.out = B)    ##
###################################################

null_hypothesis = tibble(mu = sigmas, sigma = sigmas)

file_path = function(method) paste("./data/normal_vc_", method, "_", 
                                   mu0, "_", round(eps, 2), "_", B, 
                                   ".rds", sep = "")

join_KL_grid = generate_join_grid(normal_KL_grid, null_hypothesis)
write_rds(join_KL_grid, file_path("KL"))

join_C_grid = generate_join_grid(normal_C_grid, null_hypothesis)
write_rds(join_C_grid, file_path("C"))

# Plots
figure_path = function(name, ext, extra="") paste("./figures/normal_vc_",
                                                  name, "_", mu0, "_", round(eps,2), "_", B,
                                                  extra, ext, sep = "")

join_KL_grid = read_rds(file_path("KL"))
plot_norm_grid(join_KL_grid, null_hypothesis)
ggsave(figure_path("KL", ".pdf"))
ggsave(figure_path("KL", ".png"))

join_C_grid = read_rds(file_path("C"))
plot_norm_grid(join_C_grid, null_hypothesis)
ggsave(figure_path("C", ".pdf"))
ggsave(figure_path("C", ".png"))
