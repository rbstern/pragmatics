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

sigma0 = 1
grid <- tibble()
for(sigma1 in sigmas)
{
  radius = sigma1*(1+2*eps-sigma0/sigma1-log(sigma1/sigma0))
  inter_left = mu0 - radius
  inter_right = mu0 + radius
  color = (mus > inter_left) & (mus < inter_right)
  grid %<>% rbind(tibble(mu = mus, sigma = sigma1, color = color))
}

grid %>% 
  ggplot(aes(x = mu, y = sigma, color = color)) +
  geom_point()
