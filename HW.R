library(combinat)
library(DirichletReg)
library(tidyverse)
library(magrittr)

############################################################
## Global parameters                                      ##
m0 = 20                                                   ##
p0 = 1/3                                                  ##
eps = 0.1                                                 ##
B = 200                                                   ##
thetas = seq(0.0000001, 0.99999999, length.out = B)       ##
                                                          ##
simplex = tibble(theta1 = 0, theta3 = thetas) %>%         ##
  bind_rows(tibble(theta1 = thetas, theta3 = 0)) %>%      ##
  bind_rows(tibble(theta1 = thetas, theta3 = 1 - thetas)) ##
                                                          ##
color_grid = tibble(theta1 = rep(thetas, each = B),       ##
                    theta2 = rep(thetas, B),              ##
                    color = rep(FALSE, B^2)) %>%          ##
  mutate(theta3 = 1 - theta1 - theta2) %>%                ##
  filter(theta3 > 0, theta3 < 1)                          ##
############################################################

# Generate HPD's

## Useful functions

posterior = function(data, prior_par = c(1, 1, 1)) data + prior_par

log_posterior = function(t1, t2, t3, posterior_par)
{
  ddirichlet(cbind(t1, t2, t3), posterior_par, log = TRUE)
}

# The HPD is a region with parameters
# such that the posterior is greater than v.
# This function determines v such that
# the HPD has posterior probability of 'probability'.
hpd_post_cut = function(posterior_par, probability = 0.95, B = 10^5)
{
  quantile = posterior_par %>%
    rdirichlet(n = B, alpha = .) %>%
    ddirichlet(posterior_par, log = TRUE) %>%
    quantile(probs = 1 - probability) %>%
    as.numeric()
  return(quantile)
}

## Obtain data

### Brentani data
data = "./data/brentani.csv" %>%
  read.csv() %>%
  as.tibble()

### Sample not from H-W
set.seed(1)
sample_new = sample(1:3, 1000, prob = c(0.2, 0.3, 0.5), replace = TRUE)
sample_new <- table(sample_new)
names(sample_new) = colnames(data)

### Sample from H-W
p <- 0.8
sample_new2 = sample(1:3, 1000, prob = c(p^2, 2*p*(1-p), (1-p)^2), replace = TRUE)
sample_new2 <- table(sample_new2)
names(sample_new2) = colnames(data)

### Complete data
data %<>% bind_rows(sample_new, sample_new2)

## Generate HPD's
hpd_labels = tibble(idx = 1:nrow(data), 
                    x = rep(NA, nrow(data)), 
                    y = rep(NA, nrow(data)))
hpd_chull_grid = NULL
for(ii in 1:nrow(data))
{
  counts = data[ii,]
  post = counts %>% posterior() %>% as.numeric()
  th = hpd_post_cut(post, probability = 0.8)
  this_hpd_grid = color_grid %>%
    mutate(color = (log_posterior(theta1, theta2, theta3, post) > th)) %>%
    filter(color) %>%
    select(theta1, theta3)
  this_chull = chull(this_hpd_grid)
  this_hpd_chull = this_hpd_grid[this_chull,] %>%
    mutate(idx = ii)
  hpd_chull_grid %<>% bind_rows(this_hpd_chull)
  hpd_labels[ii,] = c(ii, (post/sum(post))[1], (post/sum(post))[3])
}
hpd_chull_grid %<>% mutate(idx = as.factor(idx))

# Helpful functions

BP_diss = function(t1, t2, t3, p, m = 1)
{
  numerator = m*((t1 - p^2)^2 + (t2 - 2*p*(1-p))^2 + (t3 - (1-p)^2)^2)
  denominator = p^2*(1 - p^2) + 2*p*(1 - p)*(1 - 2*p*(1 - p)) + (1 - p)^2*(1 - (1 - p)^2)
  numerator/denominator
}

# Kullback-Leibler divergence between 
# Multinomial(t1, t2, t3) and
# Multinomial(p^2, 2p(1-p), (1-p)^2)
KL_diss = function(t1, t2, t3, p0, m = 1)
{
  this_KL = rep(NA, length(t1))
  # Multinomial probabilities for each possible value under Hardy-Weinberg and p.
  p <- c(p0^2, 2*p0*(1-p0), (1-p0)^2)
  log_mult_p = log(xsimplex(3, m, dmnom, prob = p))
  for(ii in 1:length(t1))
  {
    # Multinomial probabilities for each possible value using (t1, t2, t3).
    mult_theta = xsimplex(3, m, dmnom, prob = c(t1[ii], t2[ii], t3[ii]))
    log_mult_theta = log(mult_theta)
    # Calculate KL
    this_KL[ii] = sum(mult_theta*(log_mult_theta - log_mult_p))
  }
  this_KL
}

C_diss = function(t1, t2, t3, p0, m = 1)
{
  this_C = rep(NA, length(t1))
  # Multinomial probabilities for each possible value under Hardy-Weinberg and p.
  p <- c(p0^2, 2*p0*(1-p0), (1-p0)^2)
  mult_p = xsimplex(3, m, dmnom, prob = p)
  for(ii in 1:length(t1))
  {
    mult_theta = xsimplex(3, m, dmnom, prob = c(t1[ii], t2[ii], t3[ii]))
    this_C[ii] = 0.25*sum(abs(mult_p - mult_theta))
  }
  this_C
}

generate_grid = function(dissimilarity, p0, eps, m = 1)
{
  color_grid %>%
    mutate(diss = dissimilarity(theta1, theta2, theta3, p0, m)) %>%
    mutate(color = (diss < eps))
}

generate_join_grid = function(dissimilarity, eps, m = 1)
{
  join_grid = generate_grid(dissimilarity, thetas[1], eps, m)
  for(theta in thetas)
  {
    print(theta)
    new_grid = generate_grid(dissimilarity, theta, eps, m)
    join_grid$color = (join_grid$color | new_grid$color)
  }
  join_grid
}

# Experiments 
file_path = function(name) paste("./data/HW_", 
                                 name, "_", m0, "_",  eps, "_", B, 
                                 ".rds", sep = "")

join_grid = generate_join_grid(BP_diss, eps, m0)
write_rds(join_grid, file_path("BP"))

join_grid = generate_join_grid(KL_diss, eps, m0)
write_rds(join_grid, file_path("KL"))

join_grid = generate_join_grid(C_diss, eps, m0)
write_rds(join_grid, file_path("C"))

# Plots

## Helpful plot functions

plot_grid = function(grid, col = "dodgerblue4", alpha = 0.2)
{
  grid %<>% filter(color) %>%
    group_by(theta1)
  poly_min = grid %>% summarise(theta3 = min(theta3))
  poly_max = grid %>% 
    summarise(theta3 = max(theta3)) %>%
    arrange(desc(theta1))
  poly_grid = bind_rows(poly_min, poly_max)
  poly_grid %>%
    ggplot(aes(x = theta1, y = theta3)) +
    geom_polygon(color = col, fill = col, alpha = alpha) +
    geom_path(aes(x = theta1, y = theta3), data = simplex) +
    stat_function(fun = function(x) (1-sqrt(x))^2, n = 101, color = "red") +
    xlab(expression(theta[1])) +
    ylab(expression(theta[3]))
}

add_hpd_grid = function(plot, hpd_chull_grid, hpd_labels, col = "green4", alpha = 0.2)
{
  all_idx = unique(hpd_chull_grid$idx)
  new_plot = plot
  for(this_idx in all_idx)
  {
    this_chull = hpd_chull_grid %>% filter(idx == this_idx)
    new_plot = new_plot +
      geom_polygon(aes(x = theta1, y = theta3), data = this_chull,
                   color = col, fill = col, alpha = alpha)
  }
  new_plot +
  geom_text(aes(x = x, y = y, label = idx), data = hpd_labels)
}

figure_path = function(name, ext, extra="") paste("./figures/HW_", 
                                                  name, "_", m0, "_", eps, "_", B, 
                                                  extra, ext, sep="")

## Actual plots

### BP_diss
simple_grid = generate_grid(BP_diss, p0, eps, m0)
plot_grid(simple_grid)
ggsave("./figures/hw_bp_0.01_200_0.33.pdf")
ggsave("./figures/hw_bp_0.01_200_0.33.png")

join_grid = read_rds(file_path("BP"))
plot_grid(join_grid)
ggsave("./figures/hw_bp_0.01_200.pdf")
ggsave("./figures/hw_bp_0.01_200.png")

plot_grid(join_grid) %>% add_hpd_grid(hpd_chull_grid, hpd_labels)
ggsave("./figures/hw_bp_0.01_200_hpd.pdf")
ggsave("./figures/hw_bp_0.01_200_hpd.png")

### KL_diss
simple_grid = generate_grid(KL_diss, p0, eps)
plot_grid(simple_grid)
ggsave("./figures/hw_kl_0.01_200_0.33.pdf")
ggsave("./figures/hw_kl_0.01_200_0.33.png")

join_grid = read_rds(file_path("KL"))
plot_grid(join_grid)
ggsave("./figures/hw_kl_0.01_200.pdf")
ggsave("./figures/hw_kl_0.01_200.png")

plot_grid(join_grid) %>% add_hpd_grid(hpd_chull_grid, hpd_labels)
ggsave("./figures/hw_kl_0.01_200_hpd.pdf")
ggsave("./figures/hw_kl_0.01_200_hpd.png")

### C_diss
simple_grid = generate_grid(C_diss, p0, eps, m0)
plot_grid(simple_grid)
ggsave("./figures/hw_c_0.01_200_0.33.pdf")
ggsave("./figures/hw_c_0.01_200_0.33.png")

join_grid = read_rds(file_path("C"))
plot_grid(join_grid)
ggsave(figure_path("C", ".pdf"))
ggsave(figure_path("C", ".png"))

plot_grid(join_grid) %>% add_hpd_grid(hpd_chull_grid, hpd_labels)
ggsave(figure_path("C", ".pdf", extra = "_hpd"))
ggsave(figure_path("C", ".png", extra = "_hpd"))
