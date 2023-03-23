## -v1 models: 
# original logistic function

rm(list = ls())
library(tidyverse)
library(readr)
library(reshape2)
library(rstan)
options(mc.cores=4)
library(posterior)
library(shinystan)

# titration data
df <- read.csv('data/raw/EV71_raw_well_observations.csv', header = TRUE, sep = ',')
serumIDs <- unique(df$serumID)
serums_subset <- serumIDs
df <- df %>% filter(serumID %in% serums_subset)
df <- df %>% mutate(sample = as.numeric(as.factor(serumID)))

df_long <- df %>% dplyr::select(serumID, dilutions, outcome, sample)
head(df_long)

# pick out all data variants
bars <- df_long %>% group_by(sample) %>% 
  summarise(barcode = paste0(outcome, collapse = "")) %>% ungroup() %>% 
  pull(barcode) %>% unique() # unique barcodes

bars_lookup <- df_long %>% group_by(sample) %>% 
  summarise(barcode = paste0(outcome, collapse = ""))

df_long <- df_long %>% left_join(bars_lookup)
head(df_long)
min(df_long$sample); max(df_long$sample)

df_each_code <- df_long %>% group_by(barcode) %>% 
  summarise(sampleID = first(sample)) %>% mutate(keep = TRUE)

df_bars_summ <- bars_lookup %>% group_by(barcode) %>% count() %>% left_join(df_long)
head(df_bars_summ)
df_panels <- df_bars_summ %>% group_by(barcode) %>% summarise(sample = first(sample))

ndilutions_per_barcode <- 50

dilution1 <- 10 ^ seq(log10(4), log10(1100), length.out = ndilutions_per_barcode) 
dilution_sim <- rep(dilution1, n_distinct(df_bars_summ$barcode))
sample_sim <- unlist(map(seq(1, n_distinct(df_bars_summ$barcode)), ~ rep(., ndilutions_per_barcode)))

N_sim <- length(sample_sim)

stan_data <- list(
  N_panels = nrow(df_panels),
  phi_panel_index = df_panels$sample, # use this to extract data and plot ppc
  N = nrow(df),
  nreplicates = rep(2, nrow(df)),
  survival = df$outcome,
  dilution = df$dilutions,
  nsample = n_distinct(df$sample),
  sample = df$sample,
  N_sim = N_sim,
  sample_sim = sample_sim,
  dilution_sim = dilution_sim
)

#--------------------------------------------------------------------------------------------
# fit model

# ref - optimizing: (http://mc-stan.org/rstan/reference/stanmodel-method-optimizing.html)

model = stan_model("src/stan/E71_without_c.stan")


is_converged <- FALSE
max_tries <- 10
tries <- 0
while(!is_converged & tries <= max_tries) {
  opt <- rstan::optimizing(model, data = stan_data, init = 'random', as_vector = FALSE)
  is_converged <- if_else(opt$return_code==0, TRUE, FALSE)
  tries <- tries + 1
}
initf <- function() { 
  list( a = opt$par$a, b = opt$par$b, phi = opt$par$phi ) 
}

fit <- rstan::sampling(model, data = stan_data, chains = 4, iter = 400, init = initf)
fit
launch_shinystan(fit)

#--------------------------------------------------------------------------------------------
# analyze the posterior:

# (1)
fit <- readRDS("hpc_OXF/fitAll_model_nAll_v1.rds") 

df1 <- rstan::extract(fit, c("a", "b", "phi")) %>% as.data.frame()
df_sum <- summarize_draws(df1)
sum(df_sum$rhat > 1.01)

fit_Summary <- data.frame(rstan::summary(fit)$summary)
medians <- rstan::summary(fit)$summary[ , "50%"]
fit_Summary$median <- medians
head(fit_Summary)

# posterior predictive checks
prob <- apply(rstan::extract(fit, "prob")[[1]], 2, median)
prob_low <- apply(rstan::extract(fit, "prob")[[1]], 2, function(x) quantile(x, 0.025))
prob_high <- apply(rstan::extract(fit, "prob")[[1]], 2, function(x) quantile(x, 0.975))
indices <- df_panels$sample

df2 <- df[df$sample %in% (df_panels$sample), ] %>% 
  arrange(factor(sample, levels = indices)) %>% 
  mutate(prob = outcome / 2,
         sample = unlist(map(seq(1, n_distinct(df_bars_summ$barcode)), ~ rep(., 8)))
  )

df_sim <- tibble(dilutions = dilution_sim, sample = sample_sim, prob = prob, 
                 prob_low = prob_low, prob_high = prob_high)
ggplot(data = df2, aes(x = dilutions, y = prob)) +
  geom_point(colour = "orange") +
  geom_line(colour = "orange") +
  geom_ribbon(data = df_sim, aes(ymin = prob_low, ymax = prob_high), alpha = 0.4) +
  scale_x_log10() +
  labs(x = "Dilution", y = "Probability of cell survival") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 8),
        strip.text = element_blank(),
        panel.grid = element_blank()) +
  facet_wrap(~ sample)

g

ggsave("figures/section_3.3/EVA71/E71_ppc_E71_model1_cauchy.pdf", g, width = 6, height = 6, dpi = 300)

# (2)

post_a <- rstan::extract(fit, "a")[[1]]
hist(post_a, main = "nsamples = 1000", xlab = 'a')
post_b <- rstan::extract(fit, "b")[[1]]
hist(post_b, main = "nsamples = 1000", xlab = 'b')

# extract phi values:
phiEst <- fit_Summary[grep("phi", rownames(fit_Summary)), 'mean']
phiEst2.5 <- fit_Summary[grep("phi", rownames(fit_Summary)), 'X2.5.']
phiEst97.5 <- fit_Summary[grep("phi", rownames(fit_Summary)), 'X97.5.']
range(phiEst); length(unique(phiEst))
hist(phiEst)

# plot phis:
ggplot() + 
  geom_histogram(aes(x = phiEst), alpha = 0.5, position = 'identity', bins = 40) +
  labs(x = 'phi', y = 'Count') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 8))


# combine phi values and titers:
phis_titers <- df

for(i in 1:nrow(phis_titers)){
  phis_titers[i, 'phi'] <- phiEst[phis_titers$sample[i]] 
  phis_titers[i, 'phi2.5'] <- phiEst2.5[phis_titers$sample[i]] 
  phis_titers[i, 'phi97.5'] <- phiEst97.5[phis_titers$sample[i]] 
}

head(phis_titers)

e71 <- phis_titers[, c('serumID','titer','phi','phi2.5','phi97.5')]
e71B <- e71[!duplicated(e71$serumID), ]

# simple plot of phis vs titer:
plot(log10(e71B$titer), e71B$phi, xlab = 'log10 titer', main = "Titer vs. phi", ylab = 'phi')

# with credible intervals:
h <- ggplot(data = e71B, aes(x = log10(titer), y = log(phi))) +
  geom_pointrange(aes(ymin = log(phi2.5), ymax = log(phi97.5)), color = "orange") +
  labs(title = "EV-A71") +
  scale_x_continuous(breaks = round(seq(min(log10(e71B$titer)), 
                                        max(log10(e71B$titer)), by = 0.5),1)) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9))

h
ggsave("figures/section_3.3/EVA71/E71_logphi_titers_AllSamples-v1.pdf", h, width = 6, height = 6, dpi = 300)

phis_titers %>%
  distinct(serumID, .keep_all = TRUE) %>%
  ggplot() + 
  geom_histogram(aes(x = phi), alpha = 0.5, position = 'identity', bins = 40) +
  labs(x = 'phi', y = '') +
  theme_bw()

#---------------------------------------------------------------------------------------------
# by age and age group:

phis_titers <- phis_titers %>% mutate(ageGp = case_when(Age > 60 ~ '>60 y',
                                                        Age > 40  & Age <= 60 ~ '41-60 y',
                                                        Age > 20  & Age <= 40 ~ '21-40 y',
                                                        Age > 10  & Age <= 20 ~ '11-20 y',
                                                        Age > 5  & Age <= 10 ~ '6-10 y',
                                                        Age >= 1  & Age <= 5 ~ '1-5 y',
                                                        Age < 1  ~ '<1 y'))

phis_titers <- phis_titers[!(phis_titers$Age <= 1), ]
str(phis_titers)
levels(as.factor(phis_titers$ageGp)); levels(as.factor(phis_titers$Year)) 

j <- phis_titers %>%
  distinct(serumID, .keep_all = TRUE) %>%
  ggplot(aes(x = factor(ageGp, 
                        levels = c("1-5 y","6-10 y","11-20 y","21-40 y","41-60 y",">60 y")), 
             y = log(phi))) + 
  geom_violin(color = "orange") +
  geom_jitter(color = "orange", position = position_jitter(0.2)) +
  labs(x = "\nAge groups", y = "log(phi)") +
  theme_bw() +
  theme(axis.title = element_text(size = 11),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
        strip.text.x = element_text(size = 12)) +
  facet_grid(~Year)

j
ggsave("figures/section_3.3/EVA71/E71_logphi_byAge_AllSamples-v1.pdf", j, width = 8, height = 6, dpi = 300)

k <- phis_titers %>%
  distinct(serumID, .keep_all = TRUE) %>%
  ggplot(aes(x = factor(ageGp, 
                        levels = c("1-5 y","6-10 y","11-20 y","21-40 y","41-60 y",">60 y")), 
             y = phi)) + 
  geom_violin(color = "orange") +
  geom_jitter(color = "orange", position = position_jitter(0.2)) +
  labs(x = "\nAge groups", y = "phi") +
  theme_bw() +
  theme(axis.title = element_text(size = 11),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
        strip.text.x = element_text(size = 12)) +
  facet_grid(~Year)

k
ggsave("figures/section_3.3/EVA71/E71_phi_byAge_AllSamples-v1.pdf", k, width = 8, height = 6, dpi = 300)


phis_titers %>%
  distinct(serumID, .keep_all = TRUE) %>%
  ggplot(aes(x = Age, y = phi)) +
  geom_point(color = "orange") +
  theme_bw() +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        strip.text.x = element_text(size = 12)) +
  facet_grid(~Year)


l <- phis_titers %>%
  distinct(serumID, .keep_all = TRUE) %>%
  ggplot(aes(x = as.factor(Year), y = log(phi), group = as.factor(Year))) +
  geom_violin(color = "orange") + 
  geom_jitter(color = "orange", position = position_jitter(0.2)) +
  labs(x = 'Year') +
  theme_bw() +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 9))

l
ggsave("figures/section_3.3/EVA71/E71_logphi_byYear_AllSamples-v1.pdf", l, width = 8, height = 6, dpi = 300)

#---------------------------------------------------------------------------------------------
# plot sigmoid curve:
sigmoid <- function(a, b, phi){
  1 / (1 + exp(-(a + b*log(phi))))
}

plot(phiEst, sigmoid(fit_Summary["a", "mean"], fit_Summary["b", "mean"], phiEst), ylab='')

#---------------------------------------------------------------------------------------------
