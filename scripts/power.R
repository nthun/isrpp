# Power simulation for the confirmatory hypotheses of the ISRPP study
#   H1: Publication pressure (respondent level) predicts QRP engagement
#   H2: Country-level research capacity ("herd") predicts QRP engagement,
#       net of GDP
# Both hypotheses are tested on two outcomes (lifetime and recent engagement).
# A replicate counts as a success only if both outcomes are significant in the
# same direction.

# lme4 is loaded before the tidyverse so tidyr's nest()/expand_grid() are not masked
library(lme4)
library(tidyverse)
library(furrr)
library(patchwork)

theme_set(theme_minimal())

# ---- Settings shared by all simulations --------------------------------------

n_sims <- 500          # replicates per condition
alpha <- .05

# Survey design
min_per_country <- 5   # minimum respondents for a country to be included
n_items <- 40          # number of QRP items
n_blocks <- 4          # items are presented in blocks; people can stop after any block
p_blocks <- c(.15, .15, .15, .55)  # probability of finishing 1, 2, 3, or all 4 blocks

# Data-generating model (log-odds scale)
intercept <- qlogis(.25)  # overall prevalence of 25% for an average item/person/country
recent_shift <- 1         # "recent" engagement is rarer than lifetime engagement
sd_item <- 1              # item random intercept SD
sd_respondent <- 1.2      # respondent random intercept SD
sd_country <- .4          # country random intercept SD
rho_gdp <- .8             # correlation between the capacity indicator ("herd") and GDP
b_gdp <- .3               # share of country-level variance explained by GDP

# ---- Simulation conditions ---------------------------------------------------

# size_sd is the SD of the log-normal weights that distribute respondents across
# countries: 1 = moderately uneven country sizes, 2 = very uneven (a few large
# countries, many with only a handful of respondents).

# H1: effect of publication pressure (log-odds per SD) at the minimum (1,000),
# an intermediate (1,250), and the target (1,500) sample size, including the
# smallest effect size of interest (SESOI, odds ratio = 1.15 per SD).
# b_pressure = 0 is included to check the false-positive rate.
# The country-level herd effect is switched off here (r = 0).
h1_conditions <- expand_grid(
  hypothesis = "H1",
  k = 30,
  r = 0,
  n_total = c(1000, 1250, 1500),
  size_sd = 1,
  b_pressure = c(0, .05, .1, log(1.15), .15, .2),  # log(1.15) = SESOI (odds ratio of 1.15 per SD)
  sim = seq_len(n_sims)
)

# H2: country-level partial correlation (net of GDP) at different numbers of
# countries, at the minimum (1,000) and the target (1,500) sample size.
# Publication pressure has no effect here.
h2_conditions <- 
  expand_grid(hypothesis = "H2",
              k = c(25, 30, 40),
              r = c(0, .3, .43, .49, .54, .6),
              n_total = c(1000, 1500),
              size_sd = 1,
              b_pressure = 0,
              sim = seq_len(n_sims)
              )

# H2 sensitivity: very uneven country sizes at the target sample size
h2_skew_conditions <- 
  h2_conditions |>
  filter(n_total == 1500) |>
  mutate(size_sd = 2)

# Each row is one replicate. `full` and `null` are the fixed-effects parts of the
# two models compared by the likelihood-ratio test, `term` is the tested effect.
conditions <- 
  bind_rows(h1_conditions, h2_conditions, h2_skew_conditions) |>
  mutate(full = if_else(hypothesis == "H1", "pressure", "herd + gdp"),
         null = if_else(hypothesis == "H1", "1", "gdp"),
         term = if_else(hypothesis == "H1", "pressure", "herd")
  )

# ---- Parallel setup ----------------------------------------------------------

plan(multisession, workers = 32)

# ---- Step 1: simulate one dataset per replicate -------------------------------

# Each dataset is stored in the list-column `sim_data`, nested by outcome
# (one row for "lifetime", one for "recent"). Unnesting that column below gives
# one row per replicate x outcome, with the respondent-level data nested in `data`.
simulated <- 
  conditions |>
  mutate(
    sim_data = future_pmap(
      pick(k, r, n_total, b_pressure, size_sd),
      \(k, r, n_total, b_pressure, size_sd) {

        # Countries: a capacity indicator (herd) that correlates with GDP, and a
        # country random intercept that depends on GDP and (for H2) on the part
        # of herd that is unrelated to GDP, with partial correlation r
        countries <- 
          tibble(country = seq_len(k),
                 gdp = rnorm(k),
                 herd_unique = rnorm(k),
                 herd = rho_gdp * gdp + sqrt(1 - rho_gdp^2) * herd_unique,
                 u_country = sd_country * (b_gdp * gdp + sqrt(1 - b_gdp^2) * 
                                             (r * herd_unique + sqrt(1 - r^2) * rnorm(k)))
        ) |>
          # Standardize the predictors so effects are per SD
          mutate(across(c(gdp, herd), \(x) as.numeric(scale(x)))) |>
          # Unequal country sizes: everyone gets the minimum, the rest is
          # distributed with log-normal weights
          mutate(
            n = min_per_country +
              as.vector(rmultinom(1, n_total - k * min_per_country, rlnorm(k, sdlog = size_sd)))
          )

        # Items: each belongs to a block and has its own random intercept
        items <- tibble(
          item = seq_len(n_items),
          block = rep(seq_len(n_blocks), length.out = n_items),
          u_item = rnorm(n_items, 0, sd_item)
        )

        countries |>
          # One row per respondent
          uncount(n) |>
          mutate(
            respondent = row_number(),
            pressure = rnorm(n()),
            u_respondent = rnorm(n(), 0, sd_respondent),
            blocks_done = sample(n_blocks, n(), replace = TRUE, prob = p_blocks)
          ) |>
          # Cross with items, then drop the blocks a respondent did not complete
          crossing(items) |>
          filter(block <= blocks_done) |>
          mutate(
            # Same uniform draw for both outcomes, so anyone with a recent
            # engagement also has a lifetime engagement
            eta = intercept + u_item + u_respondent + u_country + b_pressure * pressure,
            u = runif(n()),
            lifetime = u < plogis(eta),
            recent = u < plogis(eta - recent_shift)
          ) |>
          # Collapse to the number of engaged items out of the items answered
          summarise(
            across(c(lifetime, recent), sum),
            n_answered = n(),
            .by = c(country, respondent, gdp, herd, pressure)
          ) |>
          pivot_longer(c(lifetime, recent), names_to = "outcome", values_to = "events") |>
          nest(data = -outcome)
      },
      .options = furrr_options(seed = 20261012),
      .progress = TRUE
    )
  ) |>
  unnest(sim_data)

# ---- Step 2: fit the models for every replicate x outcome ---------------------

# The fit is added with `$<-` instead of mutate() on purpose: a function defined
# inside mutate() captures all columns (including the huge `data` column) and
# furrr would copy them to every worker.
# possibly() returns NA instead of stopping the whole run if a model fails.
simulated$fit <- future_pmap(
  select(simulated, data, full, null, term),
  possibly(
    \(data, full, null, term) {

      random <- c("(1 | respondent)", "(1 | country)")

      # Binomial GLMM: successes out of the items each respondent answered
      fit_full <- glmer(
        reformulate(c(full, random), response = quote(cbind(events, n_answered - events))),
        data = data,
        family = binomial,
        control = glmerControl(calc.derivs = FALSE)
      )
      fit_null <- glmer(
        reformulate(c(null, random), response = quote(cbind(events, n_answered - events))),
        data = data,
        family = binomial,
        control = glmerControl(calc.derivs = FALSE)
      )

      # Likelihood-ratio test of the tested term
      tibble(estimate = fixef(fit_full)[[term]],
             p = anova(fit_null, fit_full)[["Pr(>Chisq)"]][2]
             )
    },
    otherwise = tibble(estimate = NA_real_, p = NA_real_)
  ),
  .options = furrr_options(seed = 20261013),
  .progress = TRUE
)

# Drop the simulated data (large) and unpack the test results
results <- 
  simulated |>
  select(-data) |>
  unnest(fit)

# ---- Power -------------------------------------------------------------------

# A replicate is a success if both outcomes are significant in the same direction
# (failed fits give NA and count as not significant)
power <- 
  results |>
  summarise(significant = isTRUE(all(p < alpha & estimate > 0) | all(p < alpha & estimate < 0)),
            .by = c(hypothesis, k, r, n_total, size_sd, b_pressure, sim)) |>
  summarise(power = mean(significant),
            mcse = sqrt(power * (1 - power) / n()),  # Monte Carlo standard error
            n_sims = n(),
            .by = c(hypothesis, k, r, n_total, size_sd, b_pressure)
            )

write_csv(results, "output/power_replicates.csv")
write_csv(power, "output/power.csv")

# ---- Plots -------------------------------------------------------------------

power_h1_plot <- 
  power |>
  filter(hypothesis == "H1") |>
  ggplot(aes(b_pressure, power, colour = factor(n_total))) +
  geom_hline(yintercept = .8, linetype = "dashed") +
  geom_hline(yintercept = .9, linetype = "dotted") +
  geom_line() +
  geom_pointrange(aes(ymin = power - 1.96 * mcse, ymax = power + 1.96 * mcse)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.1), breaks = seq(0, 1, .2)) +
  labs(x = "Effect of publication pressure (log-odds per SD)",
       y = "Power (both outcomes significant)",
       tag = "A",
       colour = "Respondents") +
  theme(legend.position = "bottom")

power_h2_plot <- 
  power |>
  filter(hypothesis == "H2", size_sd == 1, n_total == 1500) |>
  ggplot(aes(r, power, colour = factor(k))) +
  geom_hline(yintercept = .8, linetype = "dashed") +
  geom_hline(yintercept = .9, linetype = "dotted") +
  geom_line() +
  geom_pointrange(aes(ymin = power - 1.96 * mcse, ymax = power + 1.96 * mcse)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.1), breaks = seq(0, 1, .2)) +
  guides() +
  labs(x = "Country-level partial correlation (net of GDP)",
       y = "Power (both outcomes significant)",
       tag = "B",
       colour = "Countries") +
  theme(legend.position = "bottom")

power_plot <- power_h1_plot + power_h2_plot

power_plot

ggsave("output/power.png", power_plot, width = 6, height = 4, dpi = 300)

