# Power simulation for the confirmatory hypotheses of the ISRPP study
#   H1a: Publication pressure (respondent level) predicts QRP engagement
#   H1b: Having received training on QRPs (respondent level) predicts lower
#        QRP engagement
#   H2:  Country-level research capacity ("herd") predicts QRP engagement,
#        net of GDP
# The primary outcome is engagement in the past three years ("recent").
# Lifetime engagement is a secondary outcome. Power is also reported for the
# joint criterion (both outcomes significant in the same direction) for
# comparison with the earlier version of this analysis.
#
# Data are simulated and analyzed at the item level: QRP items are crossed
# with respondents, and respondents are nested in countries.

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
n_qrps <- 40           # QRP items (analyzed)
n_practices <- 60      # all practices shown (40 QRPs + 20 good practices)
block_size <- 15       # practices per block; blocks are formed after randomizing
                       # the order of all 60 practices (no stratification)
p_blocks <- c(.15, .15, .15, .55)  # probability of finishing 1, 2, 3, or all 4 blocks

# Data-generating model (log-odds scale)
# The intercept is calibrated so that the marginal lifetime prevalence, averaged
# over items, respondents, and countries, is 25%. This is somewhat lower than the
# mean item-level admission rates for the John et al. (2012) items (27-37%; John
# et al., 2012; Agnoli et al., 2017), as most Bestiary QRPs are expected to be rarer.
target_prevalence <- .25
recent_shift <- 1         # "recent" engagement is rarer than lifetime engagement
sd_item <- 1              # item random intercept SD
sd_respondent <- 1.2      # respondent random intercept SD
sd_country <- .4          # country random intercept SD
sd_slope_country <- .05   # SD of the country-specific deviations of the
                          # pressure and training slopes (assumption)
rho_gdp <- .8             # correlation between the capacity indicator ("herd") and GDP
b_gdp <- .3               # share of country-level variance explained by GDP
p_training <- .65         # share of respondents with any QRP training (assumption)

sd_total <- sqrt(sd_item^2 + sd_respondent^2 + sd_country^2)
intercept <- uniroot(
  \(b0) integrate(\(u) plogis(b0 + u) * dnorm(u, 0, sd_total), -Inf, Inf)$value - target_prevalence,
  interval = c(-10, 5)
)$root

# ---- Simulation conditions ---------------------------------------------------

# size_sd is the SD of the log-normal weights that distribute respondents across
# countries: 1 = moderately uneven country sizes, 2 = very uneven (a few large
# countries, many with only a handful of respondents).
# Effects that are switched off in a condition are set to 0.

# H1a: effect of publication pressure (log-odds per SD) at the minimum (1,000),
# an intermediate (1,250), and the target (1,500) sample size, at the minimum
# number of countries. b_pressure = 0 checks the false-positive rate.
h1a_conditions <-
  expand_grid(
    hypothesis = "H1a",
    k = 30,
    r = 0,
    n_total = c(1000, 1250, 1500),
    size_sd = 1,
    b_pressure = c(0, .05, .1, log(1.15), .15, .2),  # log(1.15) = SESOI (odds ratio of 1.15 per SD)
    b_training = 0,
    sim = seq_len(n_sims)
  )

# H1b: effect of any QRP training vs. none (log-odds) at the same sample sizes.
# Training is coded as any training reported vs. "None".
h1b_conditions <-
  expand_grid(
    hypothesis = "H1b",
    k = 30,
    r = 0,
    n_total = c(1000, 1250, 1500),
    size_sd = 1,
    b_pressure = 0,
    b_training = log(c(1, .9, .85, .8, .75, .7)),  # odds ratios of 1 to .70
    sim = seq_len(n_sims)
  )

# H2: country-level partial correlation (net of GDP) at different numbers of
# countries, at the minimum (1,000) and the target (1,500) sample size.
h2_conditions <-
  expand_grid(
    hypothesis = "H2",
    k = c(30, 35, 40),
    r = c(0, .3, .43, .5, .6),
    n_total = c(1000, 1500),
    size_sd = 1,
    b_pressure = 0,
    b_training = 0,
    sim = seq_len(n_sims)
  )

# H2 sensitivity: very uneven country sizes at the target sample size
h2_skew_conditions <-
  h2_conditions |>
  filter(n_total == 1500) |>
  mutate(size_sd = 2)

# `term` is the tested effect, `full` and `null` are the fixed-effects parts of
# the two models compared by the likelihood-ratio test, and `slope` names the
# predictor that gets a random slope across countries (none for H2).
# `effect` is the true effect in each condition and sets the expected direction.
conditions <-
  bind_rows(h1a_conditions, h1b_conditions, h2_conditions, h2_skew_conditions) |>
  mutate(
    term = case_match(hypothesis, "H1a" ~ "pressure", "H1b" ~ "training", "H2" ~ "herd"),
    full = case_match(hypothesis, "H2" ~ "herd + gdp", .default = term),
    null = case_match(hypothesis, "H2" ~ "gdp", .default = "1"),
    slope = case_match(hypothesis, "H2" ~ NA_character_, .default = term),
    effect = case_match(hypothesis, "H1a" ~ b_pressure, "H1b" ~ b_training, "H2" ~ r)
  )

# ---- Parallel setup ----------------------------------------------------------

plan(multisession, workers = 32)

# ---- Simulate and fit every replicate ----------------------------------------

# Data are simulated and analyzed in the same worker call, so the item-level
# datasets never have to be kept in memory.
# The results are added with `$<-` instead of mutate() on purpose: a function
# defined inside mutate() captures all columns and furrr would copy them to every
# worker.
# possibly() returns NA instead of stopping the whole run if a replicate fails.
conditions$fit <- future_pmap(
  select(conditions, k, r, n_total, size_sd, b_pressure, b_training, term, full, null, slope),
  possibly(
    \(k, r, n_total, size_sd, b_pressure, b_training, term, full, null, slope) {

      # Countries: a capacity indicator (herd) that correlates with GDP, and a
      # country random intercept that depends on GDP and (for H2) on the part
      # of herd that is unrelated to GDP, with partial correlation r
      countries <-
        tibble(
          country = seq_len(k),
          gdp = rnorm(k),
          herd_unique = rnorm(k),
          herd = rho_gdp * gdp + sqrt(1 - rho_gdp^2) * herd_unique,
          u_country = sd_country * (b_gdp * gdp + sqrt(1 - b_gdp^2) *
                                      (r * herd_unique + sqrt(1 - r^2) * rnorm(k))),
          v_pressure = rnorm(k, 0, sd_slope_country),
          v_training = rnorm(k, 0, sd_slope_country)
        ) |>
        # Standardize the predictors so effects are per SD
        mutate(across(c(gdp, herd), \(x) as.numeric(scale(x)))) |>
        # Unequal country sizes: everyone gets the minimum, the rest is
        # distributed with log-normal weights
        mutate(
          n = min_per_country +
            as.vector(rmultinom(1, n_total - k * min_per_country, rlnorm(k, sdlog = size_sd)))
        )

      items <- tibble(item = seq_len(n_qrps), u_item = rnorm(n_qrps, 0, sd_item))

      sim_data <-
        countries |>
        # One row per respondent
        uncount(n) |>
        mutate(
          respondent = row_number(),
          pressure = rnorm(n()),
          training = rbinom(n(), 1, p_training),
          u_respondent = rnorm(n(), 0, sd_respondent),
          blocks_done = sample(length(p_blocks), n(), replace = TRUE, prob = p_blocks)
        ) |>
        # Cross with the QRP items. Each respondent sees all 60 practices in a
        # random order, so the 40 QRPs take 40 random positions out of 60; a QRP
        # is answered if its block was reached.
        crossing(items) |>
        mutate(position = sample(n_practices, n()), .by = respondent) |>
        filter(ceiling(position / block_size) <= blocks_done) |>
        mutate(
          # Same uniform draw for both outcomes, so anyone with a recent
          # engagement also has a lifetime engagement
          eta = intercept + u_item + u_respondent + u_country +
            (b_pressure + v_pressure) * pressure +
            (b_training + v_training) * training,
          u = runif(n()),
          lifetime = as.integer(u < plogis(eta)),
          recent = as.integer(u < plogis(eta - recent_shift))
        ) |>
        select(country, respondent, item, gdp, herd, pressure, training, lifetime, recent)

      # Random-effect structures in the order they are tried: with a random slope
      # for the tested respondent-level predictor across countries, then without
      # it if that model is singular or does not converge (as in the analysis plan)
      random_sets <-
        list(
          with_slope = c(str_glue("(1 + {slope} | country)"), "(1 | country:respondent)", "(1 | item)"),
          intercepts = c("(1 | country)", "(1 | country:respondent)", "(1 | item)")
        ) |>
        keep_at(if (is.na(slope)) "intercepts" else c("with_slope", "intercepts"))

      c("lifetime", "recent") |>
        set_names() |>
        map(\(outcome) {

          fit_pair <- \(random) {
            fits <-
              list(full = full, null = null) |>
              map(\(fixed) glmer(
                reformulate(c(fixed, random), response = outcome),
                data = sim_data,
                family = binomial,
                control = glmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
              ))
            fits$ok <- fits |>
              map_lgl(\(fit) !isSingular(fit) && length(fit@optinfo$conv$lme4$messages) == 0) |>
              all()
            fits
          }

          fits <- fit_pair(random_sets[[1]])
          slope_dropped <- FALSE
          if (!fits$ok && length(random_sets) > 1) {
            fits <- fit_pair(random_sets[[2]])
            slope_dropped <- TRUE
          }

          lrt <- anova(fits$null, fits$full)

          tibble(
            estimate = fixef(fits$full)[[term]],
            chisq = lrt$Chisq[2],
            p = lrt$`Pr(>Chisq)`[2],
            converged = fits$ok,
            slope_dropped = slope_dropped
          )
        }) |>
        list_rbind(names_to = "outcome")
    },
    otherwise = tibble(outcome = c("lifetime", "recent"), estimate = NA_real_,
                       chisq = NA_real_, p = NA_real_, converged = FALSE,
                       slope_dropped = NA)
  ),
  .options = furrr_options(seed = 20261109),
  .progress = TRUE
)

results <-
  conditions |>
  unnest(fit)

# ---- Size-adjusted test for H2 -----------------------------------------------

# The planned test of H2 is a parametric-bootstrap likelihood-ratio test, which
# keeps the false-positive rate at alpha with few countries. Running a bootstrap
# inside every replicate is not feasible, so its power is approximated with a
# size-adjusted test: the critical value is the 95th percentile of the
# likelihood-ratio statistic in the replicates without an effect (same
# hypothesis, number of countries, sample size, country-size spread, and
# outcome). This is applied to all hypotheses; for H1a and H1b it should be
# close to the asymptotic test.
critical_values <-
  results |>
  filter(effect == 0) |>
  summarise(critical = quantile(chisq, 1 - alpha, na.rm = TRUE),
            .by = c(hypothesis, k, n_total, size_sd, outcome))

results <-
  results |>
  left_join(critical_values, by = join_by(hypothesis, k, n_total, size_sd, outcome)) |>
  mutate(
    # A replicate is a success if the test is significant in the direction of the
    # true effect; without a true effect, any significant result is a false positive.
    right_direction = if_else(effect == 0, TRUE, sign(estimate) == sign(effect)),
    sig_asymptotic = coalesce(p < alpha & right_direction, FALSE),
    sig_adjusted = coalesce(chisq > critical & right_direction, FALSE)
  )

# ---- Power -------------------------------------------------------------------

condition_vars <- c("hypothesis", "k", "r", "n_total", "size_sd", "b_pressure", "b_training", "effect")

# Power for each outcome separately and for the joint criterion
# (both outcomes significant). Failed fits count as not significant.
power <-
  results |>
  select(all_of(condition_vars), sim, outcome, sig_asymptotic, sig_adjusted) |>
  bind_rows(
    results |>
      summarise(across(c(sig_asymptotic, sig_adjusted), all),
                .by = c(all_of(condition_vars), sim)) |>
      mutate(outcome = "both")
  ) |>
  summarise(
    across(c(sig_asymptotic, sig_adjusted), mean, .names = "{str_remove(.col, 'sig_')}"),
    n_sims = n(),
    .by = c(all_of(condition_vars), outcome)
  ) |>
  pivot_longer(c(asymptotic, adjusted), names_to = "test", values_to = "power") |>
  mutate(mcse = sqrt(power * (1 - power) / n_sims))  # Monte Carlo standard error

# Convergence: share of fits that converged without a singular fit, and share
# of fits in which the random slope had to be dropped
convergence <-
  results |>
  summarise(
    converged = mean(converged),
    slope_dropped = mean(slope_dropped, na.rm = TRUE),
    failed = mean(is.na(p)),
    .by = c(hypothesis, k, n_total, size_sd, outcome)
  )

write_csv(results, "output/power_replicates.csv")
write_csv(power, "output/power.csv")
write_csv(convergence, "output/power_convergence.csv")

# False-positive rates of the asymptotic test (rows with effect == 0)
power |>
  filter(effect == 0, test == "asymptotic") |>
  select(hypothesis, k, n_total, size_sd, outcome, false_positive_rate = power) |>
  print(n = Inf)

# ---- Plots -------------------------------------------------------------------

# Primary outcome (recent engagement). H1a and H1b use the asymptotic test,
# H2 the size-adjusted test that approximates the parametric bootstrap.
y_scale <- scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.1), breaks = seq(0, 1, .2))
power_lines <- list(geom_hline(yintercept = .8, linetype = "dashed"),
                    geom_hline(yintercept = .9, linetype = "dotted"))

power_h1a_plot <-
  power |>
  filter(hypothesis == "H1a", outcome == "recent", test == "asymptotic") |>
  ggplot(aes(b_pressure, power, colour = factor(n_total))) +
  power_lines +
  geom_line() +
  geom_pointrange(aes(ymin = power - 1.96 * mcse, ymax = power + 1.96 * mcse)) +
  y_scale +
  labs(x = "Effect of publication pressure (log-odds per SD)",
       y = "Power",
       tag = "A",
       colour = "Respondents") +
  theme(legend.position = "bottom")

power_h1b_plot <-
  power |>
  filter(hypothesis == "H1b", outcome == "recent", test == "asymptotic") |>
  ggplot(aes(exp(b_training), power, colour = factor(n_total))) +
  power_lines +
  geom_line() +
  geom_pointrange(aes(ymin = power - 1.96 * mcse, ymax = power + 1.96 * mcse)) +
  scale_x_reverse() +
  y_scale +
  labs(x = "Odds ratio of QRP training (any vs. none)",
       y = "Power",
       tag = "B",
       colour = "Respondents") +
  theme(legend.position = "bottom")

power_h2_plot <-
  power |>
  filter(hypothesis == "H2", outcome == "recent", test == "adjusted",
         size_sd == 1, n_total == 1500) |>
  ggplot(aes(r, power, colour = factor(k))) +
  power_lines +
  geom_line() +
  geom_pointrange(aes(ymin = power - 1.96 * mcse, ymax = power + 1.96 * mcse)) +
  y_scale +
  labs(x = "Country-level partial correlation (net of GDP)",
       y = "Power",
       tag = "C",
       colour = "Countries") +
  theme(legend.position = "bottom")

power_plot <- power_h1a_plot + power_h1b_plot + power_h2_plot

power_plot

ggsave("output/power.png", power_plot, width = 9, height = 4, dpi = 300)
