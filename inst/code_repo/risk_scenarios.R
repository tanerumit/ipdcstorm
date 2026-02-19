##############################################################
# 1. Poisson probability function
##############################################################

# Calculate P(exactly k events) for a Poisson distribution
poisson_prob <- function(k, lambda) {
  (lambda^k * exp(-lambda)) / factorial(k)
}

##############################################################
# 2. Single Event Probabilities by Severity Class
##############################################################

single_event_probs <- lambda_table %>%
  mutate(
    # Probability of exactly 0, 1, 2, 3 events
    p_exactly_0 = poisson_prob(0, lambda),
    p_exactly_1 = poisson_prob(1, lambda),
    p_exactly_2 = poisson_prob(2, lambda),
    p_exactly_3 = poisson_prob(3, lambda),
    
    # Check: these should sum close to 1 for low lambda
    # (for high lambda, need more terms)
    sum_check = p_exactly_0 + p_exactly_1 + p_exactly_2 + p_exactly_3
  )

cat("\n=== SINGLE EVENT PROBABILITIES ===\n")
cat("Probability of exactly k events occurring in any given year:\n\n")
print(single_event_probs %>% 
        select(severity, lambda, starts_with("p_exactly")))


##############################################################
# 3. Multiple Event Probabilities (Same Class)
##############################################################

multiple_event_probs <- lambda_table %>%
  mutate(
    # Probability of 2 or more events
    p_2_or_more = 1 - poisson_prob(0, lambda) - poisson_prob(1, lambda),
    
    # Probability of 3 or more events
    p_3_or_more = 1 - poisson_prob(0, lambda) - poisson_prob(1, lambda) - 
      poisson_prob(2, lambda),
    
    # Expected number of years (out of 100) with 2+ events
    years_with_2plus_per_century = p_2_or_more * 100,
    
    # Return period for 2+ events (years between occurrences)
    return_period_2plus = ifelse(p_2_or_more > 0, 1 / p_2_or_more, Inf)
  )

cat("\n=== MULTIPLE EVENT PROBABILITIES (SAME CLASS) ===\n")
cat("Probability of experiencing multiple events of the SAME severity in one year:\n\n")
print(multiple_event_probs %>% 
        select(severity, lambda, p_2_or_more, p_3_or_more, 
               return_period_2plus, years_with_2plus_per_century))


##############################################################
# 4. Return Periods
##############################################################

return_periods <- lambda_table %>%
  mutate(
    # Return period = 1 / annual probability
    return_period_years = ifelse(p_at_least_one > 0, 1 / p_at_least_one, Inf),
    
    # Alternative interpretation: average years between events
    avg_years_between = ifelse(lambda > 0, 1 / lambda, Inf),
    
    # Probability of experiencing event in next N years
    p_in_10_years = 1 - (1 - p_at_least_one)^10,
    p_in_50_years = 1 - (1 - p_at_least_one)^50,
    p_in_100_years = 1 - (1 - p_at_least_one)^100
  )

cat("\n=== RETURN PERIODS ===\n")
cat("Average time between events and cumulative probabilities:\n\n")
print(return_periods %>% 
        select(severity, lambda, return_period_years, 
               p_in_10_years, p_in_50_years, p_in_100_years))


##############################################################
# 5. Joint Probabilities (Multiple Classes)
##############################################################

# Create all pairwise combinations of severity classes
severity_levels <- lambda_table$severity

# Function to calculate joint probability assuming independence
calc_joint_prob <- function(sev1, sev2, lambda_tbl) {
  p1 <- lambda_tbl %>% filter(severity == sev1) %>% pull(p_at_least_one)
  p2 <- lambda_tbl %>% filter(severity == sev2) %>% pull(p_at_least_one)
  
  if(length(p1) == 0 || length(p2) == 0) return(NA)
  
  # P(both happen) = P(A) * P(B) [assuming independence]
  return(p1 * p2)
}

# Calculate all joint probabilities
joint_probs <- expand.grid(
  severity_1 = severity_levels,
  severity_2 = severity_levels,
  stringsAsFactors = FALSE
) %>%
  filter(severity_1 != severity_2) %>%  # Exclude same-same pairs
  rowwise() %>%
  mutate(
    p_both_occur = calc_joint_prob(severity_1, severity_2, lambda_table),
    return_period_years = ifelse(p_both_occur > 0, 1 / p_both_occur, Inf)
  ) %>%
  ungroup() %>%
  arrange(desc(p_both_occur))

cat("\n=== JOINT PROBABILITIES (DIFFERENT CLASSES) ===\n")
cat("Probability of experiencing events from MULTIPLE severity classes in same year:\n")
cat("⚠️  ASSUMES INDEPENDENCE between severity classes\n")
cat("   (may overestimate - active years tend to have more storms across all classes)\n\n")
print(joint_probs)


##############################################################
# 6. Empirical Validation from Monte Carlo
##############################################################

cat("\n=== MONTE CARLO VALIDATION ===\n")
cat("Comparing theoretical Poisson probabilities with empirical simulation results:\n\n")

# Calculate empirical probabilities from the sim_counts
empirical_probs <- sim_counts %>%
  group_by(severity) %>%
  summarise(
    # Empirical P(exactly k events)
    empirical_p_0 = mean(n_events == 0),
    empirical_p_1 = mean(n_events == 1),
    empirical_p_2 = mean(n_events == 2),
    empirical_p_at_least_1 = mean(n_events >= 1),
    empirical_p_2_or_more = mean(n_events >= 2),
    
    # Empirical lambda
    empirical_lambda = mean(n_events),
    .groups = "drop"
  )

# Compare with theoretical
validation_comparison <- lambda_table %>%
  select(severity, lambda, p_at_least_one, p_zero) %>%
  left_join(empirical_probs, by = "severity") %>%
  mutate(
    # Calculate theoretical values
    theoretical_p_1 = poisson_prob(1, lambda),
    theoretical_p_2 = poisson_prob(2, lambda),
    theoretical_p_2plus = 1 - p_zero - theoretical_p_1,
    
    # Differences (should be small for large n_sim)
    diff_lambda = abs(lambda - empirical_lambda),
    diff_p1 = abs(theoretical_p_1 - empirical_p_1),
    diff_p2plus = abs(theoretical_p_2plus - empirical_p_2_or_more)
  )

print(validation_comparison %>%
        select(severity, lambda, empirical_lambda, diff_lambda,
               p_at_least_one, empirical_p_at_least_1))


##############################################################
# 7. Practical Risk Scenarios
##############################################################

cat("\n=== PRACTICAL RISK SCENARIOS ===\n\n")

# Scenario 1: What's the risk over a planning horizon?
planning_horizon <- 30  # years (e.g., infrastructure lifespan)

scenario_1 <- lambda_table %>%
  mutate(
    p_at_least_one_in_period = 1 - (1 - p_at_least_one)^planning_horizon,
    expected_count_in_period = lambda * planning_horizon
  )

cat(sprintf("Over a %d-year planning horizon:\n", planning_horizon))
print(scenario_1 %>% 
        select(severity, p_at_least_one_in_period, expected_count_in_period))

# Scenario 2: Year with ANY major event (TS or worse)
p_any_tc <- 1 - prod(lambda_table$p_zero)
cat(sprintf("\nProbability of at least one TC event (any severity): %.3f\n", p_any_tc))

# Scenario 3: Particularly active year (multiple events any class)
cat("\nProbability of experiencing 2+ events of ANY severity in one year:\n")
active_year_prob <- sim_counts %>%
  group_by(sim_year) %>%
  summarise(total = sum(n_events)) %>%
  summarise(
    p_2plus = mean(total >= 2),
    p_3plus = mean(total >= 3),
    p_5plus = mean(total >= 5)
  )
print(active_year_prob)

# Scenario 4: Worst-case year (major hurricane AND something else)
major_tc_lambda <- lambda_table %>% 
  filter(severity == "Major_TC") %>% 
  pull(lambda)

if (length(major_tc_lambda) > 0 && major_tc_lambda > 0) {
  cat("\nWorst-case scenarios involving Major TC:\n")
  worst_case <- lambda_table %>%
    filter(severity != "Major_TC") %>%
    mutate(
      p_major_and_this = calc_joint_prob("Major_TC", severity, lambda_table),
      return_period = ifelse(p_major_and_this > 0, 1 / p_major_and_this, Inf)
    )
  print(worst_case %>% select(severity, p_major_and_this, return_period))
}


##############################################################
# 8. Summary Table for Reporting
##############################################################

summary_table <- lambda_table %>%
  mutate(
    # Key metrics
    annual_prob_pct = p_at_least_one * 100,
    return_period = ifelse(p_at_least_one > 0, 1 / p_at_least_one, Inf),
    p_2plus = 1 - poisson_prob(0, lambda) - poisson_prob(1, lambda),
    annual_prob_2plus_pct = p_2plus * 100,
    
    # Expected counts over different horizons
    expected_per_decade = lambda * 10,
    expected_per_30yrs = lambda * 30,
    
    # Cumulative risk
    p_in_30yrs = 1 - (1 - p_at_least_one)^30
  ) %>%
  select(severity, lambda, annual_prob_pct, return_period, 
         annual_prob_2plus_pct, expected_per_decade, 
         expected_per_30yrs, p_in_30yrs)

cat("\n=== COMPREHENSIVE SUMMARY TABLE ===\n")
cat("Complete probability metrics for each severity class:\n\n")
print(summary_table)


##############################################################
# DONE
##############################################################