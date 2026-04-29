library(dplyr)
library(ggplot2)
test_results <- absolute_Breakpoint_Accuracy %>%
  filter(seq_depth == "100GB", seq_id == "95%") %>%
  group_by(Algorithm, Breakpoint_number) %>%
  summarise(
    p_value = wilcox.test(Distance_from_Simulated_Breakpoint, mu = 0)$p.value,
    .groups = "drop"
  )
test_results <- test_results %>%
  mutate(
    signif = case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
breakpoint10095 +
  geom_text(
    data = test_results,
    aes(x = Breakpoint_number, y = 1e8, label = signif, color = Algorithm),
    position = position_dodge(width = 0.75),
    size = 6,
    inherit.aes = FALSE
  )

ggplot(filter(absolute_Breakpoint_Accuracy, seq_depth == "100GB", seq_id == "95%"), 
       aes(sample = Distance_from_Simulated_Breakpoint)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(title = "Normal Q-Q Plot", x = "Theoretical", y = "Sample")+
  facet_grid(Algorithm ~ Breakpoint_number)


install.packages("TOSTER")
library(TOSTER)

delta <- 100  # equivalence bound (bp)

tost_results <- absolute_Breakpoint_Accuracy %>%
  filter(seq_depth == "100GB", seq_id == "95%") %>%
  group_by(Algorithm, Breakpoint_number) %>%
  summarise(
    n = n(),
    mean_val = mean(Distance_from_Simulated_Breakpoint, na.rm = TRUE),
    sd_val = sd(Distance_from_Simulated_Breakpoint, na.rm = TRUE),
    
    tost = list(
      TOSTone.raw(
        m = mean_val,
        sd = sd_val,
        n = n,
        mu = 0,
        low_eqbound = -delta,
        high_eqbound = delta
      )
    ),
    .groups = "drop"
  )


##########################
# Try with bootstrap
##########################
library(dplyr)
library(purrr)

bootstrap_tost <- function(data, B = 5000, delta = 100) {
  
  boot_means <- replicate(B, {
    resample <- data[sample(nrow(data), replace = TRUE), ]
    mean(resample$Distance_from_Simulated_Breakpoint, na.rm = TRUE)
  })
  
  obs_mean <- mean(data$Distance_from_Simulated_Breakpoint, na.rm = TRUE)
  
  ci <- quantile(boot_means, c(0.025, 0.975))
  
  equivalent <- (ci[1] > -delta) & (ci[2] < delta)
  
  tibble::tibble(
    mean = obs_mean,
    ci_lower = ci[1],
    ci_upper = ci[2],
    equivalent = equivalent
  )
}
boot_results <- absolute_Breakpoint_Accuracy %>%
  filter(seq_depth == "100GB", seq_id == "95%") %>%
  group_by(Algorithm, Breakpoint_number) %>%
  group_modify(~ bootstrap_tost(.x, B = 5000, delta = delta))


####################################
# Remove outliers
####################################
library(dplyr)
library(TOSTER)

cleaned_data <- absolute_Breakpoint_Accuracy %>%
  filter(seq_depth == "100GB", seq_id == "95%") %>%
  group_by(Algorithm, Breakpoint_number) %>%
  filter({
    # Calculate IQR bounds for the current group
    Q1 <- quantile(Distance_from_Simulated_Breakpoint, 0.25, na.rm = TRUE)
    Q3 <- quantile(Distance_from_Simulated_Breakpoint, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    
    # Keep only values within 1.5 * IQR
    Distance_from_Simulated_Breakpoint >= (Q1 - 1.5 * IQR_val) & Distance_from_Simulated_Breakpoint <= (Q3 + 1.5 * IQR_val)
  }) %>%
  ungroup()

ggplot(cleaned_data, 
       aes(sample = Distance_from_Simulated_Breakpoint)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(title = "Normal Q-Q Plot", x = "Theoretical", y = "Sample")+
  facet_grid(Algorithm ~ Breakpoint_number)

delta <- 100  # equivalence bound (bp)

tost_results <- cleaned_data %>%
  group_by(Algorithm, Breakpoint_number) %>%
  summarise(
    n = n(),
    mean_val = mean(Distance_from_Simulated_Breakpoint, na.rm = TRUE),
    sd_val = sd(Distance_from_Simulated_Breakpoint, na.rm = TRUE),
    
    tost = list(
      TOSTone.raw(
        m = mean_val,
        sd = sd_val,
        n = n,
        mu = 0,
        low_eqbound = -delta,
        high_eqbound = delta
      )
    ),
    .groups = "drop"
  )

# Bootstrap to be sure:
library(dplyr)
library(purrr)

bootstrap_tost <- function(data, B = 5000, delta = 100) {
  
  boot_means <- replicate(B, {
    resample <- data[sample(nrow(data), replace = TRUE), ]
    mean(resample$Distance_from_Simulated_Breakpoint, na.rm = TRUE)
  })
  
  obs_mean <- mean(data$Distance_from_Simulated_Breakpoint, na.rm = TRUE)
  
  ci <- quantile(boot_means, c(0.025, 0.975))
  
  equivalent <- (ci[1] > -delta) & (ci[2] < delta)
  
  tibble::tibble(
    mean = obs_mean,
    ci_lower = ci[1],
    ci_upper = ci[2],
    equivalent = equivalent
  )
}
boot_results <- cleaned_data %>%
  group_by(Algorithm, Breakpoint_number) %>%
  group_modify(~ bootstrap_tost(.x, B = 5000, delta = delta))
#############################
# Using Median
#############################

tost_results <- absolute_Breakpoint_Accuracy %>%
  filter(seq_depth == "100GB", seq_id == "95%") %>%
  group_by(Algorithm, Breakpoint_number) %>%
  summarise(
    n = n(),
    median_val = median(Distance_from_Simulated_Breakpoint, na.rm = TRUE),
    .groups = "drop"
  )

bootstrap_median_tost <- function(data, B = 5000, delta = 100) {
  
  boot_medians <- replicate(B, {
    resample <- data[sample(nrow(data), replace = TRUE), ]
    median(resample$Distance_from_Simulated_Breakpoint, na.rm = TRUE)
  })
  
  obs_median <- median(data$Distance_from_Simulated_Breakpoint, na.rm = TRUE)
  
  ci <- quantile(boot_medians, c(0.025, 0.975))
  
  tibble::tibble(
    median = obs_median,
    ci_lower = ci[1],
    ci_upper = ci[2],
    equivalent = (ci[1] > -delta & ci[2] < delta)
  )
}

boot_median_results <- absolute_Breakpoint_Accuracy %>%
  filter(seq_depth == "100GB", seq_id == "95%") %>%
  group_by(Algorithm, Breakpoint_number) %>%
  group_modify(~ bootstrap_median_tost(.x, B = 10000, delta = delta))



absolute_Breakpoint_Accuracy %>%
  filter(seq_depth == "100GB", seq_id == "95%") %>%
  group_by(Algorithm, Breakpoint_number) %>% head()
