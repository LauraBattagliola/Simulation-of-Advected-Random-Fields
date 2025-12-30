# 04_DampedFF_SpectralFilter.R
# Damped Frozen Field spectral filter visualization (Figure 7)
#
# Spectrum: S_ZZ(k, omega) = S_XX(k) * B(k, omega; v, beta)^{-alpha}
# where B(k, omega; v, beta) = (omega + k'v)^2 + (beta*omega)^2

library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)

# Frequency grids
omega_grid <- seq(-pi, pi, length.out = 400)
k1_grid <- seq(-pi/4, pi/4, length.out = 400)

# Velocity
v <- c(5, 0)

# Parameter values
alpha_vals <- c(0.5, 1.0, 1.5)  # persistence (must be > 1/2)
beta_vals <- c(0, 1.0, 2.0)     # damping

param_grid <- expand_grid(alpha = alpha_vals, beta = beta_vals)

# Compute filter B(k, omega; v, beta)^{-alpha} on the k2 = 0 plane
calculate_filter <- function(alpha_val, beta_val) {
  df <- expand_grid(omega = omega_grid, k1 = k1_grid)
  # B = (omega + k1*v1)^2 + (beta*omega)^2
  B <- (df$omega + df$k1 * v[1])^2 + (beta_val * df$omega)^2
  df$log_value <- -alpha_val * log10(B + 1e-9)
  return(df)
}

all_data <- map2_dfr(param_grid$alpha, param_grid$beta,
                     ~ calculate_filter(.x, .y) %>% mutate(alpha = .x, beta = .y))

all_data <- all_data %>%
  mutate(
    alpha_label = factor(paste("alpha ==", alpha), levels = paste("alpha ==", alpha_vals)),
    beta_label = factor(paste("beta ==", beta), levels = paste("beta ==", beta_vals))
  )

# Plot
final_plot <- ggplot(all_data, aes(x = k1, y = omega)) +
  geom_raster(aes(fill = log_value)) +
  scale_fill_viridis_c(name = expression(log[10](B(beta)^-alpha)),
                       breaks = c(-2.5, 0.0, 2.5, 5.0)) +
  geom_contour(aes(z = log_value, color = after_stat(level)),
               linewidth = 0.4, alpha = 0.8) +
  scale_color_gradient(low = "deeppink1", high = "white", guide = "none") +
  # Advection line: omega = -k'v
  geom_abline(intercept = 0, slope = -v[1], color = "black",
              linetype = "dashed", linewidth = 0.8) +
  facet_grid(rows = vars(alpha_label), cols = vars(beta_label),
             labeller = label_parsed) +
  labs(x = expression(k[1]), y = expression(omega)) +
  scale_y_continuous(breaks = c(-2.5, 0, 2.5)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(4, "cm"),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    panel.grid = element_blank()
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

print(final_plot)

pdf(file = "04_DampedFF_SpectralFilter.pdf", width = 12, height = 9)
print(final_plot)
dev.off()
