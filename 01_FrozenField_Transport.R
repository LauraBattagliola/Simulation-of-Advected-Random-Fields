# 01_FrozenField_Transport.R
# Frozen Field simulation and comparison with transport PDE solver (Figure 3)

library(ReacTran)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Grid setup
Lx <- 100; Ly <- 100
Nx <- 100; Ny <- 100

# Velocity v and diffusion (zero for pure transport)
v <- c(20, 10)
Dxx <- 0.0
Dyy <- 0.0

# Time
t_final <- 3.0
plot_times <- c(0.0, 1.0, 2.0, 3.0)
times_to_solve <- seq(0, t_final, by = 0.1)

grid.x_1D <- setup.grid.1D(x.up = 0, L = Lx, N = Nx)
grid.y_1D <- setup.grid.1D(x.up = 0, L = Ly, N = Ny)
grid_2D <- setup.grid.2D(x.grid = grid.x_1D, y.grid = grid.y_1D)

x_coords <- grid.x_1D$x.mid
y_coords <- grid.y_1D$x.mid

# Initial condition X_S(s)
initial_condition <- function(s1, s2) {
  s1_0 <- 30; s2_0 <- 30
  sigma1 <- 8; sigma2 <- 8
  exp(-((s1 - s1_0)^2 / (2 * sigma1^2) + (s2 - s2_0)^2 / (2 * sigma2^2)))
}

C_initial <- outer(x_coords, y_coords, FUN = initial_condition)

# PDE solver for transport equation
pde_model <- function(t, C, params) {
  CONC <- matrix(nrow = Nx, ncol = Ny, data = C)
  flux <- tran.2D(C = CONC, v.x = v[1], v.y = v[2], D.x = Dxx, D.y = Dyy, grid = grid_2D)
  return(list(flux$dC))
}

pde_result <- ode.2D(y = C_initial, times = times_to_solve,
                     func = pde_model, parms = NULL, 
                     dimens = c(Nx, Ny), lrw = 700000)

# Frozen Field: Z(s,t) = X_S(s - vt)
s1_matrix <- matrix(x_coords, nrow = Nx, ncol = Ny, byrow = FALSE)
s2_matrix <- matrix(y_coords, nrow = Nx, ncol = Ny, byrow = TRUE)

format_for_ggplot <- function(mat, x_vec, y_vec) {
  rownames(mat) <- x_vec
  colnames(mat) <- y_vec
  mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("s1") %>%
    pivot_longer(-s1, names_to = "s2", values_to = "concentration") %>%
    mutate(s1 = as.numeric(s1), s2 = as.numeric(s2))
}

# Extract PDE results
pde_data <- lapply(plot_times, function(t) {
  idx <- which.min(abs(pde_result[, "time"] - t))
  mat <- matrix(nrow = Nx, ncol = Ny, data = pde_result[idx, -1])
  df <- format_for_ggplot(mat, x_coords, y_coords)
  df$time <- t
  df
})
pde_plot_data <- bind_rows(pde_data)
pde_plot_data$Model <- "PDE Solver"

# Compute FF states
ff_data <- lapply(plot_times, function(t) {
  shifted_s1 <- s1_matrix - v[1] * t
  shifted_s2 <- s2_matrix - v[2] * t
  ff_state <- initial_condition(shifted_s1, shifted_s2)
  df <- format_for_ggplot(ff_state, x_coords, y_coords)
  df$time <- t
  df
})
ff_plot_data <- bind_rows(ff_data)
ff_plot_data$Model <- "FF"

all_plot_data <- bind_rows(pde_plot_data, ff_plot_data)

# RMSE between PDE and FF
rmse_data <- lapply(plot_times, function(t) {
  idx <- which.min(abs(pde_result[, "time"] - t))
  pde_mat <- matrix(nrow = Nx, ncol = Ny, data = pde_result[idx, -1])
  
  shifted_s1 <- s1_matrix - v[1] * t
  shifted_s2 <- s2_matrix - v[2] * t
  ff_mat <- initial_condition(shifted_s1, shifted_s2)
  
  rmse <- sqrt(mean((pde_mat - ff_mat)^2))
  data.frame(time = t, rmse = rmse)
})
rmse_plot_data <- bind_rows(rmse_data)

# Plotting
all_plot_data <- all_plot_data %>%
  mutate(
    time_label = factor(paste0("Z(bold(s),~t==", time, ")"),
                        levels = paste0("Z(bold(s),~t==", plot_times, ")")),
    model_label = factor(ifelse(Model == "PDE Solver", "italic('PDE Solver')", "italic('FF')"),
                         levels = c("italic('PDE Solver')", "italic('FF')"))
  )

heatmap_plot <- ggplot(all_plot_data, aes(x = s1, y = s2, fill = concentration)) +
  geom_raster() +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  coord_fixed(ratio = 1) +
  facet_grid(model_label ~ time_label, labeller = label_parsed) +
  labs(x = expression(s[1]), y = expression(s[2]), fill = "Concentration") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2.5, "cm"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    strip.background = element_blank()
  )

rmse_plot <- ggplot(rmse_plot_data, aes(x = time, y = rmse)) +
  geom_line(color = "navy", linewidth = 1) +
  geom_point(color = "navy", size = 3) +
  labs(x = "t", y = "RMSE(t)") +
  theme_bw() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = plot_times)

final_plot <- heatmap_plot + rmse_plot + plot_layout(widths = c(5, 2))

pdf(file = "01_FrozenField_Transport.pdf", width = 15, height = 8)
print(final_plot)
dev.off()
