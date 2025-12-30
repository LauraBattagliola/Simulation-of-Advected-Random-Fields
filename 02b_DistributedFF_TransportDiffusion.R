# 02b_DistributedFF_TransportDiffusion.R
# Distributed FF comparison with advection-diffusion PDE (Figure 4)

library(ReacTran)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(MASS)

# Grid
Lx <- 100; Ly <- 100
Nx <- 100; Ny <- 100

# Mean velocity and diffusion tensor
v_bar <- c(20, 10)
Dxx <- 5.0
Dyy <- 1.5
D <- matrix(c(Dxx, 0, 0, Dyy), nrow = 2)

# Simulation settings
nv_list <- c(10, 50, 200)
n_iter <- 500
set.seed(123)

t_final <- 3.0
times_to_solve <- seq(0, t_final, by = 0.1)
plot_times <- c(0.0, 1.0, 2.0, 3.0)

# Grid setup
grid.x_1D <- setup.grid.1D(x.up = 0, L = Lx, N = Nx)
grid.y_1D <- setup.grid.1D(x.up = 0, L = Ly, N = Ny)
grid_2D <- setup.grid.2D(x.grid = grid.x_1D, y.grid = grid.y_1D)
x_coords <- grid.x_1D$x.mid
y_coords <- grid.y_1D$x.mid

s1_matrix <- matrix(x_coords, nrow = Nx, ncol = Ny, byrow = FALSE)
s2_matrix <- matrix(y_coords, nrow = Nx, ncol = Ny, byrow = TRUE)

# Initial condition
initial_condition <- function(s1, s2) {
  s1_0 <- 30; s2_0 <- 30
  sigma1 <- 8; sigma2 <- 8
  exp(-((s1 - s1_0)^2 / (2 * sigma1^2) + (s2 - s2_0)^2 / (2 * sigma2^2)))
}

C_initial <- outer(x_coords, y_coords, FUN = initial_condition)

# PDE solver (advection-diffusion)
pde_model <- function(t, C, params) {
  CONC <- matrix(nrow = Nx, ncol = Ny, data = C)
  flux <- tran.2D(C = CONC, v.x = v_bar[1], v.y = v_bar[2], 
                  D.x = Dxx, D.y = Dyy, grid = grid_2D)
  return(list(flux$dC))
}

pde_result <- ode.2D(y = C_initial, times = times_to_solve, func = pde_model,
                     parms = NULL, dimens = c(Nx, Ny), lrw = 700000)

pde_matrices <- lapply(plot_times, function(t) {
  idx <- which.min(abs(pde_result[, "time"] - t))
  matrix(data = pde_result[idx, -1], nrow = Nx, ncol = Ny)
})
names(pde_matrices) <- as.character(plot_times)

# Helper for plotting
format_for_ggplot <- function(mat, x_vec, y_vec) {
  mat %>%
    as.data.frame() %>%
    `rownames<-`(x_vec) %>%
    `colnames<-`(y_vec) %>%
    tibble::rownames_to_column("s1") %>%
    pivot_longer(-s1, names_to = "s2", values_to = "concentration") %>%
    mutate(s1 = as.numeric(s1), s2 = as.numeric(s2))
}

# Visual comparison data (nv = 50)
pde_plot_data <- bind_rows(lapply(plot_times, function(t) {
  format_for_ggplot(pde_matrices[[as.character(t)]], x_coords, y_coords) %>%
    mutate(time = t)
}))
pde_plot_data$Model <- "PDE Solver"

set.seed(42)
velocities_plot <- mvrnorm(n = 50, mu = v_bar, Sigma = D)
weights_plot <- rep(1/50, 50)

dff_matrices <- lapply(plot_times, function(t) {
  dff <- matrix(0, nrow = Nx, ncol = Ny)
  for (i in 1:50) {
    vi <- velocities_plot[i, ]
    dff <- dff + weights_plot[i] * initial_condition(s1_matrix - vi[1]*t, s2_matrix - vi[2]*t)
  }
  dff
})

dff_plot_data <- bind_rows(lapply(seq_along(plot_times), function(i) {
  format_for_ggplot(dff_matrices[[i]], x_coords, y_coords) %>% mutate(time = plot_times[i])
}))
dff_plot_data$Model <- "DFF"

all_plot_data <- bind_rows(pde_plot_data, dff_plot_data) %>%
  mutate(
    time_label = factor(paste0("Z(bold(s),~t==", time, ")"),
                        levels = paste0("Z(bold(s),~t==", plot_times, ")")),
    model_label = factor(ifelse(Model == "PDE Solver", "italic('PDE Solver')", "italic('DFF')"),
                         levels = c("italic('PDE Solver')", "italic('DFF')"))
  )

# Monte Carlo for RMSE
all_rmse <- data.frame()
set.seed(123)

for (nv in nv_list) {
  cat("Running nv =", nv, "\n")
  for (k in 1:n_iter) {
    velocities <- mvrnorm(n = nv, mu = v_bar, Sigma = D)
    weights <- rep(1/nv, nv)
    
    for (t in plot_times) {
      dff <- matrix(0, nrow = Nx, ncol = Ny)
      for (i in 1:nv) {
        vi <- velocities[i, ]
        dff <- dff + weights[i] * initial_condition(s1_matrix - vi[1]*t, s2_matrix - vi[2]*t)
      }
      pde_truth <- pde_matrices[[as.character(t)]]
      rmse <- sqrt(mean((dff - pde_truth)^2))
      all_rmse <- rbind(all_rmse, data.frame(nv = nv, iter = k, time = t, rmse = rmse))
    }
  }
}

# Plots
heatmap_plot <- ggplot(all_plot_data, aes(x = s1, y = s2, fill = concentration)) +
  geom_raster() +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  facet_grid(model_label ~ time_label, labeller = label_parsed) +
  labs(x = expression(s[1]), y = expression(s[2]), fill = "Concentration") +
  theme_bw() +
  theme(legend.position = "bottom", legend.key.width = unit(2.5, "cm"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.background = element_blank())

median_data <- all_rmse %>% group_by(nv, time) %>% summarise(median_rmse = median(rmse), .groups = 'drop')
dodge <- position_dodge(width = 0.8)

rmse_plot <- ggplot(all_rmse, aes(x = factor(time), y = rmse, color = factor(nv))) +
  geom_boxplot(aes(fill = factor(nv)), alpha = 0.3, outlier.shape = NA, position = dodge) +
  geom_line(data = median_data, aes(y = median_rmse, group = factor(nv)), position = dodge, linewidth = 1) +
  geom_point(data = median_data, aes(y = median_rmse), position = dodge, size = 3) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "t", y = "RMSE(t)", color = expression(n[v]), fill = expression(n[v])) +
  theme_bw() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = "bottom")

final_plot <- heatmap_plot + rmse_plot + plot_layout(widths = c(5, 2))

pdf(file = "02b_DistributedFF_TransportDiffusion.pdf", width = 15, height = 8)
print(final_plot)
dev.off()
