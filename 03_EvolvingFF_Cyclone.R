# 03_EvolvingFF_Cyclone.R
# Evolving Frozen Field cyclone simulation (Figure 2)

library(MASS)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ReacTran)

# Spiral velocity field v(s,t) = (R(s)*cos(t)/(t+1), R(s)*sin(t)/(t+1))
# where R(s) = ||s||/2
spiral_displacement <- function(s1, s2, t, a = 1, b = 1) {
  R <- sqrt(s1^2 + s2^2) / 2
  v1 <- R * cos(a * t) / (b * (t + 1))
  v2 <- R * sin(a * t) / (b * (t + 1))
  return(cbind(v1, v2))
}

every_n <- function(x, by = 2) {
  x <- sort(x)
  x[seq(1, length(x), by = by)]
}

# Grid
Lx <- 100; Ly <- 100
Nx <- 100; Ny <- 100
n <- 100

grid.x_1D <- setup.grid.1D(x.up = 0, L = Lx, N = Nx)
grid.y_1D <- setup.grid.1D(x.up = 0, L = Ly, N = Ny)
x_coords <- grid.x_1D$x.mid
y_coords <- grid.y_1D$x.mid

# Initial condition
initial_condition <- function(s1, s2) {
  s1_0 <- 30; s2_0 <- 40
  sigma1 <- 5; sigma2 <- 5
  exp(-((s1 - s1_0)^2 / (2 * sigma1^2) + (s2 - s2_0)^2 / (2 * sigma2^2)))
}

W <- outer(x_coords, y_coords, FUN = initial_condition)

# Time
Tmax <- 9
time_vec <- 0:(Tmax - 1)

# Compute velocity field v(s,t)
V_s1 <- V_s2 <- list()
V_s1[[1]] <- V_s2[[1]] <- matrix(0, n, n)

for (t in 2:Tmax) {
  vel_s1 <- vel_s2 <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      pos_prev <- spiral_displacement(s1 = i, s2 = j, t = time_vec[t - 1])
      pos_curr <- spiral_displacement(s1 = i, s2 = j, t = time_vec[t])
      vel_s1[i, j] <- pos_curr[, 1] - pos_prev[, 1]
      vel_s2[i, j] <- pos_curr[, 2] - pos_prev[, 2]
    }
  }
  V_s1[[t]] <- vel_s1
  V_s2[[t]] <- vel_s2
}

# Velocity field plots
plot_vel <- list()
for (t in 1:Tmax) {
  Vs_df <- data.frame(
    s1 = as.numeric(melt(V_s1[[t]])$Var1),
    s2 = as.numeric(melt(V_s1[[t]])$Var2),
    v_s1 = melt(V_s1[[t]])$value,
    v_s2 = melt(V_s2[[t]])$value
  )
  
  keep_s1 <- every_n(unique(Vs_df$s1), by = 10)
  keep_s2 <- every_n(unique(Vs_df$s2), by = 10)
  Vs_sub <- filter(Vs_df, s1 %in% keep_s1 & s2 %in% keep_s2)
  
  plot_vel[[t]] <- ggplot(Vs_sub, aes(x = s1, y = s2)) +
    geom_segment(aes(xend = s1 + v_s1/2, yend = s2 + v_s2/2),
                 arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.5) +
    labs(title = bquote(bold(v)(bold(s), t == .(t - 1))),
         x = expression(v[1](bold(s), t)), y = expression(v[2](bold(s), t))) +
    theme_bw() +
    theme(legend.position = "none")
}

# Evolving FF: Z(s,t) = X_S(s - v(s,t)*t)
Z_fields <- list()
for (t in time_vec) {
  t_idx <- t + 1
  Z_fields[[t_idx]] <- matrix(NA, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      src_i <- round(i - V_s1[[t_idx]][i, j] * t)
      src_j <- round(j - V_s2[[t_idx]][i, j] * t)
      
      if (src_i > 0 && src_i <= n && src_j > 0 && src_j <= n) {
        Z_fields[[t_idx]][i, j] <- W[src_i, src_j]
      } else {
        Z_fields[[t_idx]][i, j] <- 0
      }
    }
  }
}

# Field plots
plot_titles <- list(
  expression(Z(bold(s), t == 0)), expression(Z(bold(s), t == 1)),
  expression(Z(bold(s), t == 2)), expression(Z(bold(s), t == 3)),
  expression(Z(bold(s), t == 4)), expression(Z(bold(s), t == 5)),
  expression(Z(bold(s), t == 6)), expression(Z(bold(s), t == 7)),
  expression(Z(bold(s), t == 8))
)

plot_field <- list()
for (t in 1:Tmax) {
  plot_field[[t]] <- ggplot(melt(Z_fields[[t]]), aes(x = Var1, y = Var2)) +
    geom_raster(aes(fill = value)) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    labs(x = expression(s[1]), y = expression(s[2]), title = plot_titles[[t]]) +
    theme_bw() +
    theme(legend.position = "none")
}

pdf(file = "03_EvolvingFF_Cyclone.pdf", width = 35, height = 14)
grid.arrange(
  grobs = lapply(c(plot_field[1:8], plot_vel[1:8]), function(p) {
    p + theme(axis.text = element_text(size = 20),
              axis.title = element_text(size = 22),
              plot.title = element_text(size = 24))
  }),
  nrow = 2
)
dev.off()
