# 02a_DistributedFF_Cyclone.R
# Distributed Frozen Field cyclone simulation (Figure 1)

library(MASS)
library(circular)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ReacTran)

# Rotation around center C
rotate_point <- function(s1, s2, theta, C = c(0, 0)) {
  s1_rot <- (s1 - C[1]) * cos(theta) - (s2 - C[2]) * sin(theta)
  s2_rot <- (s1 - C[1]) * sin(theta) + (s2 - C[2]) * cos(theta)
  return(c(s1_rot + C[1], s2_rot + C[2]))
}

# Grid
Lx <- 100; Ly <- 100
Nx <- 100; Ny <- 100
n <- 100

grid.x_1D <- setup.grid.1D(x.up = 0, L = Lx, N = Nx)
grid.y_1D <- setup.grid.1D(x.up = 0, L = Ly, N = Ny)
x_coords <- grid.x_1D$x.mid
y_coords <- grid.y_1D$x.mid

# Initial condition X_S(s)
initial_condition <- function(s1, s2) {
  s1_0 <- 30; s2_0 <- 40
  sigma1 <- 5; sigma2 <- 5
  exp(-((s1 - s1_0)^2 / (2 * sigma1^2) + (s2 - s2_0)^2 / (2 * sigma2^2)))
}

W <- outer(x_coords, y_coords, FUN = initial_condition)

# Parameters
v <- c(5, 0)           # base velocity
Tmax <- 7              # max time
C <- c(n/2, n/2 + 15)  # center of rotation
n_theta <- 10          # number of rotation angles
mu_theta <- 0          # mean angle
sd_theta <- 5          # sd for wrapped normal

seeds_theta <- c(45453, 46346, 6546, 675685, 49834589, 64564, 2443545, 45363)
seeds_v <- c(2374, 5345, 4567456, 464554, 93539, 353534, 6653, 5745)

plot_titles <- list(
  expression(Z(bold(s), t == 0)), expression(Z(bold(s), t == 1)),
  expression(Z(bold(s), t == 2)), expression(Z(bold(s), t == 3)),
  expression(Z(bold(s), t == 4)), expression(Z(bold(s), t == 5)),
  expression(Z(bold(s), t == 6)), expression(Z(bold(s), t == 7))
)

plot_list <- list()

for (t in 0:Tmax) {
  idx <- t + 1
  
  # Sample rotation angles from wrapped normal
  set.seed(seeds_theta[idx])
  thetas <- rwrappednormal(n_theta, mu = circular(mu_theta), sd = sd_theta,
                           control.circular = list(units = "radians"))
  p_theta <- dwrappednormal(thetas, mu = circular(mu_theta), sd = sd_theta)
  p_theta <- p_theta / sum(p_theta)
  
  # Sample velocity
  set.seed(seeds_v[idx])
  v_t <- mvrnorm(n = 1, mu = v, Sigma = diag(2))
  
  X_big <- matrix(0, n, n)
  
  if (t == 0) {
    X_big <- W
  } else {
    # Distributed FF: weighted sum over rotated/translated fields
    for (k in 1:n_theta) {
      X_rot <- matrix(0, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          s_rot <- rotate_point(i - v_t[1] * t, j - v_t[2] * t, thetas[k], C = C)
          src_i <- round(s_rot[1])
          src_j <- round(s_rot[2])
          if (src_i > 0 && src_i <= n && src_j > 0 && src_j <= n) {
            X_rot[i, j] <- W[src_i, src_j]
          }
        }
      }
      X_big <- X_big + p_theta[k] * X_rot
    }
  }
  
  # Track center position
  center_x <- C[1] + v_t[1] * t
  center_y <- C[2] + v_t[2] * t
  
  plot_list[[idx]] <- ggplot(melt(X_big), aes(x = Var1, y = Var2)) +
    geom_raster(aes(fill = value)) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    labs(x = expression(s[1]), y = expression(s[2]), title = plot_titles[[idx]]) +
    annotate("point", x = center_x, y = center_y, color = "black", 
             size = 3, shape = 4, stroke = 1.5) +
    theme_bw() +
    theme(legend.position = "none")
}

pdf(file = "02a_DistributedFF_Cyclone.pdf", width = 35, height = 7)
grid.arrange(grobs = lapply(plot_list, function(p) {
  p + theme(axis.text = element_text(size = 20),
            axis.title = element_text(size = 22),
            plot.title = element_text(size = 24))
}), nrow = 1)
dev.off()
