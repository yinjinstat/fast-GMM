# --- 0. Prepare Environment ---
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
if (!requireNamespace("ellipse", quietly = TRUE)) install.packages("ellipse")
library(MASS)
library(ellipse)

# --- Utility Functions ---
ellipse_95 <- function(Sigma, mu = c(0, 0), level = 0.95, npoints = 200) {
  ellipse(Sigma, centre = mu, level = level, npoints = npoints)
}

# --- 1. Global Settings ---
# mfrow = c(1, 2) displays two plots side-by-side
# pty = "s" ensures the plotting region is perfectly square
par(mfrow = c(1, 2), pty = "s", mar = c(2.5, 2.5, 1, 1))
set.seed(456)

n_samples <- 200
mu <- c(0, 0)
n1 <- n_samples %/% 2
n2 <- n_samples - n1

# --- Uniform Style Settings (applied to both plots) ---
pch_c1 <- 20   # cdot (solid small circle)
pch_c2 <- 4    # x (cross)
cex_val <- 0.6 # slightly smaller
col_c1 <- "black"
col_c2 <- "black"
lty_comp <- 2

# =========================================================
# Plot 1: X ~ 0.5 N(0, I2) + 0.5 N(0, 4^2 I2)
# =========================================================
Sigma1_comp1 <- diag(2) * 1        
Sigma1_comp2 <- diag(2) * (4^2)    

data1_comp1 <- mvrnorm(n = n1, mu = mu, Sigma = Sigma1_comp1)
data1_comp2 <- mvrnorm(n = n2, mu = mu, Sigma = Sigma1_comp2)
data1_mix <- rbind(data1_comp1, data1_comp2)

ell1_c1 <- ellipse_95(Sigma1_comp1, mu = mu)
ell1_c2 <- ellipse_95(Sigma1_comp2, mu = mu)

Sigma1_overall_hat <- cov(data1_mix)
ell1_overall <- ellipse_95(Sigma1_overall_hat, mu = colMeans(data1_mix))

# Plot 1 is an isotropic mixture; keeping X and Y ranges consistent looks better
all_coords1 <- rbind(data1_mix, ell1_c1, ell1_c2, ell1_overall)
lim1 <- max(abs(all_coords1)) * 1.1
range1 <- c(-lim1, lim1)

plot(0, type = "n",
     xlim = range1, ylim = range1,
     xlab = "", ylab = "", main = "") # No title

points(data1_comp1, pch = pch_c1, col = col_c1, cex = cex_val)
points(data1_comp2, pch = pch_c2, col = col_c2, cex = cex_val)

lines(ell1_c1, lty = lty_comp, lwd = 1.2, col = col_c1)
lines(ell1_c2, lty = lty_comp, lwd = 1.2, col = col_c2)

# =========================================================
# Plot 2: X ~ 0.5 N(0, I2) + 0.5 N(0, diag(1, 4^2))
# =========================================================
# Note: Here, the first element of Sigma2_comp2 is 16 (large variance on X-axis), and the fourth element is 1 (small variance on Y-axis)
Sigma2_comp1 <- matrix(c(1, 0, 0, 1), nrow = 2)   
Sigma2_comp2 <- matrix(c(4^2, 0, 0, 1), nrow = 2) 

data2_comp1 <- mvrnorm(n = n1, mu = mu, Sigma = Sigma2_comp1)
data2_comp2 <- mvrnorm(n = n2, mu = mu, Sigma = Sigma2_comp2)
data2_mix <- rbind(data2_comp1, data2_comp2)

ell2_c1 <- ellipse_95(Sigma2_comp1, mu = mu)
ell2_c2 <- ellipse_95(Sigma2_comp2, mu = mu)

Sigma2_overall_hat <- cov(data2_mix)
ell2_overall <- ellipse_95(Sigma2_overall_hat, mu = colMeans(data2_mix))

# --- Modified Range Calculation Logic for Plot 2 ---
# Calculate maximum ranges for X and Y separately instead of forcing them to be identical
all_coords2 <- rbind(data2_mix, ell2_c1, ell2_c2, ell2_overall)

max_x <- max(abs(all_coords2[, 1])) * 1.1
max_y <- max(abs(all_coords2[, 2])) * 1.1

# Here, ylim is set based on Y data only and will be much smaller than the X range
# However, because of par(pty="s"), R will stretch this smaller Y range to fill the square, solving the "flattened" issue
xlim2 <- c(-max_x, max_x)
ylim2 <- c(-max_y, max_y)

plot(0, type = "n",
     xlim = xlim2, ylim = ylim2,
     xlab = "", ylab = "", main = "") # No title

points(data2_comp1, pch = pch_c1, col = col_c1, cex = cex_val)
points(data2_comp2, pch = pch_c2, col = col_c2, cex = cex_val)

lines(ell2_c1, lty = lty_comp, lwd = 1.2, col = col_c1)
lines(ell2_c2, lty = lty_comp, lwd = 1.2, col = col_c2)