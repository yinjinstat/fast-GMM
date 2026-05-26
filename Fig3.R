set.seed(123)
par(mar = c(4.5, 4.5, 4, 2))

# Load necessary packages
library(MASS)

# Basic parameters
mu <- c(4, 2) # Separation direction (two cluster centers are at ±mu)
n_per_cluster <- 100
pi1 <- 0.5; pi2 <- 0.5

# Direction basis: u is the unit vector along mu; v is the unit vector orthogonal to u (rotated 90 degrees counter-clockwise)
mu_norm <- sqrt(sum(mu^2))
u <- mu / mu_norm
v <- c(-u[2], u[1]) # Orthogonal to u, length is 1

# Single cluster covariance S = a uu^T + b vv^T (larger in v direction), b = 4 as requested
a <- 0.5
b <- 16.0
S <- a * tcrossprod(u, u) + b * tcrossprod(v, v)

# Push the two clusters to the top-left and bottom-right (diagonal positions)
# Adjust parameters so the two clusters are in diagonal positions
delta_v <- 0  # No overall upward shift needed
sep_scale <- 2.5  # Increase the separation scale to move the clusters closer to the diagonal

m1 <- -sep_scale * mu  # Bottom-left direction
m2 <- sep_scale * mu   # Top-right direction

# Generate mixed data (covariance is S)
X1 <- mvrnorm(n = n_per_cluster, mu = m1, Sigma = S)
X2 <- mvrnorm(n = n_per_cluster, mu = m2, Sigma = S)
X_mix <- rbind(X1, X2)
X_mix <- scale(X_mix, scale = FALSE)

# Calculate actual cluster centers (after centering)
cluster_mean1 <- colMeans(X1) - colMeans(rbind(X1, X2))
cluster_mean2 <- colMeans(X2) - colMeans(rbind(X1, X2))

# Overall mean and overall covariance
m_mix <- pi1 * m1 + pi2 * m2
Sigma_mix <- S + tcrossprod(sep_scale * mu, sep_scale * mu)

# Generate matching single Gaussian data
X_gauss <- mvrnorm(n = 100, mu = m_mix, Sigma = Sigma_mix)
X_gauss <- scale(X_gauss, scale = FALSE)

# Calculate coordinate range to ensure all points in X_mix are within the bounding box
x_range <- range(X_mix[,1])
y_range <- range(X_mix[,2])
x_pad <- 0.1 * (x_range[2] - x_range[1])
y_pad <- 0.1 * (y_range[2] - y_range[1])
xlim <- c(x_range[1] - x_pad, x_range[2] + x_pad)
ylim <- c(y_range[1] - y_pad, y_range[2] + y_pad)

# Filter out points in X_gauss that fall outside the range
keep_idx <- (X_gauss[,1] >= xlim[1] & X_gauss[,1] <= xlim[2] & 
               X_gauss[,2] >= ylim[1] & X_gauss[,2] <= ylim[2])
X_gauss <- X_gauss[keep_idx, ]

# Utility function: Draw a line passing through point p=(px,py) with direction vector d, intersecting the bounding box (xlim, ylim)
draw_infinite_line <- function(px, py, dx, dy, xlim, ylim, col="black", lwd=2, lty=1) {
  # Parameterization: x = px + t*dx, y = py + t*dy
  ts <- c()
  if (abs(dx) > 1e-12) {
    ts <- c(ts, (xlim[1] - px)/dx, (xlim[2] - px)/dx)
  }
  if (abs(dy) > 1e-12) {
    ts <- c(ts, (ylim[1] - py)/dy, (ylim[2] - py)/dy)
  }
  
  # Calculate all candidate intersection points
  cand <- cbind(px + ts * dx, py + ts * dy)
  
  # Select points that fall inside or on the boundaries of the bounding box
  eps <- 1e-6
  valid <- (cand[,1] >= xlim[1]-eps & cand[,1] <= xlim[2]+eps &
              cand[,2] >= ylim[1]-eps & cand[,2] <= ylim[2]+eps)
  
  if (sum(valid) >= 2) {
    valid_points <- cand[valid, , drop = FALSE]
    # Select the two furthest points
    if (nrow(valid_points) > 2) {
      dists <- apply(valid_points, 1, function(p) (p[1]-px)^2 + (p[2]-py)^2)
      idx <- order(dists, decreasing = TRUE)[1:2]
      valid_points <- valid_points[idx, ]
    }
    segments(valid_points[1,1], valid_points[1,2], 
             valid_points[2,1], valid_points[2,2], 
             col=col, lwd=lwd, lty=lty)
  }
}

# Plotting
plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "X1", ylab = "X2",
     asp = 1, cex.axis = 1.2, cex.lab = 1.5,
     main = "")

# Hollow points: matching Gaussian (grey-blue)
points(X_gauss[,1], X_gauss[,2], pch = 1, col = "#445C8A")

# Solid points: mixed clusters (uniform black)
points(X_mix[,1], X_mix[,2], pch = 16, col = "black")

# Draw a line connecting the two cluster centers (solid black)
# Line direction vector
line_dir <- cluster_mean2 - cluster_mean1
line_dir_norm <- line_dir / sqrt(sum(line_dir^2))

# Draw the line from cluster_mean1
draw_infinite_line(cluster_mean1[1], cluster_mean1[2], 
                   line_dir_norm[1], line_dir_norm[2], 
                   xlim, ylim, col="black", lwd=2, lty=1)

# Calculate the midpoint of the line
midpoint <- (cluster_mean1 + cluster_mean2) / 2

# Perpendicular direction (rotated 90 degrees counter-clockwise)
perp_dir <- c(-line_dir_norm[2], line_dir_norm[1])

# Draw a perpendicular line from the midpoint (dashed dark grey)
draw_infinite_line(midpoint[1], midpoint[2], 
                   perp_dir[1], perp_dir[2], 
                   xlim, ylim, col="gray30", lwd=2, lty=2)

