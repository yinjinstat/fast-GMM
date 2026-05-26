#========================
# Utility: Misclassification rate (handles label switching)
#========================
err_rate <- function(pred, truth) {
  e1 <- mean(pred != truth)
  e2 <- mean((3 - pred) != truth)
  min(e1, e2)
}

#========================
# Main Experiment Flow
#========================
library(mclust)

#--- 1. Parameters ---
n <- 1000
p <- 20
num_repetitions <- 500
a <- 2
b <- 1
class_probs <- c(0.5, 0.5)

#--- 2. Data Generation ---
generate_mixture_data <- function(n, p, probs) {
  if (length(probs) != 2 || abs(sum(probs) - 1) > .Machine$double.eps^0.5) {
    stop("Probability parameter 'probs' must be a vector of length 2 that sums to 1.")
  }
  true_labels <- sample(c(1, 2), size = n, replace = TRUE, prob = probs)
  x <- matrix(0, nrow = n, ncol = p)
  n1 <- sum(true_labels == 1)
  n2 <- n - n1
  x[true_labels == 1, 1] <- rnorm(n1, mean = a, sd = b)
  x[true_labels == 2, 1] <- rnorm(n2, mean = -a, sd = b)
  if (p > 1) {
    x[, 2:p] <- rnorm(n * (p - 1), mean = 0, sd = 4)
  }
  list(data = x, labels = true_labels)
}

#--- 3. Theoretical Bayes Error Rate ---
pi1 <- class_probs[1]; pi2 <- class_probs[2]
mu1_scalar <- a; mu2_scalar <- -a; sigma_scalar <- b
decision_boundary <- (mu1_scalar + mu2_scalar) / 2 -
  (sigma_scalar^2 / (mu1_scalar - mu2_scalar)) * log(pi1 / pi2)
cat("Decision boundary (x_1):", decision_boundary, "\n")
error_type1 <- pnorm(decision_boundary, mean = mu1_scalar, sd = sigma_scalar)
error_type2 <- pnorm(decision_boundary, mean = mu2_scalar, sd = sigma_scalar, lower.tail = FALSE)
s_bayes_rate <- error_type1 * pi1 + error_type2 * pi2
cat("Theoretical Bayes error rate (s):", s_bayes_rate, "\n")

#--- 4. Main Loop ---
total_times <- numeric(num_repetitions)
total_starts <- numeric(num_repetitions)
threshold <- 0.2 * p / sqrt(n)
cat("Threshold:", threshold, "\n")
cat("Starting", num_repetitions, "repeated experiments...\n")

for (i in 1:num_repetitions) {
  start_time <- proc.time()
  random_starts_count <- 0
  
  sim <- generate_mixture_data(n, p, class_probs)
  X_data <- sim$data
  true_labels <- sim$labels
  
  repeat {
    random_starts_count <- random_starts_count + 1
    
    # --- K-means random initialization ---
    km_init <- tryCatch({
      kmeans(X_data, centers = 2, algorithm = "Lloyd", nstart = 1)
    }, error = function(e) {
      cat("K-means failed on start", random_starts_count, ", retrying.\n")
      return(NULL)
    })
    if (is.null(km_init)) next
    
    km_err <- err_rate(km_init$cluster, true_labels)
    
    # --- Convert k-means clustering results to an n x 2 z-matrix ---
    z_init <- matrix(0, nrow = n, ncol = 2)
    for (j in 1:n) {
      z_init[j, km_init$cluster[j]] <- 1
    }
    
    # --- Run EM using me() ---
    em_fit <- tryCatch({
      me(data = as.matrix(X_data),
         modelName = "EEE",
         z = z_init,
         control = emControl(tol = 1e-10, itmax = 200000000))
    }, error = function(e) {
      cat("me() failed on start", random_starts_count, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    if (is.null(em_fit)) next
    if (is.null(em_fit$z) || any(is.na(em_fit$z))) {
      cat("  Attempt", random_starts_count, ": EM result abnormal, retrying.\n")
      next
    }
    
    em_class <- map(em_fit$z)
    misclassification_rate <- err_rate(em_class, true_labels)
    
    cat(sprintf("  Attempt %d: km_err=%.4f, em_err=%.4f\n",
                random_starts_count, km_err, misclassification_rate))
    
    # --- Stopping criterion: Difference between misclassification rate and Bayes error rate < threshold ---
    if ((misclassification_rate - s_bayes_rate) < threshold) break
  }
  
  end_time <- proc.time()
  elapsed_time <- (end_time - start_time)["elapsed"]
  total_times[i] <- elapsed_time
  total_starts[i] <- random_starts_count
  
  cat(sprintf("Experiment %d/%d: Elapsed time %.2f seconds, tried %d random initializations, misclassification rate = %.4f\n",
              i, num_repetitions, elapsed_time, random_starts_count, misclassification_rate))
}

cat("\n--- Experiment completed ---\n")
cat("Total repeated experiments:", num_repetitions, "\n")
cat("Average time required:", mean(total_times), "seconds\n")
cat("Average number of random starts required:", mean(total_starts), "\n")

p_old <- par(no.readonly = TRUE); on.exit(par(p_old), add = TRUE)