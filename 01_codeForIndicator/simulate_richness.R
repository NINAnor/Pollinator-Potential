simulate_richness <- function(
    n = 50,
    x_min = 1,
    x_max = 100,
    z = 1.96,
    round_counts = TRUE,
    seed = NULL
) {

  if (x_min <= 0 || x_max <= x_min) {
    stop("Require 0 < x_min < x_max.")
  }

  if (n <= 0 || z <= 0) {
    stop("'n' and 'z' must be positive.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  mu <- (log(x_max) + log(x_min)) / 2
  sigma <- (log(x_max) - log(x_min)) / (2 * z)

  richness <- rlnorm(
    n = n,
    meanlog = mu,
    sdlog = sigma
  )

  if (round_counts) {
    richness <- round(richness)
  }

  richness
}