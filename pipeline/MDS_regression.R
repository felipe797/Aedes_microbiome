### Function for performing Multidimensional Scaling (MDS) and fitting a linear model
### on theoretical distances for a tridimensional plot


perform_mds_lm <- function(abu) {
  # Perform Multidimensional Scaling (MDS) using vegan package
  mds <- vegan::metaMDS(t(abu), k = 3)
  
  # Extract MDS coordinates
  x <- mds$species[, 1]
  y <- mds$species[, 2]
  z <- mds$species[, 3]
  
  # Fit linear model using base R functions
  fit <- stats::lm(z ~ x + y)
  
  # Predict fitted points
  fitpoints <- stats::predict(fit)
  
  # Prepare confidence intervals
  ci <- list(
    z = matrix(nrow = length(x), data = rep(0.1, 2 * length(x)))  # Set confidence interval to 0.1
  )
  
  # Return a list containing MDS coordinates, fitted points, and confidence intervals
  result <- list(
    mds_coordinates = mds$species,
    fitted_points = fitpoints,
    confidence_intervals = ci
  )
  
  return(result)
}
