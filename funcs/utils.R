library(lhs)

## Locs gen funcs ----------------------------------------
#' Generate `n` random locations in `d` dimensional spatial domain following the 
#' Latin hypercube design.
#' @param n number of locations, also the dimension of the MVN probability
#' @param d dimension of the spatial domain, where the Matern kernel is defined
#' @return an n X d matrix
latin_gen <- function(n, d) {
  as.matrix(lhs::randomLHS(n, d))
}
#' Generate `n` locations in `d` dimensional spatial domain following a grid design
#' in the unit hypercube
#' @param n number of locations, also the dimension of the MVN probability
#' @param d dimension of the spatial domain, where the Matern kernel is defined
#' @param unitDist distance between grid bands
#' @return a list containing:
#' \\item{grid}{the generated spatial locations, an n X d matrix}
#' \\item{unitDist}{the distance between grid bands}
grid_gen <- function(n, d, unitDist = NULL) {
  m <- floor(n^{
    1 / d
  }) + 1
  if (is.null(unitDist)) {
    unitDist <- 1 / m
  }
  edges <- list()
  for (i in 1:d) {
    edges[[i]] <- (1:m) * unitDist
  }
  grid <- as.matrix(expand.grid(edges))
  return(list(grid = grid[1:n, , drop = F], unit_dist = unitDist))
}
