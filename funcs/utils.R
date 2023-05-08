library(lhs)

## Locs gen funcs ----------------------------------------
latin_gen <- function(n, d) {
  as.matrix(lhs::randomLHS(n, d))
}
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
