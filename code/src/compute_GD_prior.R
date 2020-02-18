compute_GD_prior <- function(N_stick, alpha0, k){
  xi <- numeric(N_stick)
  for (j in seq(1, N_stick-1)) {
    shape1 = alpha0/N_stick + sum(k == j)
    shape2 = (N_stick - j)/N_stick * alpha0
    for (s in seq(j + 1, N_stick)) {
      shape2 = shape2 + sum(k == s)
    }
    xi[j] <- rbeta(1, shape1 = shape1, shape2 = shape2)
  }
  return(xi)
}
