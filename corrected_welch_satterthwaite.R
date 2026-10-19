cws2 <- function(s2_1, s2_2, nu_1, nu_2, w1 = 1, w2 = 1) {
  stopifnot(
    length(s2_1) == 1L, length(s2_2) == 1L,
    length(nu_1) == 1L, length(nu_2) == 1L,
    length(w1) == 1L, length(w2) == 1L,
    is.numeric(s2_1), is.numeric(s2_2),
    is.numeric(nu_1), is.numeric(nu_2),
    is.numeric(w1), is.numeric(w2),
    s2_1 >= 0, s2_2 >= 0,
    nu_1 > 0, nu_2 > 0,
    w1 > 0, w2 > 0
  )

  s2 <- w1 * s2_1 + w2 * s2_2
  den <- w1^2 * s2_1^2 / (nu_1 + 2) +
    w2^2 * s2_2^2 / (nu_2 + 2)

  if (den <= 0) return(NA_real_)
  s2^2 / den - 2
}

cws <- function(s2, nu, w = rep(1, length(s2))) {
  n <- length(s2)

  stopifnot(
    n >= 1L,
    length(nu) == n,
    length(w) == n,
    is.numeric(s2), is.numeric(nu), is.numeric(w),
    all(s2 >= 0), all(nu > 0), all(w > 0)
  )

  den <- sum(w^2 * s2^2 / (nu + 2))
  if (den <= 0) return(NA_real_)
  sum(w * s2)^2 / den - 2
}
