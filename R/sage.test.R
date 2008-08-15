## A vectorized version of the modified chisquare statistic in package
## sagenhaft. 'x' and 'y' are parallel entries in the contingency table,
## 'n1' and 'n2' are the respective outer sums:
##      |   |   |
##    --|-------|--
##      | x | y | 
##      |...|...|
##    --|---|---|--
##      |n1 | n2|

sage.test <- function (x, y, n1 = sum(x), n2 = sum(y)){
  if (any(is.na(x)) || any(is.na(y)))
    stop("missing values not allowed")
  x <- round(x)
  y <- round(y)
  if (any(x < 0) || any(y < 0))
    stop("x and y must be non-negative")
  if (length(x) != length(y))
    stop("x and y must have same length")
  n1 <- round(n1)
  n2 <- round(n2)
  if(length(n1)==1) n1 <- rep(n1, length(x))
  if(length(n2)==1) n2 <- rep(n2, length(x))
  if (length(n1) != length(n2) | length(n1) != length(x))
    stop("n1 and n2 must have same length as x and y")
  if (!missing(n1) && any(x > n1))
    stop("x cannot be greater than n1")
  if (!missing(n2) && any(y > n2))
    stop("y cannot be greater than n2")
  size <- x + y
  p.value <- rep(1, length(x))
  if (any(n1 == n2)) {
    i <- size > 0 & n1 == n2
    if (any(i)) {
      xI <- pmin(x[i], y[i])
      sizeI <- size[i]
      p.value[i] <- pbinom(xI, size = sizeI, prob = 0.5) +
        pbinom(sizeI - xI + 0.5, size = sizeI, prob = 0.5,
               lower.tail = FALSE)
    }
  }
  if (any(n1 != n2)) {
    prob <- n1/(n1 + n2)
    if (any(size > 10000 & n1 != n2)) {
      big <- size > 10000 & n1 != n2
      ibig <- (1:length(x))[big]
      for (i in ibig)
        p.value[i] <- chisq.test(matrix(c(x[i], y[i], n1[i] - x[i], n2[i] - y[i]), 2, 2))$p.value
    }
    size0 <- size[size > 0 & size <= 10000 & n1 != n2]
    prob0 <- prob[size > 0 & size <= 10000 & n1 != n2]
    mar0 <- unique(cbind(size0, prob0), MARGIN=1)
    if (nrow(mar0))
      for (ind in 1:nrow(mar0)) {
        isize <- mar0[ind,1]
        iprob <- mar0[ind,2]
        i <- (size == isize) & (prob==iprob) & n1 != n2
        p <- dbinom(0:isize, p = (prob[i])[1], size = isize)
        o <- order(p)
        cumsump <- cumsum(p[o])[order(o)]
        p.value[i] <- cumsump[x[i] + 1]
      }
  }
  p.value[p.value>1] <- 1
  return(p.value)
}

