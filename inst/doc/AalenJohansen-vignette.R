## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(AalenJohansen)

set.seed(2)

jump_rate <- function(i, t, u){
  if(i == 1){
    2 / (1+1/2*t)
  } else if(i == 2){
    3 / (1+1/2*t)
  } else{
    0
  }
}

mark_dist <- function(i, s, v){
  if(i == 1){
    c(0, 1/2, 1/2)
  } else if(i == 2){
    c(2/3, 0, 1/3)
  } else{
    0
  } 
}

lambda <- function(t){
  A <- matrix(c(2/(1+1/2*t)*mark_dist(1, t, 0), 3/(1+1/2*t)*mark_dist(2, t, 0), rep(0, 3)),
              nrow = 3, ncol = 3, byrow = TRUE)
  diag(A) <- -rowSums(A)
  A
}

## -----------------------------------------------------------------------------
n <- 1000

c <- runif(n, 0, 5)

sim <- list()
for(i in 1:n){
  sim[[i]] <- sim_path(sample(1:2, 1), rates = jump_rate, dists = mark_dist,
                       tn = c[i], bs = c(2*c[i], 3*c[i], 0))
}

## -----------------------------------------------------------------------------
sum(c == unlist(lapply(sim, FUN = function(z){tail(z$times, 1)}))) / n

## -----------------------------------------------------------------------------
fit <- aalen_johansen(sim)

## ---- fig.align = 'center', fig.height = 3, fig.width = 6---------------------
v1 <- unlist(lapply(fit$Lambda, FUN = function(L) L[2,1]))
v0 <- fit$t
p <- unlist(lapply(fit$p, FUN = function(L) L[2]))
P <- unlist(lapply(prodint(0, 5, 0.01, lambda), FUN = function(L) (c(1/2, 1/2, 0) %*% L)[2]))

par(mfrow = c(1, 2))
par(mar = c(2.5, 2.5, 1.5, 1.5))

plot(v0, v1, type = "l", lty = 2, xlab = "", ylab = "", main = "Hazard")
lines(v0, 4*log(1+1/2*v0))
plot(v0, p, type = "l", lty = 2, xlab = "", ylab = "", main = "Probability")
lines(seq(0, 5, 0.01), P)

## -----------------------------------------------------------------------------
jump_rate <- function(i, t, u){
  if(i == 1){
    2
  } else if(i == 2){
    3
  } else{
    0
  }
}

mark_dist <- function(i, s, v){
  if(i == 1){
    c(0, 1/2, 1/2)
  } else if(i == 2){
    c(2/3, 0, 1/3)
  } else{
    0
  }
}

lambda <- function(t, x){
  A <- matrix(c(2/(1+x*t)*mark_dist(1, t, 0), 3/(1+x*t)*mark_dist(2, t, 0), rep(0, 3)),
              nrow = 3, ncol = 3, byrow = TRUE)
  diag(A) <- -rowSums(A)
  A
}

n <- 10000

X <- runif(n)
c <- runif(n, 0, 5)

sim <- list()
for(i in 1:n){
  rates <- function(j, y, z){jump_rate(j, y, z)/(1+X[i]*y)}
  sim[[i]] <- sim_path(sample(1:2, 1), rates = rates, dists = mark_dist,
                       tn = c[i], bs = c(2*c[i], 3*c[i], 0))
  sim[[i]]$X <- X[i]
}

## -----------------------------------------------------------------------------
sum(c == unlist(lapply(sim, FUN = function(z){tail(z$times, 1)}))) / n

## -----------------------------------------------------------------------------
x1 <- 0.2
x2 <- 0.8

fit1 <- aalen_johansen(sim, x = x1)
fit2 <- aalen_johansen(sim, x = x2)

## ---- fig.align = 'center', fig.height = 3, fig.width = 6---------------------
v11 <- unlist(lapply(fit1$Lambda, FUN = function(L) L[2,1]))
v10 <- fit1$t
v21 <- unlist(lapply(fit2$Lambda, FUN = function(L) L[2,1]))
v20 <- fit2$t
p1 <- unlist(lapply(fit1$p, FUN = function(L) L[2]))
P1 <- unlist(lapply(prodint(0, 5, 0.01, function(t){lambda(t, x = x1)}),
FUN = function(L) (c(1/2, 1/2, 0) %*% L)[2]))
p2 <- unlist(lapply(fit2$p, FUN = function(L) L[2]))
P2 <- unlist(lapply(prodint(0, 5, 0.01, function(t){lambda(t, x = x2)}),
FUN = function(L) (c(1/2, 1/2, 0) %*% L)[2]))

par(mfrow = c(1, 2))
par(mar = c(2.5, 2.5, 1.5, 1.5))

plot(v10, v11, type = "l", lty = 2, xlab = "", ylab = "", main = "Hazard", col = "red")
lines(v10, 2/x1*log(1+x1*v10), col = "red")
lines(v20, v21, lty = 2, col = "blue")
lines(v20, 2/x2*log(1+x2*v20), col = "blue")

plot(v10, p1, type = "l", lty = 2, xlab = "", ylab = "", main = "Probability", col = "red")
lines(seq(0, 5, 0.01), P1, col = "red")
lines(v20, p2, lty = 2, col = "blue")
lines(seq(0, 5, 0.01), P2, col = "blue")

## -----------------------------------------------------------------------------
jump_rate_enlarged <- function(i, t, u){
  if(i == 1){
    2.5
  } else if(i == 2){
    4
  } else{
    0
  }
}

mark_dist_enlarged <- function(i, s, v){
  if(i == 1){
    c(0, 2/5, 2/5, 1/5)
  } else if(i == 2){
    c(2/4, 0, 1/4, 1/4)
  } else{
    0
  }
}

n <- 10000

X <- runif(n)

tn <- 5
sim <- list()
for(i in 1:n){
  rates <- function(j, y, z){jump_rate_enlarged(j, y, z)/(1+X[i]*y)}
  sim[[i]] <- sim_path(sample(1:2, 1), rates = rates, dists = mark_dist_enlarged,
                       tn = tn, abs = c(FALSE, FALSE, TRUE, TRUE),
                       bs = c(2.5*tn, 4*tn, 0, 0))
  sim[[i]]$X <- X[i]
}


## -----------------------------------------------------------------------------
sum(tn == unlist(lapply(sim, FUN = function(z){tail(z$times, 1)}))
    | 4 == unlist(lapply(sim, FUN = function(z){tail(z$states, 1)}))) / n

## -----------------------------------------------------------------------------
fit1 <- aalen_johansen(sim, x = x1, collapse = TRUE)
fit2 <- aalen_johansen(sim, x = x2, collapse = TRUE)

## ---- fig.align = 'center', fig.height = 3, fig.width = 6---------------------
v11 <- unlist(lapply(fit1$Lambda, FUN = function(L) L[2,1]))
v10 <- fit1$t
v21 <- unlist(lapply(fit2$Lambda, FUN = function(L) L[2,1]))
v20 <- fit2$t
p1 <- unlist(lapply(fit1$p, FUN = function(L) L[2]))
p2 <- unlist(lapply(fit2$p, FUN = function(L) L[2]))

par(mfrow = c(1, 2))
par(mar = c(2.5, 2.5, 1.5, 1.5))

plot(v10, v11, type = "l", lty = 2, xlab = "", ylab = "", main = "Hazard", col = "red")
lines(v10, 2/x1*log(1+x1*v10), col = "red")
lines(v20, v21, lty = 2, col = "blue")
lines(v20, 2/x2*log(1+x2*v20), col = "blue")

plot(v10, p1, type = "l", lty = 2, xlab = "", ylab = "", main = "Probability", col = "red")
lines(seq(0, 5, 0.01), P1, col = "red")
lines(v20, p2, lty = 2, col = "blue")
lines(seq(0, 5, 0.01), P2, col = "blue")

## -----------------------------------------------------------------------------
jump_rate <- function(i, t, u){
  if(i == 1){
    0.1 + 0.002*t
  } else if(i == 2){
    ifelse(u < 4, 0.29, 0.09) + 0.001*t
  } else{
    0
  }
}

mark_dist <- function(i, s, v){
  if(i == 1){
    c(0, 0.9, 0.1)
  } else if(i == 2){
    c(0, 0, 1)
  } else{
    0
  }
}

## -----------------------------------------------------------------------------
n <- 5000

c <- runif(n, 10, 40)

sim <- list()
for(i in 1:n){
  sim[[i]] <- sim_path(1, rates = jump_rate, dists = mark_dist, tn = c[i],
                       bs = c(0.1+0.002*c[i], 0.29+0.001*c[i], 0))
}

## -----------------------------------------------------------------------------
sum(c == unlist(lapply(sim, FUN = function(z){tail(z$times, 1)}))) / n

## -----------------------------------------------------------------------------
fit <- aalen_johansen(sim)

## ---- fig.align = 'center', fig.height = 4.5, fig.width = 4.5-----------------
v0 <- fit$t
p <- unlist(lapply(fit$p, FUN = function(L) L[2]))
integrand <- function(t, s){
  exp(-0.1*s-0.001*s^2)*(0.09 + 0.0018*s)*exp(-0.20*pmin(t-s, 4)-0.09*(t-s)-0.0005*(t^2-s^2))
}
P <- Vectorize(function(t){
  integrate(f = integrand, lower = 0, upper = t, t = t)$value
}, vectorize.args = "t")
plot(v0, p, type = "l", lty = 2, xlab = "", ylab = "", main = "Probability")
lines(seq(0, 40, 0.1), P(seq(0, 40, 0.1)))

## -----------------------------------------------------------------------------
landmark <- sim[unlist(lapply(sim, FUN = function(z){any(z$times <= 10
                                                         & c(z$times[-1], Inf) > 10
                                                         & z$states == 2)}))]
landmark <- lapply(landmark, FUN = function(z){list(times = z$times, states = z$states,
                                                    X = 10 - z$times[z$times <= 10
                                                                     & c(z$times[-1], Inf) > 10
                                                                     & z$states == 2])})

## -----------------------------------------------------------------------------
length(landmark) / n

## -----------------------------------------------------------------------------
u1 <- 1
u2 <- 5

fit1 <- aalen_johansen(landmark, x = u1)
fit2 <- aalen_johansen(landmark, x = u2)
fit3 <- aalen_johansen(landmark)

## ---- fig.align = 'center', fig.height = 4.5, fig.width = 4.5-----------------
v10 <- fit1$t
v20 <- fit2$t
v30 <- fit3$t
p1 <- unlist(lapply(fit1$p, FUN = function(L) L[2]))
p2 <- unlist(lapply(fit2$p, FUN = function(L) L[2]))
p3 <- unlist(lapply(fit3$p, FUN = function(L) L[2]))
P <- function(t, u){
  exp(-(t-10)*0.09-(t^2-100)*0.0005-pmax(0, pmin(t, 4-(u-10))-10)*0.20)  
} 

plot(v10, p1, type = "l", lty = 2, xlab = "", ylab = "", main = "Probability",
     col = "red", xlim = c(10, 40))
lines(seq(10, 40, 0.1), P(seq(10, 40, 0.1), u1), col = "red")
lines(v20, p2, lty = 2, col = "blue")
lines(seq(10, 40, 0.1), P(seq(10, 40, 0.1), u2), col = "blue")
lines(v30, p3, lty = 3)

