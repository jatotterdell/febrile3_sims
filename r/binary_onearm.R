stages <- 1:520
n_stages <- 1:520
n_new <- diff(c(0, n_stages))
x <- rbinom(n_new, n_new, 0.04)


theta0 <- 0.04     # Our estimate of the baseline rate
k <- 0.03          # Tolerance
theta_true <- 0.04 # The true rate under new policy
week <- 1:104
tt <- rep(1:104, each = 5)
nmax <- length(tt)
x <- rbinom(nmax, 1, theta_true)
threshold <- 0.975

in_analysis <- outer(tt, week, `<=`)
n <- apply(in_analysis, 2, function(a) sum(a))

a <- 1
b <- 1

# Earliest we could possibly stop
n_min_stop <- which.max(pbeta(theta0 + k, a, a + 1:nmax) > 0.975)
w_min_stop <- min(which(n >= n_min_stop))

# ANALYTICAL SOLUTION
prob_success <- rep(0, length(week))
for(w in week) {
  yprob <- dbinom(0:n[w], n[w], theta_true)
  p <- pbeta(theta0 + k, a + 0:n[w], b + n[w] - 0:n[w])
  success <- p > threshold
  prob_success[w] <- sum(success * yprob)
}


sims <- 100000
x <- matrix(0, sims, nmax)
y <- matrix(0, sims, tail(week, 1))
p <- matrix(0, sims, tail(week, 1))
d <- matrix(0, sims, tail(week, 1))
decision_week <- rep(0, sims)
decision_final <- rep(0, sims)

for(s in 1:sims) {
  x[s, ] <- rbinom(nmax, 1, theta_true)
  y[s, ] <- apply(in_analysis, 2, function(z) sum(x[s, z]))
  p[s, ] <- pbeta(theta0 + k, a + y[s, ], b + n - y[s, ])
  d[s, ] <- p[s, ] > threshold
  decision_week[s] <- ifelse(any(d[s, ] == 1), min(which(d[s, ] == 1)), NA)
  decision_final[s] <- tail(d[s, ], 1) 
}

mean(!is.na(decision_week))
table(decision_week)

decision_week_fill <- decision_week
decision_week_fill[is.na(decision_week)] <- tail(week, 1) + 1

prop.table(table(decision_week_fill))
cumsum(prop.table(table(decision_week_fill)))

mean(decision_week, na.rm = TRUE); median(decision_week, na.rm = TRUE)
mean(decision_week_fill)

par(mar = c(4, 6, 2, 2))
matplot(week, t(p), type = 'l', lty = 2, col = 1, ylab = expression(integral(pi^m*(theta), 0, theta[0]+k)), xlab = "Week")
lines(week, apply(t(p), 1, mean), lw = 2)
abline(h = 0.975, lty = 3)

par(mar = c(4, 6, 2, 2))
matplot(week, t(p[is.na(decision_week), ]), type = 'l', lty = 2, col = 1, ylab = expression(integral(pi^m*(theta), 0, theta[0]+k)), xlab = "Week")
lines(week, apply(t(p[is.na(decision_week), ]), 1, mean), lw = 2)
abline(h = 0.975)

plot(week, apply(d, 2, mean), type = "l")
lines(prob_success, col = "red")
plot(week, apply(apply(d, 1, cummax), 1, mean), type = 'l')



# Analytic
nmax <- 520
theta0 <- 0.07
theta_true <- 0.04
cond_p <- rep(0, nmax)
ystar <- rep(0, nmax)
for(i in 1:nmax) {
  p <- pbeta(theta0, a + 0:i, b + i - 0:i)
  ystar[i] <- min(which(p < threshold) - 1)
  if (i == 1) {
    cond_p[i] <- sum((pbeta(theta0, a + 0:1, b + 1 - 0:1) > threshold)*dbinom(0:1, 1, theta_true))
  }
  else {
    cond_p[i] <- sum((pbeta(theta0, a + ystar[i-1]:i, b + i - ystar[i-1]:i) > threshold)*dbinom(ystar[i-1]:i, i, theta_true))
  }
}

plot(cond_p, type = "l")
plot(1 - cond_p, type = 'l')
plot(cumprod(1 - cond_p), type = 'l')
