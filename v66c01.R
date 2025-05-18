## Code Figure 1
library("GoFKernel")
dev.new(width = 25, height = 9.5)
par(mfcol = c(1, 2), mai = c(0.5, 0.9, 0.3, 0.3))
set.seed(986)
x <- runif(500)
dr <- density.reflected(x, lower = 0, upper = 1)
y2 <- rep(1, length(dr$x))
plot(dr$x, dr$y, type = "l", lwd = 2.7, col = "blue", bty = "L", xlab = "", ylab = "density", 
  xlim = c(0, 1), ylim = c(0.9, 1.1), cex.axis = 0.85)
polygon(c(dr$x, rev(dr$x)), c(y2, rev(dr$y)), col = "lightgrey")
points(dr$x, y2, type = "l", col = "red", lwd = 3, lty = 3)

set.seed(4558)
x <- rnorm(50)
dr <- density.reflected(x, lower = -Inf, upper = Inf)
y2 <- dnorm(dr$x)
plot(dr$x, dr$y, type = "l", lwd = 2.7, col = "blue", bty = "L", xlab = "", ylab = "density", 
  ylim = c(0, 0.5), xlim = c(-4, 4), cex.axis = 0.85)
polygon(c(dr$x, rev(dr$x)), c(y2, rev(dr$y)), col = "lightgrey")
points(dr$x, y2, type = "l", col = "red", lwd = 3, lty = 3)


## Example 1 Fan's test
library("GoFKernel")
set.seed(125)
fan.test(runif(100), dunif, lower = 0, upper = 1)


## Example 2 Fan's test
library("GoFKernel")
f0 <- function(x) ifelse(x >= 0 & x <= 1, 2 - 2 * x, 0)
fan.test(risk76.1929, f0, lower = 0, upper = 1, kernel = "epanech")


## Example 1 dgeometric test
library("GoFKernel")
set.seed(158)
f0 <- function(x) ifelse(x >= 0 & x <= 1, 2 - 2 * x, 0)
dgeometric.test(risk76.1929, f0, lower = 0, upper = 1, n.sim = 51)


## Example 2 dgeometric test
library("GoFKernel")
library("MASS")
set.seed(1)
x <- rlnorm(200, meanlog = 1, sdlog = 1)
dgeometric.test(x, dgamma, par = lapply(fitdistr(x, "gamma")$estimate, function(i) i), 
  lower = 0, upper = Inf, n.sim = 125)


## Code Figure 2
dev.new(width = 40, height = 15)
par(mfcol = c(1, 2), mai = c(0.5, 0.9, 0.3, 0.3))
set.seed(789)
x <- runif(1000)
curve(dunif, from = 0, to = 1, ylab = "density", ylim = c(0.6, 1.2), xlab = "")
points(density(x, from = 0, to = 1), type = "l", col = "blue")

library("GoFKernel")
points(density.reflected(x, lower = 0, upper = 1), type = "l", col = "red")
segments(0, 1.2, 0.05, 1.2, col = "black")
text(0.28, 1.2, "Theoretical density: U(0,1)", cex = 0.75)
segments(0, 1.17, 0.05, 1.17, col = "red")
text(0.374, 1.17, "Kernel estimate using density.reflected", cex = 0.75)
segments(0, 1.14, 0.05, 1.14, col = "blue")
text(0.309, 1.14, "Kernel estimate using density", cex = 0.75)

set.seed(942)
f0 <- function(x) ifelse(x >= 0 & x <= 1, 2 - 2 * x, 0)
curve(f0, from = 0, to = 1, ylab = "density", ylim = c(0, 2.25), xlab = "")
x <- random.function(300, f0, 0, 1)
points(density.reflected(x, lower = 0, upper = 1), type = "l", col = "red")
segments(0, 2.25, 0.05, 2.25, col = "black")
text(0.3, 2.25, "Theoretical density: f(x)=2-2x", cex = 0.75)
segments(0, 2.13, 0.05, 2.13, col = "red")
text(0.379, 2.13, "Kernel estimate using density.reflected", cex = 0.75)


## Outputs of Section 5: A comparison of the ECDF and EKF approaches
## Tables 1 and 2
# The numbers in the tables are saved in the objects:
# SUMMARY.U.U: Block 1, Table 1. Uniform-Uniform.
# SUMMARY.N.N: Block 2, Table 1. Normal-Normal.
# SUMMARY.E.E: Block 3, Table 1. Exponential-Exponential.
# SUMMARY.L.L: Block 4, Table 1. LogNormal-LogNormal.
# SUMMARY.B.N: Block 1, Table 2. Beta-Uniform.
# SUMMARY.C.N: Block 2, Table 2. Cauchy-Normal.
# SUMMARY.G.E: Block 3, Table 2. Gamma-Exponential.
# SUMMARY.L.G: Block 4, Table 2. LogNormal-Gamma.

## FIGURES 3 to 10 The numbers to plot the figures are saved within the objects:
# unif.unif, normal.normal, exp.exp, lnorm.lnorm, beta.unif, cauchy.normal,
# gamma.exp, lnorm.gamma

## FIGURE 11 The object PTIME contains the numbers fro which to plot figure 11
# Note: The same samples analysed in section 5 are also analysed in section 6

# Packages
library("GoFKernel") # 2.0.3
library("ADGofTest") # 0.3
library("truncgof")  # 0.6-0

# Seed fixed
set.seed(1234)

# Processing Times
PTIME <- proc.time()


# -------------------------------------------------------
# Table 1: Block 2
# Normal. Actual distribution N(0,1). Null hypothesis N(0,1) #

# Data to be analysed sample size =10
normal_10 <- matrix(nrow = 1000, ncol = 10)
for (i in 1:1000) normal_10[i, ] <- rnorm(10)

# sample size =20
normal_20 <- matrix(nrow = 1000, ncol = 20)
for (i in 1:1000) normal_20[i, ] <- rnorm(20)

# sample size =50
normal_50 <- matrix(nrow = 1000, ncol = 50)
for (i in 1:1000) normal_50[i, ] <- rnorm(50)

# sample size =100
normal_100 <- matrix(nrow = 1000, ncol = 100)
for (i in 1:1000) normal_100[i, ] <- rnorm(100)

# sample size =200
normal_200 <- matrix(nrow = 1000, ncol = 200)
for (i in 1:1000) normal_200[i, ] <- rnorm(200)

# sample size =500
normal_500 <- matrix(nrow = 1000, ncol = 500)
for (i in 1:1000) normal_500[i, ] <- rnorm(500)

# Functions for hypothesis testing (we retain p-values)
ks.test.n <- function(x) return(tryCatch(stats::ks.test(x, "pnorm")$p.value, error = function(e) NA))
ad.test1 <- function(x) return(tryCatch(ADGofTest::ad.test(x, pnorm)$p.value, error = function(e) NA))
ad.test0 <- function(x) return(tryCatch(truncgof::ad.test(x, "pnorm", list(mean = 0, 
  sd = 1))$p.value, error = function(e) NA))
ad2.test.n <- function(x) return(tryCatch(ad2.test(x, "pnorm", list(mean = 0, sd = 1))$p.value, 
  error = function(e) NA))
ks.test.1 <- function(x) return(tryCatch(truncgof::ks.test(x, "pnorm", list(mean = 0, 
  sd = 1))$p.value, error = function(e) NA))
v.test.n <- function(x) return(tryCatch(v.test(x, "pnorm", list(mean = 0, sd = 1))$p.value, 
  error = function(e) NA))
w2.test.n <- function(x) return(tryCatch(w2.test(x, "pnorm", list(mean = 0, sd = 1))$p.value, 
  error = function(e) NA))
fan.n <- function(x) return(tryCatch(fan.test(x, dnorm)$p.value, error = function(e) NA))
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dnorm, n.sim = 101)$p.value, 
  error = function(e) NA))

# Object to save the p-values
tests <- c("ks.test.n", "ad.test1", "ad.test0", "ad2.test.n", "ks.test.1", "v.test.n", 
  "w2.test.n", "fan.n", "dgeometric.n")
normal.normal <- matrix(nrow = 1000, ncol = 54)
colnames(normal.normal) <- c(paste(tests, "_10", sep = ""), paste(tests, "_20", sep = ""), 
  paste(tests, "_50", sep = ""), paste(tests, "_100", sep = ""), paste(tests, "_200", 
    sep = ""), paste(tests, "_500", sep = ""))
PTIME <- c(PTIME, proc.time())

# Computing p-values
for (i in 1:9) {
  normal.normal[, i] <- apply(normal_10, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  normal.normal[, i + 9] <- apply(normal_20, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  normal.normal[, i + 18] <- apply(normal_50, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  normal.normal[, i + 27] <- apply(normal_100, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  normal.normal[, i + 36] <- apply(normal_200, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  normal.normal[, i + 45] <- apply(normal_500, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}

# Summary of scenario
output <- matrix(apply(normal.normal < 0.05, 2, sum, na.rm = T)/1000, ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
SUMMARY.N.N <- t(output)
output <- matrix(apply(is.na(normal.normal), 2, sum), ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
no.convergence.N.N <- t(output)


# -------------------------------------------------------
# Table 2: Block 2
# Cauchy-Normal. Actual distribution Ca(0,1). Null hypothesis N(0,sd(x)) #

# Data to be analysed sample size =10
cauchy_10 <- matrix(nrow = 1000, ncol = 10)
for (i in 1:1000) cauchy_10[i, ] <- rcauchy(10)

# sample size =20
cauchy_20 <- matrix(nrow = 1000, ncol = 20)
for (i in 1:1000) cauchy_20[i, ] <- rcauchy(20)

# sample size =50
cauchy_50 <- matrix(nrow = 1000, ncol = 50)
for (i in 1:1000) cauchy_50[i, ] <- rcauchy(50)

# sample size =100
cauchy_100 <- matrix(nrow = 1000, ncol = 100)
for (i in 1:1000) cauchy_100[i, ] <- rcauchy(100)

# sample size =200
cauchy_200 <- matrix(nrow = 1000, ncol = 200)
for (i in 1:1000) cauchy_200[i, ] <- rcauchy(200)

# sample size =500
cauchy_500 <- matrix(nrow = 1000, ncol = 500)
for (i in 1:1000) cauchy_500[i, ] <- rcauchy(500)

# Functions for hypothesis testing (we retain p-values)
ks.test.n <- function(x) return(tryCatch(stats::ks.test(x, "pnorm", sd = sd(x))$p.value, 
  error = function(e) NA))
ad.test1 <- function(x) return(tryCatch(ADGofTest::ad.test(x, pnorm, sd = sd(x))$p.value, 
  error = function(e) NA))
ad.test0 <- function(x) return(tryCatch(truncgof::ad.test(x, "pnorm", list(mean = 0, 
  sd = sd(x)))$p.value, error = function(e) NA))
ad2.test.n <- function(x) return(tryCatch(ad2.test(x, "pnorm", list(mean = 0, sd = sd(x)))$p.value, 
  error = function(e) NA))
ks.test.1 <- function(x) return(tryCatch(truncgof::ks.test(x, "pnorm", list(mean = 0, 
  sd = sd(x)))$p.value, error = function(e) NA))
v.test.n <- function(x) return(tryCatch(v.test(x, "pnorm", list(mean = 0, sd = sd(x)))$p.value, 
  error = function(e) NA))
w2.test.n <- function(x) return(tryCatch(w2.test(x, "pnorm", list(mean = 0, sd = sd(x)))$p.value, 
  error = function(e) NA))
fan.n <- function(x) return(tryCatch(fan.test(x, dnorm, list(mean = 0, sd = sd(x)))$p.value, 
  error = function(e) NA))
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dnorm, list(mean = 0, 
  sd = sd(x)), n.sim = 101)$p.value, error = function(e) NA))

# Object to save the p-values
tests <- c("ks.test.n", "ad.test1", "ad.test0", "ad2.test.n", "ks.test.1", "v.test.n", 
  "w2.test.n", "fan.n", "dgeometric.n")
cauchy.normal <- matrix(nrow = 1000, ncol = 54)
colnames(cauchy.normal) <- c(paste(tests, "_10", sep = ""), paste(tests, "_20", sep = ""), 
  paste(tests, "_50", sep = ""), paste(tests, "_100", sep = ""), paste(tests, "_200", 
    sep = ""), paste(tests, "_500", sep = ""))
PTIME <- c(PTIME, proc.time())

# Computing p-values
for (i in 1:9) {
  cauchy.normal[, i] <- apply(cauchy_10, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  cauchy.normal[, i + 9] <- apply(cauchy_20, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  cauchy.normal[, i + 18] <- apply(cauchy_50, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  cauchy.normal[, i + 27] <- apply(cauchy_100, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  cauchy.normal[, i + 36] <- apply(cauchy_200, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  cauchy.normal[, i + 45] <- apply(cauchy_500, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}

# Summary of scenario
output <- matrix(apply(cauchy.normal < 0.05, 2, sum, na.rm = T)/apply(!is.na(cauchy.normal), 
  2, sum), ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
SUMMARY.C.N <- t(output)
output <- matrix(apply(is.na(cauchy.normal), 2, sum), ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
no.convergence.C.N <- t(output)


# -------------------------------------------------------
# Table 1: Block 1
# Uniform. Actual distribution U(0,1). Null hypothesis U(0,1) #

# Data to be analysed sample size =10
unif_10 <- matrix(nrow = 1000, ncol = 10)
for (i in 1:1000) unif_10[i, ] <- runif(10)

# sample size =20
unif_20 <- matrix(nrow = 1000, ncol = 20)
for (i in 1:1000) unif_20[i, ] <- runif(20)

# sample size =50
unif_50 <- matrix(nrow = 1000, ncol = 50)
for (i in 1:1000) unif_50[i, ] <- runif(50)

# sample size =100
unif_100 <- matrix(nrow = 1000, ncol = 100)
for (i in 1:1000) unif_100[i, ] <- runif(100)

# sample size =200
unif_200 <- matrix(nrow = 1000, ncol = 200)
for (i in 1:1000) unif_200[i, ] <- runif(200)

# sample size =500
unif_500 <- matrix(nrow = 1000, ncol = 500)
for (i in 1:1000) unif_500[i, ] <- runif(500)

# Functions for hypothesis testing (we retain p-values)
ks.test.n <- function(x) return(tryCatch(stats::ks.test(x, "punif")$p.value, error = function(e) NA))
ad.test1 <- function(x) return(tryCatch(ADGofTest::ad.test(x, punif)$p.value, error = function(e) NA))
ad.test0 <- function(x) return(tryCatch(truncgof::ad.test(x, "punif", list(min = 0, 
  max = 1))$p.value, error = function(e) NA))
ad2.test.n <- function(x) return(tryCatch(ad2.test(x, "punif", list(min = 0, max = 1))$p.value, 
  error = function(e) NA))
ks.test.1 <- function(x) return(tryCatch(truncgof::ks.test(x, "punif", list(min = 0, 
  max = 1))$p.value, error = function(e) NA))
v.test.n <- function(x) return(tryCatch(v.test(x, "punif", list(min = 0, max = 1))$p.value, 
  error = function(e) NA))
w2.test.n <- function(x) return(tryCatch(w2.test(x, "punif", list(min = 0, max = 1))$p.value, 
  error = function(e) NA))
fan.n <- function(x) return(tryCatch(fan.test(x, dunif, lower = 0, upper = 1)$p.value, 
  error = function(e) NA))
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dunif, lower = 0, 
  upper = 1, n.sim = 101)$p.value, error = function(e) NA))

# Object to save the p-values
tests <- c("ks.test.n", "ad.test1", "ad.test0", "ad2.test.n", "ks.test.1", "v.test.n", 
  "w2.test.n", "fan.n", "dgeometric.n")
unif.unif <- matrix(nrow = 1000, ncol = 54)
colnames(unif.unif) <- c(paste(tests, "_10", sep = ""), paste(tests, "_20", sep = ""), 
  paste(tests, "_50", sep = ""), paste(tests, "_100", sep = ""), paste(tests, "_200", 
    sep = ""), paste(tests, "_500", sep = ""))
PTIME <- c(PTIME, proc.time())

# Computing p-values
for (i in 1:9) {
  unif.unif[, i] <- apply(unif_10, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  unif.unif[, i + 9] <- apply(unif_20, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  unif.unif[, i + 18] <- apply(unif_50, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  unif.unif[, i + 27] <- apply(unif_100, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  unif.unif[, i + 36] <- apply(unif_200, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  unif.unif[, i + 45] <- apply(unif_500, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}

# Summary of scenario
output <- matrix(apply(unif.unif < 0.05, 2, sum, na.rm = T)/1000, ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
SUMMARY.U.U <- t(output)
output <- matrix(apply(is.na(unif.unif), 2, sum), ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
no.convergence.U.U <- t(output)


# -------------------------------------------------------
# Table 2: Block 1 #
# Beta-Uniform. Actual distribution Be(1.3,1.3). Null hypothesis U(0,1) #

# Data to be analysed sample size =10
beta_10 <- matrix(nrow = 1000, ncol = 10)
for (i in 1:1000) beta_10[i, ] <- rbeta(10, 1.3, 1.3)

# sample size =20
beta_20 <- matrix(nrow = 1000, ncol = 20)
for (i in 1:1000) beta_20[i, ] <- rbeta(20, 1.3, 1.3)

# sample size =50
beta_50 <- matrix(nrow = 1000, ncol = 50)
for (i in 1:1000) beta_50[i, ] <- rbeta(50, 1.3, 1.3)

# sample size =100
beta_100 <- matrix(nrow = 1000, ncol = 100)
for (i in 1:1000) beta_100[i, ] <- rbeta(100, 1.3, 1.3)

# sample size =200
beta_200 <- matrix(nrow = 1000, ncol = 200)
for (i in 1:1000) beta_200[i, ] <- rbeta(200, 1.3, 1.3)

# sample size =500
beta_500 <- matrix(nrow = 1000, ncol = 500)
for (i in 1:1000) beta_500[i, ] <- rbeta(500, 1.3, 1.3)

# Functions for hypothesis testing (we retain p-values)
ks.test.n <- function(x) return(tryCatch(stats::ks.test(x, "punif")$p.value, error = function(e) NA))
ad.test1 <- function(x) return(tryCatch(ADGofTest::ad.test(x, punif)$p.value, error = function(e) NA))
ad.test0 <- function(x) return(tryCatch(truncgof::ad.test(x, "punif", list(min = 0, 
  max = 1))$p.value, error = function(e) NA))
ad2.test.n <- function(x) return(tryCatch(ad2.test(x, "punif", list(min = 0, max = 1))$p.value, 
  error = function(e) NA))
ks.test.1 <- function(x) return(tryCatch(truncgof::ks.test(x, "punif", list(min = 0, 
  max = 1))$p.value, error = function(e) NA))
v.test.n <- function(x) return(tryCatch(v.test(x, "punif", list(min = 0, max = 1))$p.value, 
  error = function(e) NA))
w2.test.n <- function(x) return(tryCatch(w2.test(x, "punif", list(min = 0, max = 1))$p.value, 
  error = function(e) NA))
fan.n <- function(x) return(tryCatch(fan.test(x, dunif, lower = 0, upper = 1)$p.value, 
  error = function(e) NA))
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dunif, lower = 0, 
  upper = 1, n.sim = 101)$p.value, error = function(e) NA))

# Object to save the p-values
tests <- c("ks.test.n", "ad.test1", "ad.test0", "ad2.test.n", "ks.test.1", "v.test.n", 
  "w2.test.n", "fan.n", "dgeometric.n")
beta.unif <- matrix(nrow = 1000, ncol = 54)
colnames(beta.unif) <- c(paste(tests, "_10", sep = ""), paste(tests, "_20", sep = ""), 
  paste(tests, "_50", sep = ""), paste(tests, "_100", sep = ""), paste(tests, "_200", 
    sep = ""), paste(tests, "_500", sep = ""))
PTIME <- c(PTIME, proc.time())

# Computing p-values
for (i in 1:9) {
  beta.unif[, i] <- apply(beta_10, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  beta.unif[, i + 9] <- apply(beta_20, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  beta.unif[, i + 18] <- apply(beta_50, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  beta.unif[, i + 27] <- apply(beta_100, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  beta.unif[, i + 36] <- apply(beta_200, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  beta.unif[, i + 45] <- apply(beta_500, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}

# Summary of scenario
output <- matrix(apply(beta.unif < 0.05, 2, sum, na.rm = T)/1000, ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
SUMMARY.B.U <- t(output)
output <- matrix(apply(is.na(beta.unif), 2, sum), ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
no.convergence.B.U <- t(output)


# -------------------------------------------------------
# Table 1: Block 3
# Exponential. Actual distribution Exp(1). Null hypothesis Exp(1) #

# Data to be analysed sample size =10
exp_10 <- matrix(nrow = 1000, ncol = 10)
for (i in 1:1000) exp_10[i, ] <- rexp(10)

# sample size =20
exp_20 <- matrix(nrow = 1000, ncol = 20)
for (i in 1:1000) exp_20[i, ] <- rexp(20)

# sample size =50
exp_50 <- matrix(nrow = 1000, ncol = 50)
for (i in 1:1000) exp_50[i, ] <- rexp(50)

# sample size =100
exp_100 <- matrix(nrow = 1000, ncol = 100)
for (i in 1:1000) exp_100[i, ] <- rexp(100)

# sample size =200
exp_200 <- matrix(nrow = 1000, ncol = 200)
for (i in 1:1000) exp_200[i, ] <- rexp(200)

# sample size =500
exp_500 <- matrix(nrow = 1000, ncol = 500)
for (i in 1:1000) exp_500[i, ] <- rexp(500)

# Functions for hypothesis testing (we retain p-values)
ks.test.n <- function(x) return(tryCatch(stats::ks.test(x, "pexp")$p.value, error = function(e) NA))
ad.test1 <- function(x) return(tryCatch(ADGofTest::ad.test(x, pexp)$p.value, error = function(e) NA))
ad.test0 <- function(x) return(tryCatch(truncgof::ad.test(x, "pexp", list(rate = 1))$p.value, 
  error = function(e) NA))
ad2.test.n <- function(x) return(tryCatch(ad2.test(x, "pexp", list(rate = 1))$p.value, 
  error = function(e) NA))
ks.test.1 <- function(x) return(tryCatch(truncgof::ks.test(x, "pexp", list(rate = 1))$p.value, 
  error = function(e) NA))
v.test.n <- function(x) return(tryCatch(v.test(x, "pexp", list(rate = 1))$p.value, 
  error = function(e) NA))
w2.test.n <- function(x) return(tryCatch(w2.test(x, "pexp", list(rate = 1))$p.value, 
  error = function(e) NA))
fan.n <- function(x) return(tryCatch(fan.test(x, dexp, lower = 0)$p.value, error = function(e) NA))
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dexp, lower = 0, n.sim = 101)$p.value, 
  error = function(e) NA))

# Object to save the p-values
tests <- c("ks.test.n", "ad.test1", "ad.test0", "ad2.test.n", "ks.test.1", "v.test.n", 
  "w2.test.n", "fan.n", "dgeometric.n")
exp.exp <- matrix(nrow = 1000, ncol = 54)
colnames(exp.exp) <- c(paste(tests, "_10", sep = ""), paste(tests, "_20", sep = ""), 
  paste(tests, "_50", sep = ""), paste(tests, "_100", sep = ""), paste(tests, "_200", 
    sep = ""), paste(tests, "_500", sep = ""))
PTIME <- c(PTIME, proc.time())

# Computing p-values
for (i in 1:9) {
  exp.exp[, i] <- apply(exp_10, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  exp.exp[, i + 9] <- apply(exp_20, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  exp.exp[, i + 18] <- apply(exp_50, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  exp.exp[, i + 27] <- apply(exp_100, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  exp.exp[, i + 36] <- apply(exp_200, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  exp.exp[, i + 45] <- apply(exp_500, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}

# Summary of scenario
output <- matrix(apply(exp.exp < 0.05, 2, sum, na.rm = T)/1000, ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
SUMMARY.E.E <- t(output)
output <- matrix(apply(is.na(exp.exp), 2, sum), ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
no.convergence.E.E <- t(output)


# -------------------------------------------------------
# Table 2: Block 3
# Gamma-Exponential. Actual distribution Ga(0.9,0.9). Null hypothesis Exp(1) #

# Data to be analysed sample size =10
gamma_10 <- matrix(nrow = 1000, ncol = 10)
for (i in 1:1000) gamma_10[i, ] <- rgamma(10, 0.9, 0.9)

# sample size =20
gamma_20 <- matrix(nrow = 1000, ncol = 20)
for (i in 1:1000) gamma_20[i, ] <- rgamma(20, 0.9, 0.9)

# sample size =50
gamma_50 <- matrix(nrow = 1000, ncol = 50)
for (i in 1:1000) gamma_50[i, ] <- rgamma(50, 0.9, 0.9)

# sample size =100
gamma_100 <- matrix(nrow = 1000, ncol = 100)
for (i in 1:1000) gamma_100[i, ] <- rgamma(100, 0.9, 0.9)

# sample size =200
gamma_200 <- matrix(nrow = 1000, ncol = 200)
for (i in 1:1000) gamma_200[i, ] <- rgamma(200, 0.9, 0.9)

# sample size =500
gamma_500 <- matrix(nrow = 1000, ncol = 500)
for (i in 1:1000) gamma_500[i, ] <- rgamma(500, 0.9, 0.9)

# Functions for hypothesis testing (we retain p-values)
ks.test.n <- function(x) return(tryCatch(stats::ks.test(x, "pexp")$p.value, error = function(e) NA))
ad.test1 <- function(x) return(tryCatch(ADGofTest::ad.test(x, pexp)$p.value, error = function(e) NA))
ad.test0 <- function(x) return(tryCatch(truncgof::ad.test(x, "pexp", list(rate = 1))$p.value, 
  error = function(e) NA))
ad2.test.n <- function(x) return(tryCatch(ad2.test(x, "pexp", list(rate = 1))$p.value, 
  error = function(e) NA))
ks.test.1 <- function(x) return(tryCatch(truncgof::ks.test(x, "pexp", list(rate = 1))$p.value, 
  error = function(e) NA))
v.test.n <- function(x) return(tryCatch(v.test(x, "pexp", list(rate = 1))$p.value, 
  error = function(e) NA))
w2.test.n <- function(x) return(tryCatch(w2.test(x, "pexp", list(rate = 1))$p.value, 
  error = function(e) NA))
fan.n <- function(x) return(tryCatch(fan.test(x, dexp, lower = 0)$p.value, error = function(e) NA))
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dexp, lower = 0, n.sim = 101)$p.value, 
  error = function(e) NA))

# Object to save the p-values
tests <- c("ks.test.n", "ad.test1", "ad.test0", "ad2.test.n", "ks.test.1", "v.test.n", 
  "w2.test.n", "fan.n", "dgeometric.n")
gamma.exp <- matrix(nrow = 1000, ncol = 54)
colnames(gamma.exp) <- c(paste(tests, "_10", sep = ""), paste(tests, "_20", sep = ""), 
  paste(tests, "_50", sep = ""), paste(tests, "_100", sep = ""), paste(tests, "_200", 
    sep = ""), paste(tests, "_500", sep = ""))
PTIME <- c(PTIME, proc.time())

# Computing p-values
for (i in 1:9) {
  gamma.exp[, i] <- apply(gamma_10, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  gamma.exp[, i + 9] <- apply(gamma_20, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  gamma.exp[, i + 18] <- apply(gamma_50, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  gamma.exp[, i + 27] <- apply(gamma_100, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  gamma.exp[, i + 36] <- apply(gamma_200, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  gamma.exp[, i + 45] <- apply(gamma_500, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}

# Summary of scenario
output <- matrix(apply(gamma.exp < 0.05, 2, sum, na.rm = T)/1000, ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
SUMMARY.G.E <- t(output)
output <- matrix(apply(is.na(gamma.exp), 2, sum), ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
no.convergence.G.E <- t(output)


# -------------------------------------------------------
# Table 1: Block 4 #
# LogNormal. Actual distribution logN(1,1). Null hypothesis LogNormal(1,1)#
library("MASS")

# Data to be analysed sample size =10
lnorm_10 <- matrix(nrow = 1000, ncol = 10)
for (i in 1:1000) lnorm_10[i, ] <- rlnorm(10, 1, 1)

# sample size =20
lnorm_20 <- matrix(nrow = 1000, ncol = 20)
for (i in 1:1000) lnorm_20[i, ] <- rlnorm(20, 1, 1)

# sample size =50
lnorm_50 <- matrix(nrow = 1000, ncol = 50)
for (i in 1:1000) lnorm_50[i, ] <- rlnorm(50, 1, 1)

# sample size =100
lnorm_100 <- matrix(nrow = 1000, ncol = 100)
for (i in 1:1000) lnorm_100[i, ] <- rlnorm(100, 1, 1)

# sample size =200
lnorm_200 <- matrix(nrow = 1000, ncol = 200)
for (i in 1:1000) lnorm_200[i, ] <- rlnorm(200, 1, 1)

# sample size =500
lnorm_500 <- matrix(nrow = 1000, ncol = 500)
for (i in 1:1000) lnorm_500[i, ] <- rlnorm(500, 1, 1)

# Functions for hypothesis testing (we retain p-values) functions
ks.test.n <- function(x) return(tryCatch(stats::ks.test(x, "plnorm", meanlog = 1, 
  sdlog = 1)$p.value, error = function(e) NA))
ad.test1 <- function(x) return(tryCatch(ADGofTest::ad.test(x, plnorm, meanlog = 1, 
  sdlog = 1)$p.value, error = function(e) NA))
ad.test0 <- function(x) return(tryCatch(truncgof::ad.test(x, "plnorm", list(meanlog = 1, 
  sdlog = 1))$p.value, error = function(e) NA))
ad2.test.n <- function(x) return(tryCatch(ad2.test(x, "plnorm", list(meanlog = 1, 
  sdlog = 1))$p.value, error = function(e) NA))
ks.test.1 <- function(x) return(tryCatch(truncgof::ks.test(x, "plnorm", list(meanlog = 1, 
  sdlog = 1))$p.value, error = function(e) NA))
v.test.n <- function(x) return(tryCatch(v.test(x, "plnorm", list(meanlog = 1, sdlog = 1))$p.value, 
  error = function(e) NA))
w2.test.n <- function(x) return(tryCatch(w2.test(x, "plnorm", list(meanlog = 1, sdlog = 1))$p.value, 
  error = function(e) NA))
fan.n <- function(x) return(tryCatch(fan.test(x, dlnorm, list(meanlog = 1, sdlog = 1), 
  lower = 0)$p.value, error = function(e) NA))
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dlnorm, list(meanlog = 1, 
  sdlog = 1), lower = 0, n.sim = 101)$p.value, error = function(e) NA))

# Object to save the p-values
tests <- c("ks.test.n", "ad.test1", "ad.test0", "ad2.test.n", "ks.test.1", "v.test.n", 
  "w2.test.n", "fan.n", "dgeometric.n")
lnorm.lnorm <- matrix(nrow = 1000, ncol = 54)
colnames(lnorm.lnorm) <- c(paste(tests, "_10", sep = ""), paste(tests, "_20", sep = ""), 
  paste(tests, "_50", sep = ""), paste(tests, "_100", sep = ""), paste(tests, "_200", 
    sep = ""), paste(tests, "_500", sep = ""))
PTIME <- c(PTIME, proc.time())

# Computing p-values
for (i in 1:9) {
  lnorm.lnorm[, i] <- apply(lnorm_10, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  lnorm.lnorm[, i + 9] <- apply(lnorm_20, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  lnorm.lnorm[, i + 18] <- apply(lnorm_50, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  lnorm.lnorm[, i + 27] <- apply(lnorm_100, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  lnorm.lnorm[, i + 36] <- apply(lnorm_200, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  lnorm.lnorm[, i + 45] <- apply(lnorm_500, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}

# Summary of scenario
output <- matrix(apply(lnorm.lnorm < 0.05, 2, sum, na.rm = T)/1000, ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
SUMMARY.L.L <- t(output)
output <- matrix(apply(is.na(lnorm.lnorm), 2, sum), ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
no.convergence.L.L <- t(output)


# -------------------------------------------------------
# Table 2: Block 4
# LogNormal-Gamma. Actual distribution logN(1,1). Null hypothesis
# Ga(est(alfa),est(beta))#

# Data to be analysed, the same generated in the previous scenario

# Functions for hypothesis testing (we retain p-values) functions
ks.test.n <- function(x) return(tryCatch(stats::ks.test(x, "pgamma", fitdistr(x, 
  "gamma")$estimate[1], fitdistr(x, "gamma")$estimate[2])$p.value, error = function(e) NA))
ad.test1 <- function(x) return(tryCatch(ADGofTest::ad.test(x, pgamma, fitdistr(x, 
  "gamma")$estimate[1], fitdistr(x, "gamma")$estimate[2])$p.value, error = function(e) NA))
ad.test0 <- function(x) return(tryCatch(truncgof::ad.test(x, "pgamma", lapply(fitdistr(x, 
  "gamma")$estimate, function(i) i))$p.value, error = function(e) NA))
ad2.test.n <- function(x) return(tryCatch(ad2.test(x, "pgamma", lapply(fitdistr(x, 
  "gamma")$estimate, function(i) i))$p.value, error = function(e) NA))
ks.test.1 <- function(x) return(tryCatch(truncgof::ks.test(x, "pgamma", lapply(fitdistr(x, 
  "gamma")$estimate, function(i) i))$p.value, error = function(e) NA))
v.test.n <- function(x) return(tryCatch(v.test(x, "pgamma", lapply(fitdistr(x, "gamma")$estimate, 
  function(i) i))$p.value, error = function(e) NA))
w2.test.n <- function(x) return(tryCatch(w2.test(x, "pgamma", lapply(fitdistr(x, 
  "gamma")$estimate, function(i) i))$p.value, error = function(e) NA))
fan.n <- function(x) return(tryCatch(fan.test(x, dgamma, lapply(fitdistr(x, "gamma")$estimate, 
  function(i) i), lower = 0)$p.value, error = function(e) NA))
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dgamma, lapply(fitdistr(x, 
  "gamma")$estimate, function(i) i), lower = 0, n.sim = 101)$p.value, error = function(e) NA))

# Object to save the p-values
tests <- c("ks.test.n", "ad.test1", "ad.test0", "ad2.test.n", "ks.test.1", "v.test.n", 
  "w2.test.n", "fan.n", "dgeometric.n")
lnorm.gamma <- matrix(nrow = 1000, ncol = 54)
colnames(lnorm.gamma) <- c(paste(tests, "_10", sep = ""), paste(tests, "_20", sep = ""), 
  paste(tests, "_50", sep = ""), paste(tests, "_100", sep = ""), paste(tests, "_200", 
    sep = ""), paste(tests, "_500", sep = ""))
PTIME <- c(PTIME, proc.time())

# Computing p-values
for (i in 1:9) {
  lnorm.gamma[, i] <- apply(lnorm_10, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  lnorm.gamma[, i + 9] <- apply(lnorm_20, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  lnorm.gamma[, i + 18] <- apply(lnorm_50, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  lnorm.gamma[, i + 27] <- apply(lnorm_100, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  lnorm.gamma[, i + 36] <- apply(lnorm_200, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}
for (i in 1:9) {
  lnorm.gamma[, i + 45] <- apply(lnorm_500, 1, tests[i])
  PTIME <- c(PTIME, proc.time())
}

# Summary of scenario
output <- matrix(apply(lnorm.gamma < 0.05, 2, sum, na.rm = T)/apply(!is.na(lnorm.gamma), 
  2, sum), ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
SUMMARY.L.G <- t(output)
output <- matrix(apply(is.na(lnorm.gamma), 2, sum), ncol = 9, byrow = T)
rownames(output) <- c(10, 20, 50, 100, 200, 500)
colnames(output) <- c("stats::ks.test", "ADGofTest::ad.test", "truncgof::ad.test", 
  "ad2.test", "truncgof::ks.test", "v.test", "w2.test", "Fan.test", "dgeometric.test")
no.convergence.L.G <- t(output)


## Code Figure 3
tests <- c("ks.test {stats}", "ad.test {ADGofTest}", "ad.test {truncgof}", "ad2.test {truncgof}", 
  "ks.test {truncgof}", "v.test {truncgof}", "w2.test {truncgof}", "fan.test {GoFKernel}", 
  "dgeometric.test {GoFKernel}")
dev.new(width = 18, height = 10.3)
par(mfrow = c(3, 3), mai = c(0.5, 0.5, 0.1, 0.1))
for (i in 1:9) {
  den.est <- density.reflected(unif.unif[, 9 + i], 0, 1)
  y.max <- max(den.est$y)
  plot(den.est, xlim = c(0, 1), ylim = c(0, y.max + 0.2), main = "", xlab = "")
  mtext(tests[i], 1, line = 2.5, cex = 0.75)
  abline(h = 1, col = "red", lty = 3)
}


## Code Figure 4
tests <- c("ks.test {stats}", "ad.test {ADGofTest}", "ad.test {truncgof}", "ad2.test {truncgof}", 
  "ks.test {truncgof}", "v.test {truncgof}", "w2.test {truncgof}", "fan.test {GoFKernel}", 
  "dgeometric.test {GoFKernel}")
dev.new(width = 18, height = 10.3)
par(mfrow = c(3, 3), mai = c(0.5, 0.5, 0.1, 0.1))
for (i in 1:9) {
  den.est <- density.reflected(normal.normal[, 18 + i], 0, 1)
  y.max <- max(den.est$y)
  plot(den.est, xlim = c(0, 1), ylim = c(0, y.max + 0.2), main = "", xlab = "")
  mtext(tests[i], 1, line = 2.5, cex = 0.75)
  abline(h = 1, col = "red", lty = 3)
}


## Code Figure 5
tests <- c("ks.test {stats}", "ad.test {ADGofTest}", "ad.test {truncgof}", "ad2.test {truncgof}", 
  "ks.test {truncgof}", "v.test {truncgof}", "w2.test {truncgof}", "fan.test {GoFKernel}", 
  "dgeometric.test {GoFKernel}")
dev.new(width = 18, height = 10.3)
par(mfrow = c(3, 3), mai = c(0.5, 0.5, 0.1, 0.1))
for (i in 1:9) {
  den.est <- density.reflected(exp.exp[, 27 + i], 0, 1)
  y.max <- max(den.est$y)
  plot(den.est, xlim = c(0, 1), ylim = c(0, y.max + 0.2), main = "", xlab = "")
  mtext(tests[i], 1, line = 2.5, cex = 0.75)
  abline(h = 1, col = "red", lty = 3)
}


## Code Figure 6
tests <- c("ks.test {stats}", "ad.test {ADGofTest}", "ad.test {truncgof}", "ad2.test {truncgof}", 
  "ks.test {truncgof}", "v.test {truncgof}", "w2.test {truncgof}", "fan.test {GoFKernel}", 
  "dgeometric.test {GoFKernel}")
dev.new(width = 18, height = 10.3)
par(mfrow = c(3, 3), mai = c(0.5, 0.5, 0.1, 0.1))
for (i in 1:9) {
  den.est <- density.reflected(lnorm.lnorm[, 36 + i], 0, 1)
  y.max <- max(den.est$y)
  plot(den.est, xlim = c(0, 1), ylim = c(0, y.max + 0.2), main = "", xlab = "")
  mtext(tests[i], 1, line = 2.5, cex = 0.75)
  abline(h = 1, col = "red", lty = 3)
}


## Code Figure 7
tests <- c("ks.test {stats}", "ad.test {ADGofTest}", "ad.test {truncgof}", "ad2.test {truncgof}", 
  "ks.test {truncgof}", "v.test {truncgof}", "w2.test {truncgof}", "fan.test {GoFKernel}", 
  "dgeometric.test {GoFKernel}")
dev.new(width = 18, height = 10.3)
par(mfrow = c(3, 3), mai = c(0.5, 0.5, 0.1, 0.1))
for (i in 1:9) {
  den.est <- density.reflected(beta.unif[, 27 + i], 0, 1)
  y.max <- max(den.est$y)
  plot(den.est, xlim = c(0, 1), ylim = c(0, y.max + 0.2), main = "", xlab = "")
  mtext(tests[i], 1, line = 2.5, cex = 0.75)
  f <- function(x) dbeta(x, 1, 100)
  curve(f, 0, 1, add = T, col = "red", lty = 3)
}


## Code Figure 8
tests <- c("ks.test {stats}", "ad.test {ADGofTest}", "ad.test {truncgof}", "ad2.test {truncgof}", 
  "ks.test {truncgof}", "v.test {truncgof}", "w2.test {truncgof}", "fan.test {GoFKernel}", 
  "dgeometric.test {GoFKernel}")
dev.new(width = 18, height = 10.3)
par(mfrow = c(3, 3), mai = c(0.5, 0.5, 0.1, 0.1))
for (i in 1:9) {
  x <- c(na.omit(cauchy.normal[, 0 + i]))
  den.est <- density.reflected(x, 0, 1)
  y.max <- max(den.est$y)
  plot(den.est, xlim = c(0, 1), ylim = c(0, y.max + 0.2), main = "", xlab = "")
  mtext(tests[i], 1, line = 2.5, cex = 0.75)
  f <- function(x) dbeta(x, 1, 100)
  curve(f, 0, 1, add = T, col = "red", lty = 3)
}


## Code Figure 9
tests <- c("ks.test {stats}", "ad.test {ADGofTest}", "ad.test {truncgof}", "ad2.test {truncgof}", 
  "ks.test {truncgof}", "v.test {truncgof}", "w2.test {truncgof}", "fan.test {GoFKernel}", 
  "dgeometric.test {GoFKernel}")
dev.new(width = 18, height = 10.3)
par(mfrow = c(3, 3), mai = c(0.5, 0.5, 0.1, 0.1))
for (i in 1:9) {
  x <- c(na.omit(gamma.exp[, 9 + i]))
  den.est <- density.reflected(x, 0, 1)
  y.max <- max(den.est$y)
  plot(den.est, xlim = c(0, 1), ylim = c(0, y.max + 0.2), main = "", xlab = "")
  mtext(tests[i], 1, line = 2.5, cex = 0.75)
  f <- function(x) dbeta(x, 1, 100)
  curve(f, 0, 1, add = T, col = "red", lty = 3)
}


## Code Figure 10
tests <- c("ks.test {stats}", "ad.test {ADGofTest}", "ad.test {truncgof}", "ad2.test {truncgof}", 
  "ks.test {truncgof}", "v.test {truncgof}", "w2.test {truncgof}", "fan.test {GoFKernel}", 
  "dgeometric.test {GoFKernel}")
dev.new(width = 18, height = 10.3)
par(mfrow = c(3, 3), mai = c(0.5, 0.5, 0.1, 0.1))
for (i in 1:9) {
  x <- c(na.omit(lnorm.gamma[, 45 + i]))
  den.est <- density.reflected(x, 0, 1)
  y.max <- max(den.est$y)
  plot(den.est, xlim = c(0, 1), ylim = c(0, y.max + 0.2), main = "", xlab = "")
  mtext(tests[i], 1, line = 2.5, cex = 0.75)
  f <- function(x) dbeta(x, 1, 100)
  curve(f, 0, 1, add = T, col = "red", lty = 3)
}


## Code Figure 11
tiempos <- matrix(PTIME[seq(8, 2203, 5)], ncol = 8)
tiempos <- tiempos[-1, ] - tiempos[-nrow(tiempos), ]
tiempos <- tiempos[, c(3, 4, 1, 2, 5:8)]/1000  # Average time
TP <- matrix(NA, nrow = 4 * 6, ncol = 8)
for (j in 1:6) {
  TP[1 + (j - 1) * 4, ] <- apply(tiempos[(1 + (j - 1) * 9):(2 + (j - 1) * 9), ], 
    2, mean)
  TP[2 + (j - 1) * 4, ] <- apply(tiempos[(3 + (j - 1) * 9):(7 + (j - 1) * 9), ], 
    2, mean)
  TP[3 + (j - 1) * 4, ] <- tiempos[8 + (j - 1) * 9, ]
  TP[4 + (j - 1) * 4, ] <- tiempos[9 + (j - 1) * 9, ]
}
Escenarios <- c("U(0,1) scenarios", "Be(1.3,1.3)-U(0,1) scenarios", "N(0,1) scenarios", 
  "", "Exp(1) scenarios", "Ga(0.9,0.9)-Exp(1) scenarios", "LogN(1,1) scenarios", 
  "")
texto <- c(1:3, 5:7)
tamanyos <- c(10, 20, 50, 100, 200, 500)
colores <- rep(c("yellow", "blue", "darkgreen", "red"), times = 2)
grupo <- c("sup-tests       ", "trunc-tests      ", "fan.test             ", "dgeometric-test", 
  "sup-tests       ", "trunc-tests      ", "fan.test             ", "dgeometric-test")
dev.new(width = 18, height = 16, units = "cm")
par(mfrow = c(4, 2), mai = c(0.55, 0.5, 0.1, 0.1))
for (i in 1:8) {
  datos <- TP[, i]
  plot(c(0, 500), c(0, 1.05 * max(datos)), type = "n", xlab = "", ylab = "")
  mtext("Sample size", side = 1, line = 1.8, cex = 0.6)
  mtext("Average Time per test (secs)", side = 2, line = 2.2, cex = 0.53)
  for (j in sample(1:4)) {
    y <- datos[c(j, j + 4, j + 8, j + 12, j + 16, j + 20)]
    points(tamanyos, y, col = colores[j], type = "l")
    points(tamanyos, y, col = colores[j], pch = 20)
  }
  points(0, max(datos), col = colores[i], cex = 0.85, pch = 20)
  text(77, max(datos), grupo[i], cex = 1)
  if (i %in% texto) {
    mtext(Escenarios[i], 1, line = 3.2, cex = 0.75)
  } else {
    if (i == 4) {
      mtext(expression(paste("Ca(0,", widehat(sigma), ")-", "N(0,1) scenarios"), 
        sep = ""), 1, line = 3.2, cex = 0.75)
    } else {
      mtext(expression(paste("logN(1,1)-Ga(", hat(alpha), ",", hat(beta), ") scenarios"), 
        sep = ""), 1, line = 3.2, cex = 0.75)
    }
  }
}


## Code for output of Table 3 (section 6)
# The same samples of section 5 are analysed
# p-values dgeometric.test with n.sim=1000 saved in object p.values1000
# The object resumen contains the numbers in the table
library("GoFKernel")
set.seed(235)
p.values1000 <- matrix(nrow = 1000, ncol = 48)
tamanys <- rep(c(10, 20, 50, 100, 200, 500), 8)
esc <- rep(c("normal", "cauchy", "uniform", "beta", "exponencial", "gamma", "logN", 
  "logN2"), each = 6)
colnames(p.values1000) <- paste(esc, "_", tamanys, sep = "")

# Normal
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dnorm, n.sim = 1000)$p.value, 
  error = function(e) NA))
p.values1000[, 1] <- apply(normal_10, 1, dgeometric.n)
p.values1000[, 2] <- apply(normal_20, 1, dgeometric.n)
p.values1000[, 3] <- apply(normal_50, 1, dgeometric.n)
p.values1000[, 4] <- apply(normal_100, 1, dgeometric.n)
p.values1000[, 5] <- apply(normal_200, 1, dgeometric.n)
p.values1000[, 6] <- apply(normal_500, 1, dgeometric.n)


# Cauchy
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dnorm, list(mean = 0, 
  sd = sd(x)), n.sim = 1000)$p.value, error = function(e) NA))
p.values1000[, 7] <- apply(cauchy_10, 1, dgeometric.n)
p.values1000[, 8] <- apply(cauchy_20, 1, dgeometric.n)
p.values1000[, 9] <- apply(cauchy_50, 1, dgeometric.n)
p.values1000[, 10] <- apply(cauchy_100, 1, dgeometric.n)
p.values1000[, 11] <- apply(cauchy_200, 1, dgeometric.n)
p.values1000[, 12] <- apply(cauchy_500, 1, dgeometric.n)


# Uniform
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dunif, lower = 0, 
  upper = 1, n.sim = 1000)$p.value, error = function(e) NA))
p.values1000[, 13] <- apply(unif_10, 1, dgeometric.n)
p.values1000[, 14] <- apply(unif_20, 1, dgeometric.n)
p.values1000[, 15] <- apply(unif_50, 1, dgeometric.n)
p.values1000[, 16] <- apply(unif_100, 1, dgeometric.n)
p.values1000[, 17] <- apply(unif_200, 1, dgeometric.n)
p.values1000[, 18] <- apply(unif_500, 1, dgeometric.n)

# beta
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dunif, lower = 0, 
  upper = 1, n.sim = 1000)$p.value, error = function(e) NA))
p.values1000[, 19] <- apply(beta_10, 1, dgeometric.n)
p.values1000[, 20] <- apply(beta_20, 1, dgeometric.n)
p.values1000[, 21] <- apply(beta_50, 1, dgeometric.n)
p.values1000[, 22] <- apply(beta_100, 1, dgeometric.n)
p.values1000[, 23] <- apply(beta_200, 1, dgeometric.n)
p.values1000[, 24] <- apply(beta_500, 1, dgeometric.n)

# Exponencial
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dexp, lower = 0, n.sim = 1000)$p.value, 
  error = function(e) NA))
p.values1000[, 25] <- apply(exp_10, 1, dgeometric.n)
p.values1000[, 26] <- apply(exp_20, 1, dgeometric.n)
p.values1000[, 27] <- apply(exp_50, 1, dgeometric.n)
p.values1000[, 28] <- apply(exp_100, 1, dgeometric.n)
p.values1000[, 29] <- apply(exp_200, 1, dgeometric.n)
p.values1000[, 30] <- apply(exp_500, 1, dgeometric.n)

# gamma
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dexp, lower = 0, n.sim = 1000)$p.value, 
  error = function(e) NA))
p.values1000[, 31] <- apply(gamma_10, 1, dgeometric.n)
p.values1000[, 32] <- apply(gamma_20, 1, dgeometric.n)
p.values1000[, 33] <- apply(gamma_50, 1, dgeometric.n)
p.values1000[, 34] <- apply(gamma_100, 1, dgeometric.n)
p.values1000[, 35] <- apply(gamma_200, 1, dgeometric.n)
p.values1000[, 36] <- apply(gamma_500, 1, dgeometric.n)

# logN
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dlnorm, list(meanlog = 1, 
  sdlog = 1), lower = 0, n.sim = 1000)$p.value, error = function(e) NA))
p.values1000[, 37] <- apply(lnorm_10, 1, dgeometric.n)
p.values1000[, 38] <- apply(lnorm_20, 1, dgeometric.n)
p.values1000[, 39] <- apply(lnorm_50, 1, dgeometric.n)
p.values1000[, 40] <- apply(lnorm_100, 1, dgeometric.n)
p.values1000[, 41] <- apply(lnorm_200, 1, dgeometric.n)
p.values1000[, 42] <- apply(lnorm_500, 1, dgeometric.n)

# logN.gamma
dgeometric.n <- function(x) return(tryCatch(dgeometric.test(x, dgamma, lapply(fitdistr(x, 
  "gamma")$estimate, function(i) i), lower = 0, n.sim = 1000)$p.value, error = function(e) NA))
p.values1000[, 43] <- apply(lnorm_10, 1, dgeometric.n)
p.values1000[, 44] <- apply(lnorm_20, 1, dgeometric.n)
p.values1000[, 45] <- apply(lnorm_50, 1, dgeometric.n)
p.values1000[, 46] <- apply(lnorm_100, 1, dgeometric.n)
p.values1000[, 47] <- apply(lnorm_200, 1, dgeometric.n)
p.values1000[, 48] <- apply(lnorm_500, 1, dgeometric.n)

seleccionar <- c(9, 18, 27, 36, 45, 54)
p100 <- cbind(normal.normal[, seleccionar], cauchy.normal[, seleccionar], unif.unif[, 
  seleccionar], beta.unif[, seleccionar], exp.exp[, seleccionar], gamma.exp[, seleccionar], 
  lnorm.lnorm[, seleccionar], lnorm.gamma[, seleccionar])
alfa.1 <- (p.values1000 < 0.1) == (p100 < 0.1)
alfa.05 <- (p.values1000 < 0.05) == (p100 < 0.05)
alfa.01 <- (p.values1000 < 0.01) == (p100 < 0.01)
resumen <- matrix(NA, 48, 6)
resumen[, 1] <- apply(alfa.1, 2, sum)/1000
resumen[, 2] <- apply(alfa.05, 2, sum)/1000
resumen[, 3] <- apply(alfa.01, 2, sum)/1000
resumen[, 4] <- diag(cor(p.values1000, p100))
diferencia <- abs(p100 - p.values1000)
resumen[, 5] <- apply(diferencia, 2, mean)
resumen[, 6] <- apply(diferencia, 2, sd)
resumen[is.na(resumen)] <- 1  # scenarios without variability in p.values
resumen <- rbind(cbind(resumen[13:18, ], resumen[19:24, ]), cbind(resumen[1:6, ], 
  resumen[7:12, ]), cbind(resumen[25:30, ], resumen[31:36, ]), cbind(resumen[37:42, 
  ], resumen[43:48, ]))
