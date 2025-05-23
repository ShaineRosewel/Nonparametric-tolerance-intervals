
# Basics

**Questions answered by tolerance intervals**  Question (1) leads to a
two-sided interval; questions (2) and (3) lead to one-sided intervals.  
- What interval will contain percent of the population measurements?  
- What interval guarantees that percent of population measurements will
not fall below a lower limit?  
- What interval guarantees that percent of population measurements will
not exceed an upper limit?

## Univariate Tolerance Intervals

### Tolerance factor

**2-sided - manual calculation**

``` r
## Compute the approximate and exact k2 factor.

## Compute the approximate k2 factor for a two-sided tolerance interval. 
## For this example, the standard deviation is computed from the sample,
## so the degrees of freedom are nu = N - 1.
N = 25
nu = N - 1
p = 0.90 # proportion
g = 0.99 # alpha
z2 = (qnorm((1+p)/2))**2
c2 = qchisq(1-g,nu)
k2 = sqrt(nu*(1 + 1/N)*z2/c2)
k2
```

    ## [1] 2.494063

**2-sided tolerance factor thru `tolerance` package**

``` r
## Compute the exact k2 factor for a two-sided tolerance interval using 
## the K.factor function in the tolerance library
library(tolerance)
```

    ## tolerance package, version 3.0.0, Released 2024-04-18
    ## This package is based upon work supported by the Chan Zuckerberg Initiative: Essential Open Source Software for Science (Grant No. 2020-255193).

``` r
K2 = K.factor(n=N, f=nu, alpha=1-g, P=p, side=2, method="EXACT", m=100)
K2
```

    ## [1] 2.505927

### Example

#### The data

``` r
## "Direct" calculation of a tolerance interval.

## Read data and name variables.

# URL of the dataset
#url <- "https://www.itl.nist.gov/div898/handbook/datasets/MPC62.DAT"

# Download the file (you can specify a path where you want to save it)
#download.file(url, destfile = "MPC62.DAT")

# Load the dataset into R
mdat = read.table("MPC62.DAT",header=FALSE, skip=50)
colnames(mdat) = c("cr", "wafer", "mo", "day", "h", "min", "op", 
                 "hum", "probe", "temp", "y", "sw", "df")
```

#### Parametric

**`tolerance` package - parametric**

``` r
## Attach tolerance library and call function.
library(tolerance)
pout=normtol.int(mdat$y, alpha=0.01, P=.90, side=2, method = "EXACT",m=100)
pout
```

    ##   alpha   P    x.bar 2-sided.lower 2-sided.upper
    ## 1  0.01 0.9 97.06984      97.00269      97.13699

**Manually calculated tolerance factor - parametric**

``` r
pout$x.bar-K2*sd(mdat$y); pout$x.bar+K2*sd(mdat$y)
```

    ## [1] 97.00269

    ## [1] 97.13699

#### Nonparametric

**`tolerance` package - nonparametric**

``` r
npout = nptol.int(mdat$y, alpha = 0.01, P = 0.90, side = 2, method = "WILKS",
                  upper = NULL, lower = NULL)
npout
```

    ##   alpha   P 2-sided.lower 2-sided.upper
    ## 1  0.01 0.9        97.014        97.114

``` r
npout
```

    ##   alpha   P 2-sided.lower 2-sided.upper
    ## 1  0.01 0.9        97.014        97.114

**Manually calculated tolerance factor - nonparametric**

``` r
#1 <= r < s < 25


#let s = n−r+1
# n-r+1-r = m 
# (n-m+1)/2 = r

# Parameters
n <- 25    # Number of trials
gamma <- 0.90 # Probability of success in each trial
alpha <- 0.01 # Significance level (1 - confidence level)

# Find the smallest m such that P(Y <= m-1) >= 1 - alpha
m <- qbinom(1 - alpha, size = n, prob = gamma)


r = (n-m+1)/2
s = m+r
```

``` r
sort(mdat$y)[ceiling(r)]; sort(mdat$y)[floor(n-r+1)]
```

    ## [1] 97.014

    ## [1] 97.114

#### Plots

``` r
plottol(pout, mdat$y, plot.type = "both", side = "two", x.lab = "X")
```

![](initial_notes_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
plottol(npout, mdat$y, plot.type = "both", side = "two", x.lab = "X")
```

![](initial_notes_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

-   Wilk’s approach: 2 sided interval

    -   sets values $L(\mathcal{X}) = X_{(r)}$ and
        $U(\mathcal{X}) = X_{(s)}$ where $r < s$, as the limits of the
        $\gamma$-content $100 \left( 1 - \alpha\right)\%$-confidence
        tolerance interval
        $\left( L(\mathcal{X}), U(\mathcal{X}) \right)$
    -   The values of $r$ and $s$ are chosen to satisfy
        $1 \le r \le s \le n$ and $s-r=m$ where $m$ is the smallest
        value for which $P(Y \le m-1) \ge 1-\alpha$, and
        $Y \sim \text{Binomial}(n, \gamma)$. It is customary to take the
        value of $s$ to be equal to $n-r+1$, as suggested by Wilks
        et.al. so that
    -   the nonparametric interval is given by $(X_{(r)}, X_{(n-r+1)})$

[Ref](https://www.itl.nist.gov/div898/handbook/prc/section2/prc263.htm)

#### Aside

``` r
set.seed(5)
simulated <- rlnorm(1000, 0, 1)

hist(simulated)
```

![](initial_notes_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
find_silver_bw <- function(sample){
  n <- length(sample)
  iqr <- IQR(sample)
  sd <- sd(sample)
  return(0.9*n^(-1/5)*min(c(sd, iqr/1.34)))
}

get_grid <- function(sample){
  lb <- min(sort(sample)) - sd(sample)
  ub <- max(sort(sample)) + sd(sample)
  return(seq(from=lb, to=ub, by=0.0001))
}

gaussian_kernel <- function(bandwidth, grid, datapoint) {
  return((1/(bandwidth*sqrt(2*pi))) * exp(-0.5*((grid-datapoint)/bandwidth)^2))
#  return((1 / sqrt(2 * pi)) * exp(-0.5 * u^2))
}


# fcn to estimate 
kde_manual <- function(grid, datapoints, bandwidth) {
  n <- length(datapoints)
  sapply(grid, function(grid_val) {
    sum(gaussian_kernel(bandwidth, grid_val, datapoints) / n)
  })
}

kde_manual_presum <- function(grid, datapoints, bandwidth) {
  n <- length(datapoints)
  sapply(grid, function(grid_val) {
    gaussian_kernel(bandwidth, grid_val, datapoints)
  })
}


# KDE for CDF estimation using Gaussian kernel
kde_cdf <- function(grid, datapoints, bandwidth) {
  sapply(grid, function(grid_val) {
    mean(pnorm((grid_val - datapoints) / bandwidth))
  })
}
```

``` r
estimated_cdf_sim <- kde_cdf(get_grid(simulated), simulated, find_silver_bw(simulated))


plot(get_grid(simulated), estimated_cdf_sim, type = "l", lwd = 2, col = "blue",
     xlab = "x", ylab = "Estimated CDF",
     main = "Estimated CDFs Comparison")
```

![](initial_notes_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
estimated_pdf_sim <- kde_manual(get_grid(simulated), simulated, find_silver_bw(simulated))


plot(get_grid(simulated), estimated_pdf_sim, type = "l", lwd = 2, col = "blue",
     xlab = "x", ylab = "Estimated PDF",
     main = "Estimated CDFs Comparison")
```

![](initial_notes_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
library("GoFKernel")
```

    ## Loading required package: KernSmooth

    ## KernSmooth 2.23 loaded
    ## Copyright M. P. Wand 1997-2009

``` r
kde_cdf_value <- function(x, datapoints, bandwidth) {
  mean(pnorm((x - datapoints) / bandwidth))
}

kde_cdf_inverse <- function(probabilities, datapoints, bandwidth, lower = -10, upper = 10) {
  sapply(probabilities, function(p) {
    root_func <- function(x) kde_cdf_value(x, datapoints, bandwidth) - p
    uniroot(root_func, lower = lower, upper = upper)$root
  })
}

kde_cdf_function <- function(datapoints, bandwidth) {
  function(x) mean(pnorm((x - datapoints) / bandwidth))
}




# =================

U_mat <- kde_cdf(simulated, simulated, find_silver_bw(simulated))
#U_mat <- sapply(p_data, function(i) kde_cdf(i, p_data, h))

# Step 4: Compute Y_ij = max{U_ij, 1 - U_ij}
Y_mat <- pmax(U_mat, 1 - U_mat)

#u<-data.frame(original = transformed_data, complement = 1 - transformed_data)
#y1 <- apply(u, 1, max)
npout = nptol.int(Y_mat, alpha = alpha, P = gamma, side = 1, method = "WILKS",
                upper = NULL, lower = NULL)
p <- npout$`1-sided.upper`

f.inv <- inverse(kde_cdf_function(simulated,find_silver_bw(simulated)),lower=-100,upper=100)
  
#f.inv(p);f.inv(1-p); qlnorm(1-p); exp(qnorm(1-p))


pnorm(log(f.inv(p)),0,1) - pnorm(log(f.inv(1-p)),0,1)
```

    ## [1] 0.9274381

``` r
plnorm(f.inv(p),0,1) - plnorm(f.inv(1-p),0,1)
```

    ## [1] 0.9274381

``` r
f.inv <- inverse(function(i)kde_cdf(i, simulated,find_silver_bw(simulated)),lower=-1000,upper=1000)
  
p <- 0.97
f.inv(p); qnorm(p)
```

    ## [1] 6.773139

    ## [1] 1.880794

``` r
f <- function(x) pbeta(x, shape1=2, shape2=3)
f.inv <- inverse(f,lower=0,upper=1)
f.inv(.2)
```

    ## [1] 0.2123161

``` r
qbeta(p=0.2, shape1=2, shape2=3)
```

    ## [1] 0.2123171

## Kernel density estimation

[Reference
1](https://medium.com/analytics-vidhya/kernel-density-estimation-kernel-construction-and-bandwidth-optimization-using-maximum-b1dfce127073)
[Plotting KDE](https://rpubs.com/mcocam12/kdf_byhand)

Kernel - Kernel functions are used to estimate density of random
variables and as weighing function in non-parametric regression.

``` r
x_data <- sort(c(65, 75, 67, 79, 81, 91))
x_data_sd <- sd(x_data)
x_data_iqr <- IQR(x_data)
x_data_n <- length(x_data)
h_silver <- 0.9*x_data_n^(-1/5)*min(c(x_data_sd, x_data_iqr/1.34))


xi <- 65#x_data[1]
h <- 5.5 #h_silver
min_val <- 40#min(x_data)-round(h)
max_val <- 120#max(x_data)+round(h)
x <- seq(from = min_val, to = max_val, by = 1) # grid
```

### manual

``` r
# Estimate PDF
outmat <- t(kde_manual_presum(x, x_data, h))
estimated_density <- kde_manual(x, x_data, h)


# Estimate CDF
estimated_cdf <- kde_cdf(x, x_data, h)
```

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(ggplot2)

rownames(outmat) <- x
colnames(outmat) <- x_data
mat <- cbind(grid = rownames(outmat), outmat)
rownames(outmat) <- 1:nrow(outmat)

df_long <- as.data.frame(mat,as.numeric) %>%
  pivot_longer(cols = colnames(outmat), names_to = "sample", values_to = "k")

df_long$grid <- as.numeric(df_long$grid)
df_long$sample <- as.factor(df_long$sample)
df_long$k <- as.numeric(df_long$k)
df_long <- df_long %>% mutate(norm_kern = k/6)
comp <- df_long %>% group_by(grid) %>% summarize(composite=sum(k)/x_data_n)
ggplot() +
  geom_line(data = df_long, aes(x = grid, y = norm_kern, color = sample, group = sample), alpha = 0.75) +
  geom_line(data = comp, aes(x = grid, y = composite), alpha = 0.5) +
  xlim(min_val, max_val) +
  geom_vline(xintercept = c(65, 67, 75, 79, 81, 91), linetype = "dotted", size = 0.2) +
  theme_minimal() +
  ggtitle("Kernels at different data points")
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](initial_notes_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
#df_long %>% filter(sample==65) %>% select(k) %>% pull() %>% hist()

#estimated_pdf <- function(n, kernel_matrix){
#  return(rowSums(kernel_matrix)/n)
#}
```

``` r
x_data <- sort(c(65, 75, 67, 79, 81, 91))
x_data_sd <- sd(x_data)
x_data_iqr <- IQR(x_data)
x_data_n <- length(x_data)
h_silver <- 0.9*x_data_n^(-1/5)*min(c(x_data_sd, x_data_iqr/1.34))


xi <- 65#x_data[1]
h <- 5.5 #h_silver
min_val <- 40#min(x_data)-round(h)
max_val <- 120#max(x_data)+round(h)
x <- seq(from = min_val, to = max_val, by = 1) # grid
```

### manual

``` r
# Assume: x_grid and pdf_vals are already defined and same length

# Approximate CDF using left Riemann sums
dx <- diff(x)            # Spacing between x points (assumes equally spaced)
dx <- c(dx, tail(dx, 1))      # Pad to same length

cdf_from_pdf <- cumsum(comp[['composite']] * dx)

plot(x, cdf_from_pdf, type = "l", lwd = 2, col = "blue",
     xlab = "x", ylab = "Estimated CDF",
     main = "Estimated CDFs Comparison")

lines(x, estimated_cdf, lwd = 2, col = "red", lty = 2)

# Add empirical CDF for reference
lines(x, ecdf(x_data)(x), col = "darkgreen", lwd = 2, lty = 3)
legend("bottomright", legend = c("Numerical Integration", "Kernel CDF (pnorm)", "Empirical CDF"),
c(1, 2, 3), lwd = 2)
```

![](initial_notes_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
#legend("bottomright", legend = c("From PDF (numerical integration)", "From pnorm (kernel CDF)"),
    #   col = c("blue", "red"), lty = c(1, 2), lwd = 2)
```

``` r
f <- function(x) pnorm(x, mean = 0, sd = 1)
f.inv <- inverse(f,lower=-100,upper=100)
f.inv(.2)
```

    ## [1] -0.8416276

## DGP

``` r
set.seed(1029)
uninorm <- rlnorm(100,0,1)
```

``` r
npout = nptol.int(uninorm, alpha = 0.05, P = 0.95, side = 2, method = "WILKS",
                  upper = NULL, lower = NULL)
npout
```

    ##   alpha    P 2-sided.lower 2-sided.upper
    ## 1  0.05 0.95     0.1768112      5.782058

``` r
# STEP 1
library("compositions")
```

    ## Welcome to compositions, a package for compositional data analysis.
    ## Find an intro with "? compositions"

    ## 
    ## Attaching package: 'compositions'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     anova, cor, cov, dist, var

    ## The following object is masked from 'package:graphics':
    ## 
    ##     segments

    ## The following objects are masked from 'package:base':
    ## 
    ##     %*%, norm, scale, scale.default

``` r
gamma <- 0.95
alpha <- 0.05
n<-500
MyVar <- matrix(c(
1,0.95,0.95,
0.95,1,0.95,
0.95,0.95,1),byrow=TRUE,nrow=3)
MyMean <- c(0,0,0)
p <- 3

library("MASS")
```

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
set.seed(1724532)
datas <- rlnorm.rplus(n,MyMean,MyVar)
#datas <- mvrnorm(n,MyMean,MyVar, tol = 1e-06, empirical = FALSE)
colnames(datas) <- c("p_1", "p_2", "p_3")
mat <- cbind(grid = rownames(datas), datas)
rownames(datas) <- sapply(1:n, function(i) paste0("n_", sprintf("%03d", i)))

h <- find_silver_bw(datas[,1])
x <- get_grid(datas[,1])

#p1_outmat <-sapply(datas[,1], function(xi) gaussian_kernel(xi = xi, x = x, h = h))
#dx <- diff(x)
#dx <- c(dx, tail(dx, 1))

estimated_cdf <- kde_cdf(x, datas[,1], h)

#cdf_vals <- cumsum(estimated_pdf(n ,p1_outmat) * dx) #ESTIMATED CDF

plot(x, estimated_cdf, type = "l", lwd = 2,
     xlab = "x", ylab = "Estimated CDF",
     main = "Estimated CDF")
```

![](initial_notes_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
plot(x, estimated_cdf, type = "p",
     xlab = "x", ylab = "Estimated CDF",
     main = "Estimated CDF")
```

![](initial_notes_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

``` r
estimated_density <- kde_manual(x, datas[,1], h)
plot(x, estimated_density, type = "l", lwd = 2,
      xlab = "x", ylab = "Estimated CDF",
      main = "Estimated PDF")
```

![](initial_notes_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
hist(estimated_density)
```

![](initial_notes_files/figure-gfm/unnamed-chunk-29-2.png)<!-- -->

``` r
get_upper_limit <- function(p_data){
  h <- find_silver_bw(p_data)
  #x <- get_grid(p_data)
  #estimated_cdf <- kde_cdf(x, p_data, h)
  #transformed_data <- transform_data(p_data, estimated_cdf, x)
  

  U_mat <- kde_cdf(p_data, p_data, h)
  #U_mat <- sapply(p_data, function(i) kde_cdf(i, p_data, h))
  
  # Step 4: Compute Y_ij = max{U_ij, 1 - U_ij}
  Y_mat <- pmax(U_mat, 1 - U_mat)
  
  #u<-data.frame(original = transformed_data, complement = 1 - transformed_data)
  #y1 <- apply(u, 1, max)
  npout = nptol.int(Y_mat, alpha = alpha, P = gamma, side = 1, method = "WILKS",
                  upper = NULL, lower = NULL)
  return(npout$`1-sided.upper`)
}

get_upper_limit(datas[,1])
```

    ## [1] 0.9661319

``` r
# step 5
k_j <- apply(datas, 2, get_upper_limit)
```

``` r
k_j
```

    ##       p_1       p_2       p_3 
    ## 0.9661319 0.9661984 0.9677460

``` r
# step 6
max(k_j)
```

    ## [1] 0.967746

``` r
# step 7

#p_data <- datas#[,1]

# Custom CDF function

get_tolerance_interval <- function(p_data, maxk=max(k_j)){
  h <- find_silver_bw(p_data)
  
  #f.inv <- inverse(function(i) kde_cdf(i, p_data, h),#kde_cdf(i, p_data, h),
  #                 lower=-100,upper=100)
  f.inv <- inverse(kde_cdf_function(p_data,h),lower=-100,upper=100)
  lower_bound <- f.inv(1 - maxk)
  upper_bound <- f.inv(maxk)

  #print(f.inv(1 - maxk))
  return(c(lower_bound, upper_bound))
}
```

# SIM

``` r
# STEP 1
library("compositions")
gamma <- 0.95
alpha <- 0.05
n<-500
MyVar <- matrix(c(
1,0.95,0.95,
0.95,1,0.95,
0.95,0.95,1),byrow=TRUE,nrow=3)
MyMean <- c(0,0,0)
p <- 3

set.seed(7779)
datas <- rlnorm.rplus(n,MyMean,MyVar)


hist(datas[,1], col = rgb(1, 0, 0, 0.4), xlim = range(datas), main = "Overlaid Histograms", xlab = "Values")
hist(datas[,2], col = rgb(0, 1, 0, 0.4), add = TRUE)
hist(datas[,3], col = rgb(0, 0, 1, 0.4), add = TRUE)
legend("topright", legend = c("1", "2", "3"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 1, 0, 0.4), rgb(0, 0, 1, 0.4)))
```

![](initial_notes_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

``` r
colnames(datas) <- c("p_1", "p_2", "p_3")
mat <- cbind(grid = rownames(datas), datas)
rownames(datas) <- sapply(1:n, function(i) paste0("n_", sprintf("%03d", i)))
```

``` r
par(mfrow = c(1, 3))
hist(datas[,1], main = "1", xlab = "Values", col = "skyblue")
hist(datas[,2], main = "2", xlab = "Values", col = "salmon")
hist(datas[,3], main = "3", xlab = "Values", col = "lightgreen")
```

![](initial_notes_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

``` r
par(mfrow = c(1, 1))
```

``` r
k_j <- apply(datas, 2, get_upper_limit)
k_j
```

    ##       p_1       p_2       p_3 
    ## 0.9672254 0.9669883 0.9666109

``` r
maxkj <- max(k_j)
max(k_j)
```

    ## [1] 0.9672254

``` r
tols <- apply(datas, 2, function(s) get_tolerance_interval(s,maxk=maxkj))
tols
```

    ##              p_1         p_2         p_3
    ## [1,] -0.02632234 -0.05404223 -0.03109261
    ## [2,]  6.34726861  6.83613998  7.22415607

``` r
p_data <- datas[,1]
h <- find_silver_bw(p_data)

#f.inv <- inverse(function(i) kde_cdf(i, p_data, h),#kde_cdf(i, p_data, h),
#                 lower=-100,upper=100)
f.inv <- inverse(kde_cdf_function(p_data,h),lower=-100,upper=100)
lower_bound <- f.inv(1 - maxkj)
upper_bound <- f.inv(maxkj)
lower_bound;upper_bound 
```

    ## [1] -0.02632234

    ## [1] 6.347269

``` r
step3 <- function(i){
  return(plnorm(tols[2,i],MyMean[i],MyVar[i,i]) - plnorm(tols[1,i],MyMean[i],MyVar[i,i]))
}

step4 <- function(i){
  return(tols[2,i] - tols[1,i])
}


I <- min(sapply(1:3, step3)) >= gamma
I
```

    ## [1] TRUE

``` r
Lj <- sapply(1:3, step4)
```

``` r
# Bonferoni

bonferonni <- function(p_data) {
  
  if (length(p_data) >= 119) {
    to <- nptol.int(p_data, alpha = alpha/3, P = gamma, side = 2, method = "WILKS",
                  upper = NULL, lower = NULL)

    lower_bound <- to$`2-sided.lower`
    upper_bound <- to$`2-sided.upper`
      
  }
  
  else {
    lower_bound <- min(p_data)
    upper_bound <- max(p_data)
  }

  return(c(lower_bound, upper_bound))
}
```

``` r
results <- replicate(5000, {
  datas <- rlnorm.rplus(n, MyMean, MyVar)
  
 # library(MASS)
  
  ## Step 1: Generate multivariate normal samples
  #Y <- mvrnorm(n = n, mu = MyMean, Sigma = MyVar)
  
  ## Step 2: Exponentiate to get multivariate lognormal samples
  #datas <- exp(Y)  
  
  
  colnames(datas) <- c("p_1", "p_2", "p_3")
  mat <- cbind(grid = rownames(datas), datas)
  rownames(datas) <- sapply(1:n, function(i) paste0("n_", sprintf("%03d", i)))
  
  k_j <- apply(datas, 2, get_upper_limit)
  tols <- apply(datas, 2, function(s) get_tolerance_interval(s, maxk = max(k_j)))
  bonf <- apply(datas, 2, bonferonni)
  
  step3 <- function(i){
    vb <- plnorm(tols[2, i], MyMean[i], MyVar[i, i]) - plnorm(tols[1, i], MyMean[i], MyVar[i, i])
    #print(vb)
    return(vb)
  }
  
  step4 <- function(i){
    return(tols[2, i] - tols[1, i])
  }

  step3b <- function(i){
    vb <- plnorm(bonf[2, i], MyMean[i], MyVar[i, i]) - plnorm(bonf[1, i], MyMean[i], MyVar[i, i])
    #print(vb)
    return(vb)
  }
  
  step4b <- function(i){
    return(bonf[2, i] - bonf[1, i])
  }

  
  #print(min(sapply(1:3, step3)))
  I <- min(sapply(1:3, step3)) >= gamma
  #print(I)
  Lj <- sapply(1:3, step4)
  
  Ib <- min(sapply(1:3, step3b)) >= gamma
  #print(I)
  Ljb <- sapply(1:3, step4b)
  
  list(Lj = Lj, I = I, Ljb = Ljb, Ib = Ib)
}, simplify = FALSE)

# If you want to extract all I values and all Lj values in a matrix form:
I_vals <- sapply(results, function(res) res$I)
Lj_vals <- t(sapply(results, function(res) res$Lj))
Ib_vals <- sapply(results, function(res) res$Ib)
Ljb_vals <- t(sapply(results, function(res) res$Ljb))


mean(I_vals)
```

    ## [1] 0.9502

``` r
colMeans(Lj_vals)
```

    ##      p_1      p_2      p_3 
    ## 6.431427 6.435920 6.439542

``` r
mean(Ib_vals)
```

    ## [1] 0.9832

``` r
colMeans(Ljb_vals)
```

    ##      p_1      p_2      p_3 
    ## 9.232833 9.206014 9.233567

``` r
mk <- get_upper_limit(uninorm)
get_tolerance_interval(uninorm,mk)
```

    ## [1] -0.2003739  4.9905292

``` r
hist(log(uninorm))
```

![](initial_notes_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

``` r
nptol.int(uninorm, alpha = 0.05, P = 0.95, side = 2, method = "WILKS",
                  upper = NULL, lower = NULL)
```

    ##   alpha    P 2-sided.lower 2-sided.upper
    ## 1  0.05 0.95     0.1768112      5.782058

``` r
hist((uninorm))
```

![](initial_notes_files/figure-gfm/unnamed-chunk-43-2.png)<!-- -->

``` r
hist(kde_cdf(datas[,1], datas[,1], h), probability = TRUE, 
     main = "Histogram of PIT-transformed values",
     xlab = "PIT values")
abline(h = 1, col = "red", lwd = 2) 
```

![](initial_notes_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

[Probability Integral
Transform](https://matthewfeickert.github.io/Statistics-Notes/notebooks/Introductory/probability-integral-transform.html)

# Nonparametric simultaneous tolerance intervals for small dimensions based on kernel density estimates

## Abstract

-   This study aims to contribute to the existing knowledge on the
    computation of STIs by developing procedures to compute
    nonparametric STIs with accurate coverage.
-   The tolerance interval’s property of containing a specified
    proportion of sampled population values with high degree of
    confidence makes its computation meaningful.
-   Whenever there are several populations, it may be of interest to
    perform simultaneous inference. For this reason, this study proposes
    methods to construct simultaneous tolerance intervals (STIs).
-   Furthermore, since in many cases of practical applications the
    underlying distribution is unknown, the proposed STIs are derived
    under a nonparametric setting.
-   The proposed procedures and algorithms are then assessed through
    performance metrics, such as estimated coverage probabilities and
    expected lengths.
-   The performance of the proposed methodology is also compared with
    the Bonferroni correction approach.
-   The proposed methods show accurate results, with coverage
    probabilities close to the nominal level.
-   The nonparametric STIs computed using the proposed methods are
    generally better than the ones obtained through a
    Bonferroni-corrected approach.
-   A real-life application on the assessment of liver function is
    presented to illustrate the use of the proposed method.

## Definitions

-   **Tolerance intervals**

    -   an interval that contains at least a specified proportion
        $\gamma$ of the population, with a specified degree of
        confidence, $100 \left( 1 - \alpha\right)\%$

    -   Let $\mathcal{X} = \left\{ X_1, X_2, \dots, X_n \right\}$ be a
        random sample from univariate distribution $F(\cdot)$ and
        suppose $X \sim F(\cdot)$, where $X$ is independent of the
        random sample, then a subset of $\mathbb{R}$, say
        $T = T (\mathcal{X})$, computed from the random sample, is a
        $\gamma$-content $100 \left( 1 - \alpha\right)\%$-confidence
        tolerance interval for $X$ if it satisfies the ff: (The
        probability (over all possible samples $\mathcal{X}$) that the
        tolerance interval $T$ constructed from the sample will capture
        at least a proportion $\gamma$ of the population is equal to
        $1-\alpha$)

    -   *If I keep collecting samples and building intervals over and
        over, then 95% of the time, those intervals will contain at
        least 90% of the population.* - That’s a ($\gamma$, $1-\alpha$)
        tolerance interval.

    $$
    P_{\mathcal{X}} \left\{ P_X \left( X \in T \mid \mathcal{X} \right) \geq \gamma \right\} = 1 - \alpha
    $$

    -   set $T$ could be of the form:

        -   $\left( L(\mathcal{X}), U(\mathcal{X}) \right)$ for a
            two-sided tolerance interval
        -   $\left(-\infty, U(\mathcal{X}) \right)$ for a one-sided
            upper tolerance interval
        -   $\left( L(\mathcal{X}), +\infty) \right)$ for a one-sided
            lower tolerance interval
        -   $U(\mathcal{X})$ and $L(\mathcal{X})$: upper and lower
            tolerance limits

    -   sample use case: deriving appropriate process capability limits
        ($U(\mathcal{X})$ and $L(\mathcal{X})$) for a manufactured
        product ($X$), so that with a given level of confidence
        ($1 - \alpha$), they contain the capability measurements of at
        least a specified proportion of units ($\gamma$) from the
        sampled manufacturing process (ie, gusto mo may confidence ka na
        yung tolerance interval mo, contains $\gamma$ proportion nung
        samples)

-   **Simultaneous inference**

    -   Whenever there are several populations, it is of interest to
        construct simultaneous statistical intervals from the sample
        data.
    -   Having several simultaneous interval estimates that are naïvely
        constructed only increases the likelihood that at least one of
        such inferential statements does not hold true.
    -   For this reason, this study considers integrating the criteria
        for simultaneous inference in computing statistical tolerance
        intervals so that the resulting coverage probabilities are close
        to the nominal confidence level.

-   **Current studies**

    -   The STIs that are available in literature rely on the assumption
        of normality. However, in many applications, this assumption is
        unwarranted. For instance, Wright and Royston (1999) point out
        that positive skewness is common among laboratory measurements,
        which are the bases for the construction of reference intervals.
        Thus, there is a need to establish STIs through nonparametric
        means.
    -   While approaches to compute nonparametric tolerance regions are
        available in some existing studies \[see, for instance Young and
        Mathew (2020) and Lucagbo (2021)\], these studies do not use the
        criterion for simultaneous inference.

-   **Nonparametric tolerance intervals**

    -   Suppose $X_{(1)}, X_{(2)}, \dots, X_{(n)}$ are order statistics
        of a random sample from a population with a univariate
        continuous distribution function $F(x)$.

    -   Wilk’s approach: 2 sided interval

        -   sets values $L(\mathcal{X}) = X_{(r)}$ and
            $U(\mathcal{X}) = X_{(s)}$ where $r < s$, as the limits of
            the $\gamma$-content
            $100 \left( 1 - \alpha\right)\%$-confidence tolerance
            interval $\left( L(\mathcal{X}), U(\mathcal{X}) \right)$
        -   The values of $r$ and $s$ are chosen to satisfy
            $1 \le r \le s \le n$ and $s-r=m$ where $m$ is the smallest
            value for which $P(Y \le m-1) \ge 1-\alpha$, and
            $Y \sim \text{Binomial}(n, \gamma)$. It is customary to take
            the value of $s$ to be equal to $n-r+1$, as suggested by
            Wilks et.al. so that
        -   the nonparametric interval is given by
            $(X_{(r)}, X_{(n-r+1)})$

    -   Wilk’s approach: 1 sided interval: Lower tolerance limit

        -   sets $L(\mathcal{X}) = X_{(r)}$ as the lower limit of the
            $\gamma$-content $100 \left( 1 - \alpha\right)\%$-confidence
            tolerance interval.
        -   The value of $r$ is taken to be the largest integer for
            which $P(Y \ge r \mid n, 1- \gamma) \ge 1 - \alpha$ where
            $Y \sim \text{Binomial}(n, 1-\gamma)$.
        -   The non parametric one-sided lower tolerance interval is
            given by $(X_{(r)}, \infty)$

    -   Wilk’s approach: 1 sided interval: Upper tolerance limit

        -   sets $U(\mathcal{X}) = X_{S}$, where it can be shown that
        -   $s=n-r+1$ and $r$ is the value derived for the lower
            tolerance limit.
        -   The nonparametric one-sided upper tolerance interval is
            given by $(-\infty, X_{(n-r+1)})$

    -   There is a minimum sample size requirement in computing the
        nonparametric tolerance intervals. This depends on the values of
        $\gamma$ and $1-\alpha$

    -   The coverage probabilities associated with nonparametric
        tolerance intervals are known to be conservative.

    -   Although in this study we shall only be using Wilks’ approach,
        we mention that recently some developments aimed at improving
        coverage probabilities of nonparametric tolerance intervals have
        been available in literature. For

-   **Kernel density estimation**

    -   The methodologies proposed in this study involve estimating the
        Cumulative Distribution Function (CDF) through kernel density
        estimation (KDE), which is a generalization of density
        estimation through histograms.

    -   A weight function in histogram construction is simply replaced
        by another function, called the kernel function, denoted by
        $K(\cdot)$, which satisfies the ff conditions: $K(\cdot) \ge 0$
        and $\int^{\infty}_{-\infty} K(x) dx = 1$

    -   This function is often taken as a symmetric density function,
        such as a standard normal density function.

    -   KDE, as a nonparametric statistical tool, estimates the unknown
        Probability Density Function (PDF) or CDF.

    -   KDE Methodology

        -   Suppose that $X_1, X_2, \dots, X_n$ is a random sample from
            some unknown PDF $f(x)$ and CDF $F(x)$. Then, KDE estimates
            the PDF as
            $\hat{f}(x) = \frac{1}{n} \sum^n_{i=1} \frac{1}{h} K \left( \frac{x-X_i}{h} \right)$
            where
        -   $h$ is the bandwidth that acts as a smoothing parameter.
            Here, larger values of $h$ result in smoother density
            estimates. As it gets smaller, the density gets rougher.
            Optimal choices for $h$ include Silverman’s rule of thumb
            given by:
            $h = 0.9 n^{-\frac{1}{5}} min \left\{ S, \frac{IQR}{1.34}\right\}$
            where $S$ is the standard deviation of the sample data and
            IQR is its interquartile range.
        -   Common choices for $K(\cdot)$ include Gaussian kernel and
            the Epanichikov kernel.
        -   This study uses the Gaussian kernel and estimtes the CDF by
            plugging the PDF of the standard normal distribution into
            $K(\cdot)$ and integrating the resulting quantity, giving
            $\hat{F}(x)=\frac{1}{n} \sum^n_{i=1} \Phi \left( \frac{x-X_i}{h} \right)$
            where $\hat{F}(x)$ is the estimate of $F(x)$ and $\Phi$ is
            the CDF of the standard normal distribution.

## Methodology

### Data layout & statistical criterion

Here discusses the data to be used in computing STIs and the criterion
to be followed.

**Data**: random sample
$\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n$ of multivariate
measurements with dimension $p$, coming from an unknown continuous
distributions, say $F_X(\cdot)$. Each
$\mathbf{X}_i = \left( X_{i1}, X_{i2}, \dots, X_{ip}\right)'$ is a
$(p \times 1)$ column vector of measurements taken from the $i$th
subject, $i = 1, 2, \dots, n$

#### 2 sided STIs

**Objective**: We want to construct two-sided STIs for each component of
$\mathbf{X} = (X_1, X_2, \dots, X_p)'$. That is, we want to find a
region of the form
$(c_1, d_1) \times (c_2, d_2) \times \dots \times (c_p, d_p)$ where
$c_j$ and $d_j$, $j = 1, 2, \dots, p$, are functions of the random
sample $\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n$ such that

$$
P_{\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n} \left\{ P_{X_j} \left(c_j < X_j <d_j \mid \mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \right) \ge \gamma, \forall j = 1, 2, \dots, p \right\} = 1 - \alpha
$$

The quantities $\gamma$ and $\left( 1-\alpha\right)$ are between 0 and 1
and are referred to as the **content probability** and **confidence
level**.

**Criterion**: With a confidence level of $100(1-\alpha)\%$, each
marginal interval $(c_j, d_j), j = 1, 2, \dots, p$, should have a
content of at least $\gamma$.

#### One sided STI

**Objective**: We want to construct 1-sided STIs.

These only have either lower limits only or upper limits only. That is,
we want to find a region of the forms

$$
(c_1, \infty) \times (c_2, \infty) \times \dots, (c_p, \infty)
$$

and

$$
(-\infty, d_1) \times (-\infty, d_2) \times \dots (-\infty, d_p)
$$

where $c_j$ and $d_j$, $j=1, 2, \dots p$, are still functions of the
random sample $\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n$, such
that

$$
P_{\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n} \left\{ P_{X_j} \left(X_j > c_j \mid \mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \right) \ge \gamma, \forall j = 1, 2, \dots, p \right\} = 1 - \alpha
$$

and

$$
P_{\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n} \left\{ P_{X_j} \left(X_j < d_j \mid \mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \right) \ge \gamma, \forall j = 1, 2, \dots, p \right\} = 1 - \alpha
$$

**Idea**:The process of obtaining simultaneous tolerance limits for the
random vector $X = \left( X_1, X_2, \dots, X_p \right)$ with $p$-variate
distribution $F_X(\cdot)$ can also be thought of as constructing
simultaneous tolerance limits for the $p$ populations represented by the
$p$ components of $\mathbf{X}$.

**Notes**: Methodologies here make no assumptions about the correlation
structure of $\mathbf{X}$ and so procedures also apply to case where
$X_1, X_2, \dots, X_p$ are independent.

### Method A: nonparametric two-sided STIs

Reference: Lucagbo (2021), who develops nonparametric rectangular
tolerance regions. The main idea of the proposed methods is to transform
each component in $\mathbf{X}$.

Let $F_j(\cdot)$ be the CDF of the continuous random variable $X_j$, the
$j$th component of $\mathbf{X}$. Moreover, let
$k \equiv k(\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n)$ and
$k' \equiv k'(\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n)$ be some
functions of the random sample. The $\gamma$-content, $100(1-\alpha)\%$
confidence two-sided STIs for $\mathbf{X}$ are set to be of the
following form:

$$
\prod^p_{j=1} \left( F^{-1}_j (k'), F^{-1}_j(k) \right)
$$

For each $j = 1, 2, \dots, p$, the end points of the $j$th interval are
expressed as quantiles of $X_j$. And aligned with the criterion, the
values of $k$ and $k'$ are computed to satisfy the following condition:

$$
P_{\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n} \left\{ P_{X_j} \left(F_j^{-1} \left(k'\right) < X_j < F_j^{-1} \left(k\right) \mid \mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \right) \ge \gamma, \forall j = 1, 2, \dots, p \right\} = 1 - \alpha
$$

Equivalently,

$$
P_{\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n} \left\{ P_{X_j} \left( k' < F \left( X_j \right) < k \mid \mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \right) \ge \gamma, \forall j = 1, 2, \dots, p \right\} = 1 - \alpha
$$

Note that $F_j(X_j), j = 1, 2, \dots, p$ are identically distributed as
$U(0, 1)$ random variables. Hence, the choice of common $k$ and $k'$ for
the $p$ components.

Moreover, using the symmetry property of the uniform distribution, we
can choose $k'$ as $1-k$ (ex: Since $k' < k$, let k = 0.7. Since dist is
uniform, its max is 1. Hence,
$k' = 1 - k = 1-0.7 = 0.3 \rightarrow k'=0.3 < k=0.7$). Thus, the
criterion can be written as

$$
P_{\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n} \left\{ P_{X_j} \left( 1-k < F_j \left( X_j \right) < k \mid \mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \right) \ge \gamma, \forall j = 1, 2, \dots, p \right\} = 1 - \alpha
\\
\Rightarrow
\\
P_{\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n} \left\{ \min\limits_{1 \le j \le p} P_{X_j} \left( 1-k < F_j \left( X_j \right) < k \mid \mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \right) \ge \gamma \right\} = 1 - \alpha \
\\\text{take the min prob across j for w/c the prob that at least proportion gamma is within the interval}
\\
\text{this is mostly correct/acceptable prob for all j than taking the max one}
$$

Since the CDF $F_j(\cdot), j = 1, 2, \dots p$ are unknown, we estimate
them marginally via KDE, through the procedure in [definitions
section](#definitions)

$$
P_{\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n} \left\{ \min\limits_{1 \le j \le p} P_{X_j} \left( 1-k < \hat{F}_j \left( X_j \right) < k \mid \mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \right) \ge \gamma \right\} = 1 - \alpha
$$

is also equivalent to

$$
P_{\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n} \left\{ \min\limits_{1 \le j \le p} P_{X_j} \left( \max \left\{ \hat{F}_j(X_j) , 1 - \hat{F}_j(X_j)\right\} \left( X_j \right) < k \mid \mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \right) \ge \gamma \right\} = 1 - \alpha
$$

Explanation:

$$
1-k < \hat{F}_j \left( X_j \right) < k
\\
\equiv
\\
1-k < \hat{F}_j\left( X_j \right) \ \text{and} \  \hat{F}_j \left( X_j \right) < k
\\
\equiv
\\
1-\hat{F}_j\left( X_j \right) < k \ \text{and} \  \hat{F}_j \left( X_j \right) < k
\\
\text{this shows that both 1-F and F are less than k and is mathematically equal to saying}
\\
\max \left\{ \hat{F}_j(X_j) , 1 - \hat{F}_j(X_j)\right\} < k
$$

Define $Y_j = \max \left\{ \hat{F}_j(X_j) , 1 - \hat{F}_j(X_j)\right\}$
and disregard the “minimum” condition for now.

Notice $k$ is consistent with the definition of upper tolerance limit of
the random variable $Y_j$ (see: [1-sided STI](#one-sided-sti)). So, if
we focus only on one component $Y_j$, we could estimate $k$ by the upper
tolerance limit of $Y_j$ to be computed from the data using [this
procedure](#definitions) (nonparametric tolerance intervals) by
computing the nonparametric upper tolerance limit based on
$Y_{1j}, Y_{2j}, \dots, Y_{nj}$, where
$Y_{ij} = \max \left\{ \hat{F}_j(X_{ij}) , 1 - \hat{F}_j(X_{ij})\right\}, i=1,2,\dots,n$.

Since $p>1$ and to account for the “minimum” condition (This means that
the weakest (i.e., smallest) probability across all $p$ dimensions
should still meet the threshold $\gamma$), we obtain an estimate of $k$,
say $\hat{k}$, by taking the maximum of the marginal upper tolerance
limits of the $Y_j$s.

$$
\hat{k} = \max \left\{ k_1, k_1, \dots, k_p\right\}
$$

where $k_j$ is a marginal upper tolerance limit for the probability
condition in dimension $j$.

Using this conservative estimate $\hat{k}$ ensures that the minimum of
the inner probabilities as shown below is at least $\gamma$; i.e., in at
least $(1-\alpha)$ fraction of all datasets, the weakest probability
among all $p$ dimensions meets of exceeds $\gamma$

$$
P_{\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n} \left\{ \min\limits_{1 \le j \le p} P_{X_j} \left( \max \left\{ \hat{F}_j(X_j) , 1 - \hat{F}_j(X_j)\right\}  < \hat{k} \mid \mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \right) \ge \gamma \right\} = 1 - \alpha
$$

Now, going back to the original form:

$$
P_{\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n} \left\{ P_{X_j} \left(F_j^{-1} (1-\hat{k}) < X_j < F_j^{-1} (\hat{k}) \mid \mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \right) \ge \gamma, \forall j = 1, 2, \dots, p \right\} = 1 - \alpha
$$

We say, the $\gamma$-content, $100(1-\alpha)\%$ confidence two sided
STIs for the data are given by

$$
\prod^p_{j=1}\left( c_j, d_j \right) = \prod^p_{j=1}\left( \hat{F}_j^{-1}(1-\hat{k}), \hat{F}_j^{-1}(\hat{k}) \right)
$$

------------------------------------------------------------------------

**Algorithm: Nonparametric two-sided STIs**

1.  Data are given by
    $\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n \overset{iid}{\sim}F_{\mathbf{X}}(\mathbf{x})$
    where $\mathbf{X}_i = \left( X_{i1}, X_{i2}, \dots, X_{ip}\right)'$
    is a $(p \times 1)$ vector of measurements taken from the $i$th
    subject, $i = 1, 2, \dots, n$. Also $F_{\mathbf{X}}(\mathbf{x})$ is
    some unknown (continuous) distribution, where $F_j(x)$ is the
    unknown distribution of the $j$th component

2.  For each $j = 1, 2, \dots, p$, obtain the estimate $\hat{F}_j(x)$
    via KDE using $X_{1j}, X_{2j}, \dots, X_{nj}$.

3.  Compute $U_{ij} = \hat{F}_j(X_{ij}), i = 1, 2, \dots, n$ and
    $j = 1, 2, \dots, p$

4.  Compute
    $Y_{ij} = \max \left\{ U_{ij}, 1-U_{ij} \right\}, i = 1, 2, \dots, n$
    and $j = 1, 2, \dots, p$

5.  For each $j=1, 2, \dots, p$, compute the $\gamma$-content,
    $100(1-\alpha)\%$-confidence nonparametric upper tolerance limit
    based on $Y_{1j}, Y_{2j}, \dots, Y_{nj}$ using the methodology
    [here](#definitions) (nonparametric tolerance intervals). Call this
    $k_j, j = 1, 2, \dots, p$

6.  Get the maximum of the $k_j$s. Call this
    $\hat{k} = \max\{k_1, k_2, \dots, k_p\}$

7.  For each $j=1, 2, \dots, p$, compute
    $c_j =\hat{F}_j^{-1}(1-\hat{k})$ and $d_j = \hat{F}_j^{-1}(\hat{k})$

8.  The $\gamma$-content, $100(1-\alpha)\%$-confidence nonparametric
    two-sided STIs are given by $(c_j, d_j), j =1, 2, \dots, p$

common $\hat{k}$ is to be used across all estimated distributions
$\hat{F}_j$

### Simulation A: nonparametric two-sided STIs

This is aimed to evaluate the performance of the proposed methodology.

-   To reflect the potential skewness of the observations, the data are
    generated from a multivariate lognormal distribution with mean
    vector $\mathbf{0}$ and covariance matrix $\Sigma$ in logarithmic
    scale.
-   DGP will be done using “compositions” R package.
-   The confidence level and content are set at $1-\alpha=0.95$ and
    $\gamma = 0.95$.
-   We use the sample sizes $n =100,n=300$, and $n=500$ for the
    simulations.
-   The covariance matrix $\Sigma$ is also varied to represent different
    correlation structures. We run simulations for the values of
    $\Sigma$ of the form $(1-\rho)\mathbf{I}_p + \rho1_p1_p'$. Here
    $\mathbf{I}_p$ is the $(p \times p)$ identity matrix and $1_p$ is a
    $(p \times 1)$ column vector of 1s. This choice of $\Sigma$ is a
    correlation matrix that assumes an exchangeable correlation
    structure. We run simulations for $p=2$ and $p=3$ for the ff
    covariance matrices:

$$
\Sigma_1 = 0.05\mathbf{I}_p + 0.95\ 1_p1'_p =
\begin{bmatrix}
1 & 0.95 & \ldots & 0.95\\
0.95 & 1 & \ldots & 0.95\\
0.95 & 0.95 & \ddots & 0.95\\
\vdots & \vdots & \vdots & \vdots\\
0.95 & 0.95 & \ldots & 1\\
\end{bmatrix}
$$

$$
\Sigma_2 = 0.5\mathbf{I}_p + 0.5\ 1_p1'_p =
\begin{bmatrix}
1 & 0.5 & \ldots & 0.5\\
0.5 & 1 & \ldots & 0.5\\
0.5 & 0.5 & \ddots & 0.5\\
\vdots & \vdots & \vdots & \vdots\\
0.5 & 0.5 & \ldots & 1\\
\end{bmatrix}
$$

$$
\Sigma_3 = 0.8\mathbf{I}_p + 0.2\ 1_p1'_p =
\begin{bmatrix}
1 & 0.2 & \ldots & 0.2\\
0.2 & 1 & \ldots & 0.2\\
0.2 & 0.2 & \ddots & 0.2\\
\vdots & \vdots & \vdots & \vdots\\
0.2 & 0.2 & \ldots & 1\\
\end{bmatrix}
$$

In addition, we also run simulations for the following non-exchangeable
structures $\Sigma$:

$$
\Sigma_4=
\begin{bmatrix}
1 & -0.95 & -0.95 \\
-0.95 & 1 & 0.95 \\
-0.95 & 0.95 & 1
\end{bmatrix}
$$

$$
\Sigma_5=
\begin{bmatrix}
1 & -0.5 & -0.5 \\
-0.5 & 1 & 0.5 \\
-0.5 & 0.5 & 1
\end{bmatrix}
$$

$$
\Sigma_6=
\begin{bmatrix}
1 & -0.2 & -0.2 \\
-0.2 & 1 & 0.2 \\
-0.2 & 0.2 & 1
\end{bmatrix}
$$

$$
\Sigma_7=
\begin{bmatrix}
1 & 0.9 & 0.5 \\
0.9 & 1 & 0.1 \\
0.5 & 0.1 & 1
\end{bmatrix}
$$

$$
\Sigma_8=
\begin{bmatrix}
1 & 0.5 & 0.5 \\
0.5 & 1 & 0.95 \\
0.5 & 0.95 & 1
\end{bmatrix}
$$

$$
\Sigma_9=
\begin{bmatrix}
1 & 0.2 & 0.2 \\
0.2 & 1 & 0.95 \\
0.2 & 0.95 & 1
\end{bmatrix}
$$

$$
\Sigma_{10}=
\begin{bmatrix}
1 & -0.5 & -0.5 \\
-0.5 & 1 & 0.95 \\
-0.5 & 0.95 & 1
\end{bmatrix}
$$

The choices for the covariance matrices $\Sigma_1$ to $\Sigma_{10}$ are
made so that they include as wide variety as possible of the correlation
values in terms of both size and direction.

Aside: An exchangeable variance-covariance matrix is a structured
covariance matrix where all off-diagonal elements are equal, meaning
that each pair of variables has the same covariance.

**Univariate case**

Furthermore, although this study is primarily concerned with
multivariate measurements, we shall also investigate the performance of
the KDE-based method in the univariate case and compare it with the
performance of the standard approach to compute nonparametric tolerance
intervals, which is that of Wilks (1941).

-   For the univariate case, we generate the simulated samples from the
    univariate lognormal distribution with log-scale mean 0.
-   The log-scale variances used are $\sigma^2 = 0.95, \sigma^2=0.5$,
    and $\sigma^2 = 0.20$.

Note that since the proposed methods apply mainly to the multivariate
case, in the univariate care we refrain from saying “proposed method”,
and instead say “KDE-based method” even though the univariate case is a
special case of the multivariate case.

The Gaussian kernel is used in obtaining the KDE. Moreover, the
bandwidth $h$ used in the simulations is Silverman’s rule of thumb,
which is the preferred bandwidth for the Gaussian kernel. Some numerical
simulations (not reported here) have also been performed using a
different bandwidth choice and the results are very similar.

In Step 7 of Algorithm A, the inverse of the estimated distribution
functions are calculated using the `GoFKernel` R package. In addition,
nonparametric tolerance limits described are implemented in R through
the `tolerance` R package. To estimate the coverage probability, this
study uses $M=5000$ simulated samples in running Monte Carlo
simulations.

### Simulation A: Performance Evaluation

1.  Look at the estimated coverage probabilities. The desired or nominal
    confidence level is 0.95.

For method A, also compute the expected lengths of the toleerance
intervals for each component.

**Algorithm: Performance evaluation of the Nonparametric two-sided STIs
obtained from Algo A through estimated coverage probability and expected
lengths**

1.  Generate the random sample
    $\mathbf{X}_1, \mathbf{X}_2, \dots, \mathbf{X}_n$ from the
    multivariate lognormal distribution with log-scale mean
    $\boldsymbol{\mu}$ and given log-scale covariance matrix $\Sigma$.
    Let $\mu_j$ be the $j$th element of $\boldsymbol{\mu}$ and
    $\sigma^2_j$ be the $j$th diagonal element of $\Sigma$.

2.  Compute the 2-sided tolerance limits $c_j$ and $d_j$,
    $j = 1, 2, \dots, p$ using the procedure in Algo A.

3.  Compute
    $\min_{1 \leq j \leq p} \left\{ \Phi \left( \log d_j; \mu_j, \sigma^2_j\right) - \Phi \left( \log c_j; \mu_j, \sigma^2_j\right) \right\}$
    where $\Phi \left( x; \mu_j, \sigma^2_j\right)$ denotes the CDF of
    the $N \left( \mu_j, \sigma_j^2 \right)$

## Basics

### Tolerance intervals

**Parametric and univariate**

Given a sample size $n$, sample mean $\bar{x}$, and sample standard
deviation $s$, the two-sided normal tolerance interval is:

$$
\bar{x} \pm k \cdot s
$$

where $k$ is a factor based on:

$n$, confidence level, $(1-\alpha)$, and population converage $\gamma$

Assumptions:

-   normally distributed data
-   random sampling
-   if unknown mean and variance: We estimate both the mean and the
    standard deviation from the data. Why it matters: That’s what
    distinguishes this from a theoretical tolerance interval based on
    known parameters. It affects the derivation of $\mathcal{k}$
-   moderate sample size: What it means: While larger n is better, a
    sample of 20 is acceptable for parametric TIs, especially if the
    distribution is normal. Why it matters: Small samples lead to more
    variability in s, making the TI wider or less reliable.

**Example**

Case for unknown $\mu$ and $\sigma$: Estimate the mean $\bar{x}$ and
standard deviation, $\sigma$

``` r
#set.seed(777)
#v = rnorm(20, mean = 10, sd = 0.25)

hist(mdat$y)
```

![](initial_notes_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

``` r
# Q-Q plot
qqnorm(mdat$y)
qqline(mdat$y, col = "blue")
```

![](initial_notes_files/figure-gfm/unnamed-chunk-45-2.png)<!-- -->

``` r
shapiro.test(mdat$y)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  mdat$y
    ## W = 0.97084, p-value = 0.6663

``` r
#W: the test statistic
#p-value > 0.05 → data is likely normally distributed
#p-value ≤ 0.05 → data likely not normally distributed
```

Aside:

In R, to find the smallest value $m$ such that
$P(Y \le m-1) \ge 1-\alpha$ for a binomial distribution
$Y \sim \text{Binomial}(n, \gamma)$, you can use the cumulative
distribution function (CDF) for the Binomial distribution, which is
available in the `pbinom()` function.

The condition $P(Y \le m-1) \ge 1-\alpha$ suggests that you are looking
for the quantile $m$ such that the cumulative probability is at least
$1-\alpha$.
