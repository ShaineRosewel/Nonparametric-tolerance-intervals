## Nonparametric simultaneous tolerance intervals

This is related to [Nonparametric simultaneous tolerance intervals for small dimensions based on kernel density estimates](https://www.tandfonline.com/doi/full/10.1080/03610918.2025.2458573). Please visit the associated [webpage](https://shainerosewel.github.io/Nonparametric-tolerance-intervals/initial_notes.html)

### Problems encountered:
- Negative values in the calculated TI due to KDE
- Coding algorithm 18 for 2 sided intervals

### Recent developments:
- Comparison of metrics across different number of variables.

### To do:

- [ ] Bootstrap calibration
- [ ] Persist data replicates for reproducibility or systematically use `set.seed()`
- [ ] Multivariate 1-sided intervals: Try reproducing for other means and variance covariance matrices
- [ ] Rewrite the codes to make them readable
- [ ] Write a separate source code (package, if possible)
- [x] Review lognormal and normal  
- [x] Multivariate 2-sided intervals: Try reproducing for other means and variance covariance matrices
  
### References:
- [Nonparametric simultaneous tolerance intervals for small dimensions based on kernel density estimates](https://www.tandfonline.com/doi/full/10.1080/03610918.2025.2458573)
- [Improving coverage probabilities for parametric tolerance intervals via bootstrap calibration](https://onlinelibrary.wiley.com/doi/10.1002/sim.8537)
- [An Introduction to the Bootstrap](https://www.taylorfrancis.com/books/mono/10.1201/9780429246593/introduction-bootstrap-bradley-efron-tibshirani)
- [Probability Integral Transform](https://matthewfeickert.github.io/Statistics-Notes/notebooks/Introductory/probability-integral-transform.html)
- [Kernel Density Estimation](https://medium.com/analytics-vidhya/kernel-density-estimation-kernel-construction-and-bandwidth-optimization-using-maximum-b1dfce127073)
- [Kernel Density Estimation Applet](https://mathisonian.github.io/kde/)
- [Understanding Gaussian Kernel Density](https://rpubs.com/mcocam12/kdf_byhand)
- [Tolerance intervals for a normal distribution](https://www.itl.nist.gov/div898/handbook/prc/section2/prc263.htm)
