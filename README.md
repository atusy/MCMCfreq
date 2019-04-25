# MCMCfreq
## Overviwe
The MCMCfreq package provides Bayesian Gaussian Process regression with Stan for one dimentional data. Althoug this package was created to correctly estimate the age frequency spectrum of detrital zircon, it can be widely applied to various data instead of histogram.

## Installation

In advance, you need a working environment of Stan.

```r
install.packages("rstan")
install.packages("rstantools")
install.packages("devtools")
library(devtools)
install_github("Tan-Furukawa/MCMCfreq")
```

## Examples

The function `auto_cmpt_freq()` is effective for detrital zircon data.

```r
library(MCMCfreq)
d <- Osayama
e <- auto_cmpt_freq(d)
freq_graph(e, hist = T, ylab = "Frequency")
```

The function `cmpt_freq()` is effective for all one dimentional data.
If WAIC from one pair of rho and sigma is smaller than another pair, it means that is better model;
therefore, we have to serch better rho and sigma to make WAIC smaller.


```r
library(MCMCfreq)
d <- Osayama
e <- cmpt_freq(d, rho = 3, sigma = 2)
freq_graph(e, hist = T, ylab = "Frequency")
WAIC <- e$WAIC
print(WAIC)
```
```r
>print(WAIC)
[1] 275.4525
```

## Author
Tan Furukawa 

e-mail: rpackagetan@gmail.com




