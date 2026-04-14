<div align="center">
  <img src="img/logo.png" alt="BISAM Logo" width="200"/>

<div id="user-content-toc">
    <ul align="center" style="list-style: none; margin-bottom: 0;">
        <summary>
        <h1>BISAM</h1>
        <p>Bayesian Indicator Saturation Model in R</p>
        </summary>
    </ul>
</div>

---

  <!-- <a href="https://github.com/Avriox/BISAM">
    <img src="https://img.shields.io/badge/Status-Under%20Development-orange.svg" alt="Development Status"/>
  </a>
  <a href="https://github.com/Avriox/BISAM/actions/workflows/R-CMD-check.yaml">
    <img src="https://github.com/Avriox/BISAM/actions/workflows/R-CMD-check.yaml/badge.svg?branch=candidate" alt="R CMD Check Status"/>
  </a>
  <a href="https://cran.r-project.org/package=BISAM">
    <img src="https://www.r-pkg.org/badges/version/BISAM" alt="CRAN Status"/>
  </a>
  <a href="LICENSE"> <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"/> </a>
  <img src="https://img.shields.io/badge/Depends-Rcpp-blue.svg" alt="Depends Rcpp"/> -->
</div>

> [!WARNING]
> BISAM is currently under **active development**. The API may change without notice, and the package is **not yet available on CRAN**. The most up-to-date development version is available in the `candidate` branch. Please use with caution in production environments.

---

## Overview

**BISAM** provides an R implementation of Bayesian Indicator Saturation Models. This statistical technique is designed to robustly identify and model multiple structural breaks, outliers, or other forms of parameter instability within regression models or time series data.

The core computations are accelerated using **Rcpp** for improved performance, making it suitable for analysing larger datasets or more complex models.

---

## 🚀 Installation

You can install BISAM using one of the following methods:

**1. From GitHub:**

The latest development version (currently recommended) is available on the `candidate` branch. You will need to clone the content of this branch to your local disk.

```R
library(Rcpp)
setwd("path/to/BISAM")
Rcpp::compileAttributes(verbose = TRUE)
install.packages("./", repos = NULL, type = "source", verbose = TRUE, clean=TRUE)
```

**2. From GitHub Releases (Planned):**

Once official releases are tagged, you will be able to download the source package (`BISAM_x.y.z.tar.gz`) from the [Releases Page](https://github.com/Avriox/BISAM/releases).

```R
# Replace 'path/to/BISAM_x.y.z.tar.gz' with the actual path and filename
install.packages("path/to/BISAM_x.y.z.tar.gz", repos = NULL, type = "source")
```

**3. From CRAN (Planned):**

BISAM is planned for submission to the Comprehensive R Archive Network (CRAN) in the future. Once accepted, you will be able to install it easily:

```R
# install.packages("BISAM") # Not yet available!
```

---

## ✨ Quick Start

Here's a basic example of how to use BISAM (replace with your actual primary function 

```R
# Load the package
library(BISAM)

# Generate or load some sample data
# set.seed(123)
# y <- arima.sim(list(order = c(1,0,0), ar = 0.7), n = 100)
# y[50:55] <- y[50:55] + 5 # Introduce a break/outliers
# x <- rnorm(100)
# data <- data.frame(y = y, x = x)

# Fit the BISAM model (replace with actual function and arguments)
# model_fit <- bisam(y ~ x, data = data, ...) # Fictional function call

# Explore the results
# summary(model_fit)
# plot(model_fit)

# Please replace the above placeholder code with a minimal working example
# demonstrating the core functionality of your package.
```

---

## 🛠️ Dependencies

* R (>= Your required R version, e.g., 4.0.0)
* Rcpp
* *(List other key R package dependencies here, e.g., stats, graphics)*

---

## Acknowledgements & Foundational Work

This package builds upon the indicator saturation methodology. Key references include:

* **Original Methodology:** *(Add citation for the paper that introduced the strategy, e.g., Author(s), Year, Title, Journal/Link)*
* **Computational Enhancements:** *(Add citation for your master thesis, e.g., Your Name, Year, Thesis Title, University/Link)*

---

## License

This project is licensed under the **[Your License Name]** License. See the [LICENSE](LICENSE) file for details.

*(Common choices for R packages include MIT, GPL-2, GPL-3. Make sure you add a `LICENSE` file to your repository root containing the full license text).*

---

## Citation (Package)

TODO