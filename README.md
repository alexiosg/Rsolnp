
[![R-CMD-check](https://github.com/alexiosg/Rsolnp/actions/workflows/rcmdcheck.yaml/badge.svg)](https://github.com/alexiosg/Rsolnp/actions/workflows/rcmdcheck.yaml)
[![Last-changedate](https://img.shields.io/badge/last%20change-2025--07--05-yellowgreen.svg)](/commits/master)
[![packageversion](https://img.shields.io/badge/Package%20version-2.0.2-orange.svg?style=flat-square)](commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/Rsolnp)](https://cran.r-project.org/package=Rsolnp)

# Rsolnp

Since version 2.0.0, a new function, `csolnp`, represents an faster
version of solnp based on a full translation to C++. Additionally,
options for passing analytic gradient and Jacobian functions is now
implemented. Other enhancements include a new test suite based on the
Hock-Schittkowski 306 problems.

A new multi-start solver has also been implemented in the `csolnp_ms`
function (see vignette for more details).

## Installation

The package can be installed from CRAN or the github repo:

``` r
remotes::install_github("alexiosg/Rsolnp", dependencies = TRUE)
```
