
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="man/figures/logo.png" align="left" height="139" alt="" /> WData

<!-- badges: start -->

[![R-CMD-check](https://github.com/noeliasanchmrt/WData/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/noeliasanchmrt/WData/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/noeliasanchmrt/WData/graph/badge.svg)](https://app.codecov.io/gh/noeliasanchmrt/WData)
<!-- badges: end -->

Set of tools for analyzing and modeling data that may be subject to
biases in sampling. It offers functions to estimate density function and
cumulative distribution function from a biased sample of a continuous
distribution. Regarding density function estimation, the package
includes Bhattacharyya et al. (1988) and Jones (1991) density estimators
and various bandwidth selectors for the latter, enhancing the
flexibility and adaptability of density estimation to different types of
datasets and biases. For cumulative distribution function estimation,
the package includes the empirical estimator proposed by Cox (2005) and
the kernel-type estimator by Bose and Dutta (2022), along with several
bandwidth selectors for the latter. Finally, the package includes
Muttlak (1988) real length-biased dataset on shrub width as an example
dataset.

## Installation

You can install the development version of WData from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("noeliasanchmrt/WData")
library(WData)
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/tidyverse/dplyr/issues). For questions and
other discussion, please use [forum.posit.co](https://forum.posit.co/).

## Usage

### Real data

The species *Cercocarpus montanus*, commonly known as mountain mahogany,
is a deciduous shrub native to the western United States and northern
Mexico. It is typically found on slopes, canyons, and rocky, arid
formations with calcareous or other alkaline soils. This small shrub has
white flowers and oval-shaped leaves with serrated edges. It is highly
drought-tolerant and can survive even in nutrient-poor soils. Its
intricate root system prevents landslides on sloped terrains, while its
branches provide food and habitat for various animal species. These
characteristics make *Cercocarpus montanus* an ideal species for
studying wildlife recovery in a given geographic region.

During the fall semester of 1986, graduate students in a biological
sampling techniques course taught by Lyman L. McDonald at the University
of Wyoming conducted a study on the size of *Cercocarpus montanus* in an
old limestone quarry located just east of Laramie, Wyoming (United
States).

The sampling was conducted using line transect methods, which are widely
used in ecological studies to measure species abundance in a given area
and other relevant parameters. These methods involve randomly placing
parallel sampling lines (transects) across the study area. Researchers
traverse these lines, recording variables of interest.

To establish the sampling lines, a baseline was set across the study
region. Random positions were generated along this baseline following a
uniform distribution, and each sampling line was drawn perpendicularly
from these points. One limitation of this approach is the need for a
large number of transects to cover the area adequately. An alternative
method involves selecting a single random point along the baseline and
setting a fixed distance between transects. More details on variations
of this method can be found in Buckland et al. (2001).

The Laramie quarry was covered with north-south-oriented rock fissures.
Since moisture levels and vegetation density were higher in these
fissures, a baseline parallel to them was established. The transects
were drawn perpendicular to this baseline, crossing the terrain fissures
instead of running parallel to them. A distance of 41.6 meters was set
between transects, and two independent replicates (I and II) were
obtained, each with three equidistant parallel transects. In total, six
transects were surveyed from the baseline to the eastern boundary of the
quarry. Students traversed the transects, identifying *Cercocarpus
montanus* shrubs intersected by the line. For each shrub, its maximum
height, the number of main branches, and its width (the maximum distance
between two parallel tangent lines to the shrub’s contour along the
transect) were measured.

Since *Cercocarpus montanus* is a rhizomatous species and adjacent
shrubs may be interconnected via their root system, a shrub was defined
as an individual if it had a distinct cluster of stems at the base and
was at least 15 centimeters away from its nearest neighbor. For shrub
clusters, the length of their intersection with the transect was
recorded. More details on the sampling procedure and additional
measurements can be found in Muttlak (1988).

Due to the sampling method, wider shrubs had a higher probability of
being intersected by the transects. Consequently, the recorded shrub
widths represent a sample biased by longitudinal bias, meaning the bias
function is given by $w(x) = x$. The height and branch count
measurements were also subject to bias, although the bias function $w$
is more complex as it depends on the relationship between shrub width
and these respective variables.

``` r
summary(shrub.data)
summary(shrub.data$Width)
```

### Density estimation

#### `df.bhatta()`: Bhattacharyya et al. (1988) density estimator

``` r
library(WData)
par(mfrow = c(1, 3))
bhatta <- df.bhatta(shrub.data$Width, bw = "nrd0", kernel = "gaussian", from = -0.4, to = 3)
bw.ucv <- bw.ucv(shrub.data$Width, lower = 0.15, upper = 0.3)
bhatta <- df.bhatta(shrub.data$Width, bw = bw.ucv, kernel = "gaussian", from = -0.4, to = 3)
bhatta <- df.bhatta(shrub.data$Width, bw = "SJ-ste", kernel = "gaussian", from = -0.3, to = 3)
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-2-1.png" alt="@bhattacharyya1988 density estimator for shrub width." width="100%" />
<p class="caption">

Bhattacharyya et al. (1988) density estimator for shrub width.
</p>

</div>

#### `df.jones()`: Jones (1991) density estimator

The function allows different bandwidth selection methods:

- `"bw.f.BGMnrd0"`: Normal reference rule-of-thumb selector.
- `"bw.f.BGMcv"`: Cross-validation-based selector.
- `"bw.f.BGMboot1"`: Bootstrap-based selector (method 1).
- `"bw.f.BGMboot2"`: Bootstrap-based selector (method 2).

``` r
par(mfrow = c(2, 3))
jones <- df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGMnrd0", from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
jones <- df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGMcv", lower = 0.01, upper = 0.5, nh = 500L, from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
#> Interval where bandwidth is searched: [0.010000, 0.500000]
jones <- df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGMboot1", from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
#> Pilot Bandwidth for Bootstrap: 0.293207
jones <- df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGMboot1", bw0 = "Opt", from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
#> Pilot Bandwidth for Bootstrap: 0.267953
bw.f.BGMboot2 <- bw.f.BGMboot2(y = shrub.data$Width, from = 0.001, to = 3, nh = 200L, plot = F)
#> Interval where bandwidth is searched: [0.000161, 217.341159]
#> Interval where density is evaluated: [0.001000, 3.000000]
#> Pilot Bandwidth for Bootstrap: 0.075912
jones <- df.jones(shrub.data$Width, kernel = "gaussian", bw = bw.f.BGMboot2, from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="@jones1991 density estimator for shrub width." width="100%" />
<p class="caption">

Jones (1991) density estimator for shrub width.
</p>

</div>

### Distribution estimation

#### `cdf.cox()`: Cox (2005) distribution estimator

``` r
par(mfrow = c(1, 1))
plot(cdf.cox(shrub.data$Width), xlab = "", ylab = "", main = "", col = "blue", xlim = c(0, 3))
rug(shrub.data$Width)
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-4-1.png" alt="@cox2005 distribution estimator for shrub width." width="50%" />
<p class="caption">

Cox (2005) distribution estimator for shrub width.
</p>

</div>

#### Bose and Dutta (2022) distribution estimator

##### `bw.F.BD()`: Bose and Dutta (2022) local bandwidth selector

``` r
par(mfrow = c(2, 2))
bd <- cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", c_adj = rep(0.25, 512))
#> Interval for Estimation: [0.000000, 3.000000]
bd <- cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", c_adj = rep(0.5, 512))
#> Interval for Estimation: [0.000000, 3.000000]
bd <- cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", c_adj = rep(1.3, 512))
#> Interval for Estimation: [0.000000, 3.000000]
c_adj <- ifelse(seq(from = 0, to = 3, length.out = 512) <= quantile(shrub.data$Width, 0.05) |
  seq(from = 0, to = 3, length.out = 512) >= quantile(shrub.data$Width, 0.95), 0.5, 1.3)
bd <- cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", c_adj = c_adj)
#> Interval for Estimation: [0.000000, 3.000000]
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-5-1.png" alt="@bose2022 distribution estimator for shrub width using local bandwidth selector." width="100%" />
<p class="caption">

Bose and Dutta (2022) distribution estimator for shrub width using local
bandwidth selector.
</p>

</div>

##### `bw.F.SBCnrd0()`, `bw.F.SBCcv()`and `bw.F.SBCpi()`: Global bandwidth selectors

``` r
par(mfrow = c(1, 3))
bd <- cdf.bd(shrub.data$Width, from = 0, to = 3, correction = "left")
#> Interval for Estimation: [0.000000, 3.000000]
bw_cv <- bw.F.SBCcv(shrub.data$Width, lower = 0.05, upper = 0.2, nh = 250, plot = F)
#> Interval where bandwidth is searched: [0.050000, 0.200000]
bd <- cdf.bd(shrub.data$Width, correction = "left", bw = bw_cv)
#> Interval for Estimation: [0.030000, 3.130000]
bd <- cdf.bd(shrub.data$Width, from = 0, to = 3, correction = "left", bw = "bw.F.SBCpi")
#> Interval for Estimation: [0.000000, 3.000000]
#> Pilot Bandwidth: 0.214138
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-6-1.png" alt="@bose2022 distribution estimator for shrub width using global bandwidths." width="100%" />
<p class="caption">

Bose and Dutta (2022) distribution estimator for shrub width using
global bandwidths.
</p>

</div>

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-bhattacharyya1988" class="csl-entry">

Bhattacharyya, B. B., Franklin, L. A., and Richardson, G. D. (1988), “A
comparison of nonparametric unweighted and length-biased density
estimation of fibres,” *Communications in Statistics - Theory and
Methods*, Taylor & Francis, 17, 3629–3644.

</div>

<div id="ref-bose2022" class="csl-entry">

Bose, A., and Dutta, S. (2022), “Kernel based estimation of the
distribution function for length biased data,” *Metrika*, 85, 269–287.

</div>

<div id="ref-buckland2001" class="csl-entry">

Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L.,
Borchers, D. L., and Thomas, L. (2001), *Introduction to distance
sampling: Estimating abundance of biological populations*, Oxford
University Press.

</div>

<div id="ref-cox2005" class="csl-entry">

Cox, D. (2005), “Some sampling problems in technology,” in *Selected
statistical papers of sir david cox*, eds. D. Hand and A. Herzberg,
Cambridge University Press, pp. 81–92.

</div>

<div id="ref-jones1991" class="csl-entry">

Jones, M. C. (1991), “Kernel density estimation for length biased data,”
*Biometrika*, \[Oxford University Press, Biometrika Trust\], 78,
511–519.

</div>

<div id="ref-muttlak1988" class="csl-entry">

Muttlak, H. A. (1988), “Some aspects of ranked set sampling with size
biased probability of selection,” *ProQuest Dissertations and Theses*,
PhD thesis, University of Wyoming.

</div>

</div>
