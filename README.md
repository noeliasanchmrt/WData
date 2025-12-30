
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
samples and biases. For cumulative distribution function estimation, the
package includes the empirical estimator proposed by Cox (2005) and the
kernel-type estimator by Bose and Dutta (2022), along with several
bandwidth selectors for the latter. Regarding quantile function, it
includes the estimators proposed by Sen (1984) and the empirical
quantile function based on the distribution function estimator proposed
by Cox (2005).

Finally, the package includes a real length-biased dataset on shrub
width from Muttlak (1988) and fish length data from the ICES *Spanish
North Coast Bottom Trawl Survey* (DATRAS).

## Installation

You can install the development version of WData from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# devtools::install_github("noeliasanchmrt/WData")
library(WData)
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/tidyverse/dplyr/issues). For questions and
other discussion, please use [forum.posit.co](https://forum.posit.co/).

## Usage

### Real data

#### *Cercocarpus montanus* width data

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

#### ICES *Spanish North Coast Bottom Trawl Survey* (DATRAS)

DATRAS is maintained by ICES (International Council for the Exploration
of the Sea) to collate, document, and standardise trawl survey data,
ensure data quality, and facilitate access for stock assessment and fish
community studies. It primarily contains data from bottom trawl surveys
coordinated by ICES expert groups, covering the Baltic Sea, Skagerrak,
Kattegat, North Sea, English Channel, Celtic Sea, Irish Sea, Bay of
Biscay, and the eastern Atlantic from the Shetlands to Gibraltar.

A complete list of ICES DATRAS surveys is available
[here](https://datras.ices.dk/home/surveycode.aspx) and correspondent
data can be manually downloaded from the ICES (2025) [web
portal](https://datras.ices.dk/Data_products/Download/Download_Data_public.aspx)
or obtained using the `icesDatras` R package. Additionally, `icesVocab`
([CRAN](https://cran.r-project.org/web/packages/icesVocab/index.html))
provides an R interface to the [ICES Vocabularies
database](https://vocab.ices.dk/services/), facilitating the mapping of
`ValidAphiaID` (Valid WoRMS AphiaID) to `ScientificName_WoRMS`
(Scientific name according to WoRMS).

The fields available at DATRAS are:

| Field                  | Description                                  |
|------------------------|----------------------------------------------|
| `RecordType`           | Record type identifier (`HL`).               |
| `Survey`               | Survey abbreviation (`SP-NORTH`).            |
| `Quarter`              | Survey quarter.                              |
| `Country`              | Reporting country code (`ES`).               |
| `Ship`                 | Vessel identifier.                           |
| `Gear`                 | Gear type used in the haul.                  |
| `SweepLngt`            | Sweep length (metres).                       |
| `GearEx`               | Gear exploitation / expansion code.          |
| `DoorType`             | Type of trawl door.                          |
| `StNo`                 | National station code or number.             |
| `HaulNo`               | Haul number.                                 |
| `Year`                 | Calendar year (YYYY).                        |
| `SpecCodeType`         | Species coding system used.                  |
| `SpecCode`             | DATRAS species code.                         |
| `SpecVal`              | Species category or value.                   |
| `Sex`                  | Sex of the individual.                       |
| `TotalNo`              | Total number of individuals in the catch.    |
| `CatIdentifier`        | Catch category identifier.                   |
| `NoMeas`               | Number of individuals measured.              |
| `SubFactor`            | Subsampling factor (`SubWgt / CatCatchWgt`). |
| `SubWgt`               | Subsample weight (grams).                    |
| `CatCatchWgt`          | Catch weight per category (grams).           |
| `LngtCode`             | Length measurement code.                     |
| `LngtClass`            | Length (centimetres).                        |
| `HLNoAtLngt`           | Number of individuals at length `LngtClass`. |
| `DevStage`             | Development stage.                           |
| `LenMeasType`          | Length measurement type.                     |
| `ValidAphiaID`         | Valid WoRMS AphiaID.                         |
| `ScientificName_WoRMS` | Scientific name according to WoRMS.          |
| `DateofCalculation`    | Date of last update (YYYYMMDD).              |

Length distributions are recorded for all fish species and selected
commercial cephalopods and shellfish. Length measurement precision
varies by taxon (1 mm to 5 cm) and is tracked in `LngtCode` variable.
More details in the [Manual of the IBTS North Eastern Atlantic
Surveys](https://ices-library.figshare.com/articles/report/SISP_15_-_Manual_of_the_IBTS_North_Eastern_Atlantic_Surveys/19051037).

Fishing gears are selective by species and, within species, by size. The
probability that a fish is retained and hauled to the surface depends on
its body size because of mesh selectivity. Although retention is more
directly determined by girth, length is typically recorded in practice,
and the two measures are strongly linearly related. As a result, the
observed length distributions are biased toward larger individuals
because of size-dependent trawl selectivity (Wileman et al. (1996)).

We have considered *Spanish North Coast Bottom Trawl Survey* from 2024
to ensure recent sampling and then filtered the dataset to identify
species for which all hauls had complete length measurements
(`TotalNo == NoMeas`) and `LngtCode` equal to `.` or `1`, which
indicates a measurement precision higher than 1 cm.

``` r
hldata <- icesDatras::getDATRAS(
  record = "HL",
  survey = "SP-NORTH",
  year = 2024,
  quarter = 1:4
)

hldata$ScientificName_WoRMS <- sapply(
  hldata$Valid_Aphia,
  function(x) icesVocab::getCodeDetail("SpecWoRMS", x)$detail$Description
)

library(dplyr)

res_filt <- hldata %>%
  mutate(
    calc = TotalNo == NoMeas
  ) %>%
  group_by(ScientificName_WoRMS, LngtCode) %>%
  summarise(
    n_haul = n(),
    n_haul_true = sum(calc),
    n_HLNoAtLngt = sum(HLNoAtLngt),
    LngtClass_min = min(LngtClass, na.rm = TRUE),
    LngtClass_max = max(LngtClass, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(perc_true = n_haul_true / n_haul) %>%
  filter(
    perc_true == 1,
    LngtCode %in% c(".", "1"),
    n_haul >= 25,
    n_HLNoAtLngt > 100

  ) %>%
  select(-perc_true, -n_haul_true) %>%
  arrange(-n_HLNoAtLngt, -n_haul)

pander::pander(res_filt,
  caption = "Summary table for 2024 _Spanish North Coast Bottom Trawl Survey_.",
  style = "rmarkdown",
  split.tables = Inf
)
```

| ScientificName_WoRMS | LngtCode | n_haul | n_HLNoAtLngt | LngtClass_min | LngtClass_max |
|:--:|:--:|:--:|:--:|:--:|:--:|
| Conger conger | 1 | 886 | 1392 | 10 | 132 |
| Raja clavata | 1 | 903 | 1280 | 11 | 95 |
| Phycis blennoides | 1 | 190 | 585 | 9 | 51 |
| Eutrigla gurnardus | 1 | 236 | 495 | 9 | 46 |
| Raja montagui | 1 | 261 | 401 | 21 | 106 |
| Lophius budegassa | 1 | 191 | 208 | 6 | 95 |
| Lophius piscatorius | 1 | 167 | 188 | 8 | 102 |
| Chelidonichthys lucerna | 1 | 115 | 154 | 17 | 74 |
| Scomber scombrus | 1 | 75 | 143 | 17 | 37 |
| Leucoraja naevus | 1 | 121 | 138 | 22 | 64 |

Summary table for 2024 *Spanish North Coast Bottom Trawl Survey*.

*Conger conger* was deemed the most appropriate species for length
analyses, as it comprised the highest number of observations among all
species and exhibited the broadest length range.

##### *Conger conger* 2024 data

`spnorthdatras.data` includes *Conger conger* 2024 length data from the
ICES *Spanish North Coast Bottom Trawl Survey* (DATRAS). Missing values
in the original DATRAS files (`-9`) were replaced by `NA`. Some
variables were converted to factors, and several fields from the
original database were omitted; otherwise, the column names and formats
remained aligned with the original specification.

To facilitate the length frequency analysis, we transformed the data
into a long format, where each row represents a single fish measurement.
This was done by expanding rows using the `HLNoAtLngt` column, which
indicates the number of individuals in each length class.

``` r
conger_data_long <- spnorthdatras.data %>%
  select(HaulNo, LngtClass, HLNoAtLngt) %>%
  tidyr::uncount(HLNoAtLngt, .remove = FALSE) %>%
  select(HaulNo, LngtClass) %>%
  arrange(HaulNo, LngtClass)


pander::pander(
  conger_data_long %>%
    summarise(
      n = n(),
      Min = min(LngtClass, na.rm = TRUE),
      Q1 = quantile(LngtClass, 0.25, na.rm = TRUE),
      Median = median(LngtClass, na.rm = TRUE),
      Mean = mean(LngtClass, na.rm = TRUE),
      Q3 = quantile(LngtClass, 0.75, na.rm = TRUE),
      Max = max(LngtClass, na.rm = TRUE)
    ),
  caption = "Summary table for *Conger conger* 2024 data.",
  style = "rmarkdown",
  split.tables = Inf
)
```

|  n   | Min | Q1  | Median | Mean  | Q3  | Max |
|:----:|:---:|:---:|:------:|:-----:|:---:|:---:|
| 1392 | 10  | 35  |   39   | 39.96 | 45  | 132 |

Summary table for *Conger conger* 2024 data.

``` r

conger_data_LngtClass <- conger_data_long %>%
  pull(LngtClass)

plot(cdf.cox(conger_data_LngtClass), col = "blue", main = "", xlab = "Length (cm)", ylab = "Cumulative Distribution Function")
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-spnorthdatras-1.png" alt="Cumulative distribution function for *Conger conger* 2024 length data." width="100%" />
<p class="caption">

Cumulative distribution function for *Conger conger* 2024 length data.
</p>

</div>

### Density estimation

#### `df.bhatta()`: Bhattacharyya et al. (1988) density estimator

``` r
par(mfrow = c(1, 3))
bhatta <- df.bhatta(shrub.data$Width, bw = "nrd0", kernel = "gaussian", from = -0.4, to = 3)
bw.ucv <- bw.ucv(shrub.data$Width, lower = 0.15, upper = 0.3)
bhatta <- df.bhatta(shrub.data$Width, bw = bw.ucv, kernel = "gaussian", from = -0.4, to = 3)
bhatta <- df.bhatta(shrub.data$Width, bw = "SJ-ste", kernel = "gaussian", from = -0.3, to = 3)
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-bhattacharyya1988-1.png" alt="@bhattacharyya1988 density estimator for shrub width." width="100%" />
<p class="caption">

Bhattacharyya et al. (1988) density estimator for shrub width.
</p>

</div>

#### `df.jones()`: Jones (1991) density estimator

The function allows different bandwidth selection methods:

- `"bw.f.BGM.rt"`: Normal reference rule-of-thumb selector.
- `"bw.f.BGM.cv"`: Cross-validation-based selector.
- `"bw.f.BGM.boot1"`: Bootstrap-based selector (method 1).
- `"bw.f.BGM.boot2"`: Bootstrap-based selector (method 2).

``` r
par(mfrow = c(2, 3))
jones <- df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGM.rt", from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
jones <- df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGM.cv", lower = 0.01, upper = 0.5, nh = 100L, from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
#> Interval where bandwidth is searched: [0.010000, 0.500000]
jones <- df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGM.boot1", from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
#> Pilot Bandwidth for Bootstrap: 0.293207
jones <- df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGM.boot1", bw0 = "PI", from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
#> Pilot Bandwidth for Bootstrap: 0.267953
bw.f.BGM.boot2 <- bw.f.BGM.boot2(y = shrub.data$Width, from = 0.001, to = 3, nh = 100L, plot = F)
#> Interval where bandwidth is searched: [0.000161, 217.341159]
#> Interval where density is evaluated: [0.001000, 3.000000]
#> Pilot Bandwidth for Bootstrap: 0.075912
jones <- df.jones(shrub.data$Width, kernel = "gaussian", bw = bw.f.BGM.boot2, from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-jones1991-1.png" alt="@jones1991 density estimator for shrub width." width="100%" />
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

<img src="man/figures/README-cox2005-1.png" alt="@cox2005 distribution estimator for shrub width." width="50%" />
<p class="caption">

Cox (2005) distribution estimator for shrub width.
</p>

</div>

#### Bose and Dutta (2022) distribution estimator

##### `bw.F.BD()`: Bose and Dutta (2022) local bandwidth selector

``` r
par(mfrow = c(2, 2))
bd <- cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", cy.seq = rep(0.25, 512))
#> Interval for Estimation: [0.000000, 3.000000]
bd <- cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", cy.seq = rep(0.5, 512))
#> Interval for Estimation: [0.000000, 3.000000]
bd <- cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", cy.seq = rep(1.3, 512))
#> Interval for Estimation: [0.000000, 3.000000]
cy.seq <- ifelse(seq(from = 0, to = 3, length.out = 512) <= quantile(shrub.data$Width, 0.05) |
  seq(from = 0, to = 3, length.out = 512) >= quantile(shrub.data$Width, 0.95), 0.5, 1.3)
bd <- cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", cy.seq = cy.seq)
#> Interval for Estimation: [0.000000, 3.000000]
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-bose2022_AMSE-1.png" alt="@bose2022 distribution estimator for shrub width using local bandwidth selector." width="100%" />
<p class="caption">

Bose and Dutta (2022) distribution estimator for shrub width using local
bandwidth selector.
</p>

</div>

##### `bw.F.SBC.rt()`, `bw.F.SBC.cv()`and `bw.F.SBC.pi()`: Global bandwidth selectors

``` r
par(mfrow = c(1, 3))
bd <- cdf.bd(shrub.data$Width, from = 0, to = 3, correction = "left")
#> Interval for Estimation: [0.000000, 3.000000]
bw_cv <- bw.F.SBC.cv(shrub.data$Width, lower = 0.05, upper = 0.2, nh = 100, plot = F)
#> Interval where bandwidth is searched: [0.050000, 0.200000]
bd <- cdf.bd(shrub.data$Width, correction = "left", bw = bw_cv)
#> Interval for Estimation: [0.030000, 3.130000]
bd <- cdf.bd(shrub.data$Width, from = 0, to = 3, correction = "left", bw = "bw.F.SBC.pi")
#> Interval for Estimation: [0.000000, 3.000000]
#> Pilot Bandwidth: 0.214138
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-bose2022_AMISE-1.png" alt="@bose2022 distribution estimator for shrub width using global bandwidths." width="100%" />
<p class="caption">

Bose and Dutta (2022) distribution estimator for shrub width using
global bandwidths.
</p>

</div>

### Quantile estimation

#### `qf.sen()`: Sen (1984) quantile estimator

``` r
par(mfrow = c(1, 1))
plot(qf.sen(shrub.data$Width), xlab = "", ylab = "", main = "", col = "blue", xlim = c(0, 1))
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-sen1984-1.png" alt="@sen1984 quantile estimator for shrub width." width="50%" />
<p class="caption">

Sen (1984) quantile estimator for shrub width.
</p>

</div>

#### `qf.SBC()`: Empirical quantile estimator based on Cox (2005)

``` r
par(mfrow = c(1, 1))
plot(qf.SBC(shrub.data$Width), xlab = "", ylab = "", main = "", col = "blue", xlim = c(0, 1))
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-ecdfcox2005-1.png" alt="Empirical quantile estimator for shrub width." width="50%" />
<p class="caption">

Empirical quantile estimator for shrub width.
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

<div id="ref-ices_datras" class="csl-entry">

ICES (2025), “ICES database on trawl surveys (DATRAS),” Copenhagen,
Denmark: <https://datras.ices.dk>.

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

<div id="ref-sen1984" class="csl-entry">

Sen, P. K. (1984), “On asymptotic representations for reduced quantiles
in sampling from a length-biased distribution,” *Calcutta Statistical
Association Bulletin*, 33, 59–68.

</div>

<div id="ref-wileman1996" class="csl-entry">

Wileman, D. A., Ferro, R., Fonteyne, R., and Millar, R. (1996), *Manual
of methods of measuring the selectivity of towed fishing gears*, ICES
cooperative research reports (CRR).

</div>

</div>
