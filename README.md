
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="man/figures/logo.png" align="left" height="139" alt="" /> WData

<!-- badges: start -->

[![R-CMD-check](https://github.com/noeliasanchmrt/WData/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/noeliasanchmrt/WData/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/noeliasanchmrt/WData/graph/badge.svg)](https://app.codecov.io/gh/noeliasanchmrt/WData)
<!-- badges: end -->

Set of tools for analyzing and modeling data that may be subject to
biases in sampling. It offers functions to estimate density functions,
cumulative distribution functions, and quantiles from biased samples of
a continuous distribution. Regarding density function estimation, the
package includes Bhattacharyya et al. (1988) and Jones (1991) density
estimators and various bandwidth selectors for Jones (1991) density
estimator, enhancing the flexibility and adaptability of density
estimation to different types of datasets and biases. For cumulative
distribution function estimation, the package includes the empirical
estimator proposed by Cox (2005) and the kernel-type estimator by Bose
and Dutta (2022), along with several bandwidth selectors for the latter.
For quantile function estimation, it provides two empirical estimators
(one of which corresponds to Sen (1984) proposal) and one kernel-type
estimator. Additionally, for sparsity function estimation, the package
includes the two estimators proposed by Akbari et al. (2019) (one
plug-in and one kernel-type), as well as a third kernel-type estimator.
Moreover, it offers functions for simulating datasets given a bias
function and a reference distribution, providing users with the ability
to generate synthetic datasets for various modeling and analysis
purposes. It also includes functions for constructing boxplots and QQ
plots for biased data, facilitating graphical exploration. Finally, the
package includes Muttlak (1988) real length-biased dataset on shrub
width as an example dataset.

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
df.bhatta(shrub.data$Width, bw = "nrd0", kernel = "gaussian", from = -0.4, to = 3)
#> 
#> Call:
#>  df.bhatta(y = shrub.data$Width, bw = "nrd0", kernel = "gaussian",     from = -0.4, to = 3)
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.2011
#> 
#>        x               y            
#>  Min.   :-0.40   Min.   : 0.000343  
#>  1st Qu.: 0.45   1st Qu.: 0.022069  
#>  Median : 1.30   Median : 0.147644  
#>  Mean   : 1.30   Mean   : 0.432512  
#>  3rd Qu.: 2.15   3rd Qu.: 0.582869  
#>  Max.   : 3.00   Max.   :16.776218  
#>                  NA's   :61
bw.ucv <- bw.ucv(shrub.data$Width, lower = 0.15, upper = 0.3)
df.bhatta(shrub.data$Width, bw = bw.ucv, kernel = "gaussian", from = -0.4, to = 3)
#> 
#> Call:
#>  df.bhatta(y = shrub.data$Width, bw = bw.ucv, kernel = "gaussian",     from = -0.4, to = 3)
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.2696
#> 
#>        x               y           
#>  Min.   :-0.40   Min.   : 0.00085  
#>  1st Qu.: 0.45   1st Qu.: 0.02773  
#>  Median : 1.30   Median : 0.15356  
#>  Mean   : 1.30   Mean   : 0.44577  
#>  3rd Qu.: 2.15   3rd Qu.: 0.55339  
#>  Max.   : 3.00   Max.   :20.26203  
#>                  NA's   :61
df.bhatta(shrub.data$Width, bw = "SJ-ste", kernel = "gaussian", from = -0.3, to = 3)
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-2-1.png" alt="@bhattacharyya1988 density estimator for shrub width." width="100%" />
<p class="caption">
Bhattacharyya et al. (1988) density estimator for shrub width.
</p>

</div>

    #> 
    #> Call:
    #>  df.bhatta(y = shrub.data$Width, bw = "SJ-ste", kernel = "gaussian",     from = -0.3, to = 3)
    #> 
    #> Data: y (89 obs.);   Bandwidth 'bw' = 0.2232
    #> 
    #>        x                y            
    #>  Min.   :-0.300   Min.   : 0.000508  
    #>  1st Qu.: 0.525   1st Qu.: 0.024325  
    #>  Median : 1.350   Median : 0.149939  
    #>  Mean   : 1.350   Mean   : 0.471471  
    #>  3rd Qu.: 2.175   3rd Qu.: 0.575482  
    #>  Max.   : 3.000   Max.   :29.654158  
    #>                   NA's   :47

#### `df.jones()`: Jones (1991) density estimator

The function allows different bandwidth selection methods:

- `"bw.f.BGMnrd0"`: Normal reference rule-of-thumb selector.
- `"bw.f.BGMcv"`: Cross-validation-based selector.
- `"bw.f.BGMboot1"`: Bootstrap-based selector (method 1).
- `"bw.f.BGMboot2"`: Bootstrap-based selector (method 2).

``` r
par(mfrow = c(2, 3))
df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGMnrd0", from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
#> 
#> Call:
#>  df.jones(y = shrub.data$Width, bw = "bw.f.BGMnrd0", kernel = "gaussian",     from = -0.4, to = 3)
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.2268
#> 
#>        x               y            
#>  Min.   :-0.40   Min.   :0.0006323  
#>  1st Qu.: 0.45   1st Qu.:0.0348495  
#>  Median : 1.30   Median :0.1553392  
#>  Mean   : 1.30   Mean   :0.2930535  
#>  3rd Qu.: 2.15   3rd Qu.:0.5231650  
#>  Max.   : 3.00   Max.   :0.9160607
df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGMcv", lower = 0.01, upper = 0.5, nh = 500L, from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
#> Interval where bandwidth is searched: [0.010000, 0.500000]
#> 
#> Call:
#>  df.jones(y = shrub.data$Width, bw = "bw.f.BGMcv", kernel = "gaussian",     from = -0.4, to = 3, lower = 0.01, upper = 0.5, nh = 500L)
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.09543
#> 
#>        x               y            
#>  Min.   :-0.40   Min.   :0.0000001  
#>  1st Qu.: 0.45   1st Qu.:0.0082395  
#>  Median : 1.30   Median :0.1175688  
#>  Mean   : 1.30   Mean   :0.2935432  
#>  3rd Qu.: 2.15   3rd Qu.:0.5142717  
#>  Max.   : 3.00   Max.   :1.3333035
df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGMboot1", from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
#> Pilot Bandwidth for Bootstrap: 0.293074
#> 
#> Call:
#>  df.jones(y = shrub.data$Width, bw = "bw.f.BGMboot1", kernel = "gaussian",     from = -0.4, to = 3)
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.2294
#> 
#>        x               y            
#>  Min.   :-0.40   Min.   :0.0006558  
#>  1st Qu.: 0.45   1st Qu.:0.0356421  
#>  Median : 1.30   Median :0.1571095  
#>  Mean   : 1.30   Mean   :0.2930114  
#>  3rd Qu.: 2.15   3rd Qu.:0.5234226  
#>  Max.   : 3.00   Max.   :0.9117800
df.jones(shrub.data$Width, kernel = "gaussian", bw = "bw.f.BGMboot1", bw0 = "Opt", from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
#> Pilot Bandwidth for Bootstrap: 0.267910
#> 
#> Call:
#>  df.jones(y = shrub.data$Width, bw = "bw.f.BGMboot1", kernel = "gaussian",     from = -0.4, to = 3, bw0 = "Opt")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.2169
#> 
#>        x               y            
#>  Min.   :-0.40   Min.   :0.0005441  
#>  1st Qu.: 0.45   1st Qu.:0.0323594  
#>  Median : 1.30   Median :0.1516872  
#>  Mean   : 1.30   Mean   :0.2931922  
#>  3rd Qu.: 2.15   3rd Qu.:0.5222664  
#>  Max.   : 3.00   Max.   :0.9325705
bw.f.BGMboot2 <- bw.f.BGMboot2(y = shrub.data$Width, from = 0.001, to = 3, nh = 200L, plot = F)
#> Interval where bandwidth is searched: [0.000161, 217.341159]
#> Interval where density is evaluated: [0.001000, 3.000000]
#> Pilot Bandwidth for Bootstrap: 0.075912
df.jones(shrub.data$Width, kernel = "gaussian", bw = bw.f.BGMboot2, from = -0.4, to = 3)
#> Interval for Estimation: [-0.400000, 3.000000]
#> 
#> Call:
#>  df.jones(y = shrub.data$Width, bw = bw.f.BGMboot2, kernel = "gaussian",     from = -0.4, to = 3)
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 1.092
#> 
#>        x               y         
#>  Min.   :-0.40   Min.   :0.0470  
#>  1st Qu.: 0.45   1st Qu.:0.1444  
#>  Median : 1.30   Median :0.2596  
#>  Mean   : 1.30   Mean   :0.2297  
#>  3rd Qu.: 2.15   3rd Qu.:0.3155  
#>  Max.   : 3.00   Max.   :0.3369
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
cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", c_adj = rep(0.25, 512))
#> Interval for Estimation: [0.000000, 3.000000]
#> 
#> Call:
#>  cdf.bd(y = shrub.data$Width, bw = "bw.F.BD", from = 0, to = 3,     correction = "left", c_adj = rep(0.25, 512))
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.14540.14540.11560.10120.0920.085460.080470.076480.073190.07040.067990.065880.064020.062350.060840.059470.058220.057070.0560.055010.054090.053220.052410.051650.050930.050250.04960.048990.04840.047850.047310.046810.046320.045850.04540.044970.044550.044150.043760.043390.043030.042680.042340.042010.04170.041390.041090.04080.040520.040240.039970.039710.039460.039210.038970.038740.038510.038280.038060.037850.037640.037440.037240.037040.036850.036660.036470.036290.036120.035940.035770.035610.035440.035280.035120.034970.034820.034670.034520.034370.034230.034090.033950.033820.033680.033550.033420.03330.033170.033050.032930.032810.032690.032570.032460.032340.032230.032120.032010.031910.03180.03170.031590.031490.031390.031290.03120.03110.0310.030910.030820.030720.030630.030540.030460.030370.030280.03020.030110.030030.029940.029860.029780.02970.029620.029540.029470.029390.029310.029240.029160.029090.029020.028940.028870.02880.028730.028660.028590.028530.028460.028390.028330.028260.02820.028130.028070.0280.027940.027880.027820.027760.02770.027640.027580.027520.027460.02740.027350.027290.027230.027180.027120.027070.027010.026960.02690.026850.02680.026740.026690.026640.026590.026540.026490.026440.026390.026340.026290.026240.026190.026150.02610.026050.0260.025960.025910.025870.025820.025780.025730.025690.025640.02560.025550.025510.025470.025430.025380.025340.02530.025260.025220.025180.025130.025090.025050.025010.024970.024930.024890.024860.024820.024780.024740.02470.024660.024630.024590.024550.024520.024480.024440.024410.024370.024330.02430.024260.024230.024190.024160.024120.024090.024060.024020.023990.023950.023920.023890.023850.023820.023790.023760.023720.023690.023660.023630.02360.023570.023530.02350.023470.023440.023410.023380.023350.023320.023290.023260.023230.02320.023170.023140.023110.023080.023060.023030.0230.022970.022940.022910.022890.022860.022830.02280.022780.022750.022720.022690.022670.022640.022610.022590.022560.022530.022510.022480.022460.022430.022410.022380.022350.022330.02230.022280.022250.022230.02220.022180.022160.022130.022110.022080.022060.022030.022010.021990.021960.021940.021920.021890.021870.021850.021820.02180.021780.021750.021730.021710.021690.021660.021640.021620.02160.021580.021550.021530.021510.021490.021470.021450.021420.02140.021380.021360.021340.021320.02130.021280.021260.021240.021210.021190.021170.021150.021130.021110.021090.021070.021050.021030.021010.020990.020970.020950.020930.020920.02090.020880.020860.020840.020820.02080.020780.020760.020740.020730.020710.020690.020670.020650.020630.020610.02060.020580.020560.020540.020520.020510.020490.020470.020450.020430.020420.02040.020380.020360.020350.020330.020310.02030.020280.020260.020240.020230.020210.020190.020180.020160.020140.020130.020110.020090.020080.020060.020040.020030.020010.020.019980.019960.019950.019930.019920.01990.019880.019870.019850.019840.019820.01980.019790.019770.019760.019740.019730.019710.01970.019680.019670.019650.019640.019620.019610.019590.019580.019560.019550.019530.019520.01950.019490.019470.019460.019440.019430.019420.01940.019390.019370.019360.019340.019330.019320.01930.019290.019270.019260.019250.019230.019220.019210.019190.019180.019160.019150.019140.019120.019110.01910.019080.019070.019060.019040.019030.019020.0190.018990.018980.018960.018950.018940.018930.018910.01890.018890.018870.018860.018850.018840.018820.018810.01880.018780.018770.018760.018750.018730.018720.018710.01870.018690.018670.018660.018650.018640.018620.018610.01860.018590.018580.01856
#> 
#>        x              y           
#>  Min.   :0.00   Min.   :-0.01955  
#>  1st Qu.:0.75   1st Qu.: 0.67113  
#>  Median :1.50   Median : 0.93439  
#>  Mean   :1.50   Mean   : 0.78291  
#>  3rd Qu.:2.25   3rd Qu.: 0.99712  
#>  Max.   :3.00   Max.   : 1.00000
cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", c_adj = rep(0.5, 512))
#> Interval for Estimation: [0.000000, 3.000000]
#> 
#> Call:
#>  cdf.bd(y = shrub.data$Width, bw = "bw.F.BD", from = 0, to = 3,     correction = "left", c_adj = rep(0.5, 512))
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.29070.29070.23130.20230.1840.17090.16090.1530.14640.14080.1360.13180.1280.12470.12170.11890.11640.11410.1120.110.10820.10640.10480.10330.10190.10050.09920.097980.096810.095690.094630.093610.092630.09170.09080.089940.08910.08830.087530.086780.086060.085360.084680.084030.083390.082780.082180.08160.081030.080480.079950.079430.078920.078430.077940.077470.077010.076570.076130.07570.075280.074870.074470.074080.073690.073320.072950.072590.072230.071890.071550.071210.070890.070560.070250.069940.069630.069330.069040.068750.068460.068180.067910.067640.067370.067110.066850.066590.066340.06610.065850.065610.065380.065140.064910.064690.064470.064250.064030.063810.06360.063390.063190.062990.062780.062590.062390.06220.062010.061820.061630.061450.061270.061090.060910.060740.060560.060390.060220.060050.059890.059720.059560.05940.059240.059090.058930.058780.058630.058480.058330.058180.058030.057890.057750.057610.057470.057330.057190.057050.056920.056780.056650.056520.056390.056260.056140.056010.055880.055760.055640.055520.055390.055270.055160.055040.054920.054810.054690.054580.054460.054350.054240.054130.054020.053910.053810.05370.053590.053490.053390.053280.053180.053080.052980.052880.052780.052680.052580.052480.052390.052290.05220.05210.052010.051920.051820.051730.051640.051550.051460.051370.051280.05120.051110.051020.050940.050850.050770.050680.05060.050510.050430.050350.050270.050190.050110.050030.049950.049870.049790.049710.049630.049560.049480.04940.049330.049250.049180.049110.049030.048960.048890.048810.048740.048670.04860.048530.048460.048390.048320.048250.048180.048110.048040.047980.047910.047840.047780.047710.047640.047580.047510.047450.047380.047320.047260.047190.047130.047070.047010.046940.046880.046820.046760.04670.046640.046580.046520.046460.04640.046340.046280.046230.046170.046110.046050.0460.045940.045880.045830.045770.045720.045660.045610.045550.04550.045440.045390.045330.045280.045230.045170.045120.045070.045020.044960.044910.044860.044810.044760.044710.044660.044610.044560.044510.044460.044410.044360.044310.044260.044210.044160.044120.044070.044020.043970.043930.043880.043830.043780.043740.043690.043650.04360.043550.043510.043460.043420.043370.043330.043280.043240.04320.043150.043110.043060.043020.042980.042930.042890.042850.042810.042760.042720.042680.042640.042590.042550.042510.042470.042430.042390.042350.042310.042270.042230.042190.042150.042110.042070.042030.041990.041950.041910.041870.041830.041790.041750.041720.041680.041640.04160.041560.041530.041490.041450.041410.041380.041340.04130.041260.041230.041190.041160.041120.041080.041050.041010.040980.040940.04090.040870.040830.04080.040760.040730.040690.040660.040620.040590.040560.040520.040490.040450.040420.040390.040350.040320.040290.040250.040220.040190.040150.040120.040090.040060.040020.039990.039960.039930.039890.039860.039830.03980.039770.039740.03970.039670.039640.039610.039580.039550.039520.039490.039460.039420.039390.039360.039330.03930.039270.039240.039210.039180.039150.039120.039090.039060.039040.039010.038980.038950.038920.038890.038860.038830.03880.038780.038750.038720.038690.038660.038630.038610.038580.038550.038520.038490.038470.038440.038410.038380.038360.038330.03830.038270.038250.038220.038190.038170.038140.038110.038090.038060.038030.038010.037980.037950.037930.03790.037880.037850.037820.03780.037770.037750.037720.03770.037670.037650.037620.037590.037570.037540.037520.037490.037470.037440.037420.03740.037370.037350.037320.03730.037270.037250.037220.03720.037180.037150.03713
#> 
#>        x              y           
#>  Min.   :0.00   Min.   :-0.04929  
#>  1st Qu.:0.75   1st Qu.: 0.64373  
#>  Median :1.50   Median : 0.92950  
#>  Mean   :1.50   Mean   : 0.76847  
#>  3rd Qu.:2.25   3rd Qu.: 0.99690  
#>  Max.   :3.00   Max.   : 1.00000
cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", c_adj = rep(1.3, 512))
#> Interval for Estimation: [0.000000, 3.000000]
#> 
#> Call:
#>  cdf.bd(y = shrub.data$Width, bw = "bw.F.BD", from = 0, to = 3,     correction = "left", c_adj = rep(1.3, 512))
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.75590.75590.60130.5260.47840.44440.41850.39770.38060.36610.35350.34260.33290.32420.31640.30930.30280.29680.29120.28610.28130.27680.27260.26860.26480.26130.25790.25470.25170.24880.2460.24340.24090.23840.23610.23380.23170.22960.22760.22560.22380.22190.22020.21850.21680.21520.21370.21220.21070.20930.20790.20650.20520.20390.20270.20140.20020.19910.19790.19680.19570.19470.19360.19260.19160.19060.18970.18870.18780.18690.1860.18520.18430.18350.18260.18180.1810.18030.17950.17870.1780.17730.17660.17590.17520.17450.17380.17310.17250.17180.17120.17060.170.16940.16880.16820.16760.1670.16650.16590.16540.16480.16430.16380.16320.16270.16220.16170.16120.16070.16020.15980.15930.15880.15840.15790.15750.1570.15660.15610.15570.15530.15490.15440.1540.15360.15320.15280.15240.1520.15170.15130.15090.15050.15010.14980.14940.1490.14870.14830.1480.14760.14730.1470.14660.14630.1460.14560.14530.1450.14470.14430.1440.14370.14340.14310.14280.14250.14220.14190.14160.14130.1410.14070.14050.14020.13990.13960.13930.13910.13880.13850.13830.1380.13770.13750.13720.1370.13670.13650.13620.1360.13570.13550.13520.1350.13470.13450.13430.1340.13380.13360.13330.13310.13290.13270.13240.13220.1320.13180.13160.13130.13110.13090.13070.13050.13030.13010.12990.12970.12950.12930.1290.12880.12870.12850.12830.12810.12790.12770.12750.12730.12710.12690.12670.12650.12640.12620.1260.12580.12560.12540.12530.12510.12490.12470.12460.12440.12420.1240.12390.12370.12350.12340.12320.1230.12290.12270.12250.12240.12220.12210.12190.12170.12160.12140.12130.12110.1210.12080.12060.12050.12030.12020.120.11990.11970.11960.11940.11930.11920.1190.11890.11870.11860.11840.11830.11810.1180.11790.11770.11760.11750.11730.11720.1170.11690.11680.11660.11650.11640.11620.11610.1160.11580.11570.11560.11550.11530.11520.11510.1150.11480.11470.11460.11450.11430.11420.11410.1140.11380.11370.11360.11350.11340.11320.11310.1130.11290.11280.11270.11250.11240.11230.11220.11210.1120.11190.11170.11160.11150.11140.11130.11120.11110.1110.11090.11070.11060.11050.11040.11030.11020.11010.110.10990.10980.10970.10960.10950.10940.10930.10920.10910.1090.10890.10880.10870.10860.10850.10840.10830.10820.10810.1080.10790.10780.10770.10760.10750.10740.10730.10720.10710.1070.10690.10680.10670.10660.10650.10640.10640.10630.10620.10610.1060.10590.10580.10570.10560.10550.10540.10540.10530.10520.10510.1050.10490.10480.10470.10470.10460.10450.10440.10430.10420.10410.10410.1040.10390.10380.10370.10360.10360.10350.10340.10330.10320.10310.10310.1030.10290.10280.10270.10270.10260.10250.10240.10230.10230.10220.10210.1020.1020.10190.10180.10170.10160.10160.10150.10140.10130.10130.10120.10110.1010.1010.10090.10080.10070.10070.10060.10050.10040.10040.10030.10020.10020.10010.10.099940.099870.09980.099730.099650.099580.099510.099440.099370.09930.099230.099160.099090.099020.098960.098890.098820.098750.098680.098610.098550.098480.098410.098340.098280.098210.098140.098080.098010.097940.097880.097810.097750.097680.097620.097550.097490.097420.097360.097290.097230.097160.09710.097040.096970.096910.096850.096780.096720.096660.096590.09653
#> 
#>        x              y          
#>  Min.   :0.00   Min.   :-0.1175  
#>  1st Qu.:0.75   1st Qu.: 0.5637  
#>  Median :1.50   Median : 0.9117  
#>  Mean   :1.50   Mean   : 0.7284  
#>  3rd Qu.:2.25   3rd Qu.: 0.9956  
#>  Max.   :3.00   Max.   : 1.0000
c_adj <- ifelse(seq(from = 0, to = 3, length.out = 512) <= quantile(shrub.data$Width, 0.05) |
  seq(from = 0, to = 3, length.out = 512) >= quantile(shrub.data$Width, 0.95), 0.5, 1.3)
cdf.bd(shrub.data$Width, correction = "left", from = 0, to = 3, bw = "bw.F.BD", c_adj = c_adj)
#> Interval for Estimation: [0.000000, 3.000000]
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-5-1.png" alt="@bose2022 distribution estimator for shrub width using local bandwiths." width="100%" />
<p class="caption">
Bose and Dutta (2022) distribution estimator for shrub width using local
bandwiths.
</p>

</div>

    #> 
    #> Call:
    #>  cdf.bd(y = shrub.data$Width, bw = "bw.F.BD", from = 0, to = 3,     correction = "left", c_adj = c_adj)
    #> 
    #> Data: y (89 obs.);   Bandwidth 'bw' = 0.29070.29070.23130.20230.1840.17090.16090.1530.14640.14080.1360.13180.1280.12470.12170.11890.11640.11410.1120.110.10820.10640.10480.10330.10190.10050.09920.097980.096810.095690.094630.093610.092630.09170.09080.089940.08910.22960.22760.22560.22380.22190.22020.21850.21680.21520.21370.21220.21070.20930.20790.20650.20520.20390.20270.20140.20020.19910.19790.19680.19570.19470.19360.19260.19160.19060.18970.18870.18780.18690.1860.18520.18430.18350.18260.18180.1810.18030.17950.17870.1780.17730.17660.17590.17520.17450.17380.17310.17250.17180.17120.17060.170.16940.16880.16820.16760.1670.16650.16590.16540.16480.16430.16380.16320.16270.16220.16170.16120.16070.16020.15980.15930.15880.15840.15790.15750.1570.15660.15610.15570.15530.15490.15440.1540.15360.15320.15280.15240.1520.15170.15130.15090.15050.15010.14980.14940.1490.14870.14830.1480.14760.14730.1470.14660.14630.1460.14560.14530.1450.14470.14430.1440.14370.14340.14310.14280.14250.14220.14190.14160.14130.1410.14070.14050.14020.13990.13960.13930.13910.13880.13850.13830.1380.13770.13750.13720.1370.13670.13650.13620.1360.13570.13550.13520.1350.13470.13450.13430.1340.13380.13360.13330.13310.13290.13270.13240.13220.1320.13180.13160.13130.13110.13090.13070.13050.13030.13010.12990.12970.12950.12930.1290.12880.12870.12850.12830.12810.12790.12770.12750.12730.12710.12690.12670.12650.12640.12620.1260.12580.12560.12540.12530.12510.12490.12470.12460.12440.12420.1240.12390.12370.12350.12340.12320.1230.12290.12270.12250.12240.12220.12210.12190.12170.12160.12140.12130.12110.1210.12080.12060.12050.12030.12020.120.11990.11970.11960.11940.11930.11920.1190.11890.11870.11860.11840.11830.11810.1180.11790.11770.11760.11750.11730.11720.1170.11690.11680.11660.11650.11640.11620.11610.1160.11580.11570.11560.11550.11530.11520.11510.1150.11480.11470.11460.11450.11430.11420.11410.1140.11380.11370.11360.11350.11340.11320.11310.1130.11290.11280.11270.11250.11240.11230.11220.11210.1120.11190.11170.11160.042890.042850.042810.042760.042720.042680.042640.042590.042550.042510.042470.042430.042390.042350.042310.042270.042230.042190.042150.042110.042070.042030.041990.041950.041910.041870.041830.041790.041750.041720.041680.041640.04160.041560.041530.041490.041450.041410.041380.041340.04130.041260.041230.041190.041160.041120.041080.041050.041010.040980.040940.04090.040870.040830.04080.040760.040730.040690.040660.040620.040590.040560.040520.040490.040450.040420.040390.040350.040320.040290.040250.040220.040190.040150.040120.040090.040060.040020.039990.039960.039930.039890.039860.039830.03980.039770.039740.03970.039670.039640.039610.039580.039550.039520.039490.039460.039420.039390.039360.039330.03930.039270.039240.039210.039180.039150.039120.039090.039060.039040.039010.038980.038950.038920.038890.038860.038830.03880.038780.038750.038720.038690.038660.038630.038610.038580.038550.038520.038490.038470.038440.038410.038380.038360.038330.03830.038270.038250.038220.038190.038170.038140.038110.038090.038060.038030.038010.037980.037950.037930.03790.037880.037850.037820.03780.037770.037750.037720.03770.037670.037650.037620.037590.037570.037540.037520.037490.037470.037440.037420.03740.037370.037350.037320.03730.037270.037250.037220.03720.037180.037150.03713
    #> 
    #>        x              y           
    #>  Min.   :0.00   Min.   :-0.04929  
    #>  1st Qu.:0.75   1st Qu.: 0.63255  
    #>  Median :1.50   Median : 0.92561  
    #>  Mean   :1.50   Mean   : 0.76489  
    #>  3rd Qu.:2.25   3rd Qu.: 0.99690  
    #>  Max.   :3.00   Max.   : 1.00000

##### `bw.F.SBCnrd0()`, `bw.F.SBCcv()`and `bw.F.SBCpi()`: Global bandwidth selectors

``` r
par(mfrow = c(1, 3))
cdf.bd(shrub.data$Width, from = 0, to = 3, correction = "left")
#> Interval for Estimation: [0.000000, 3.000000]
#> 
#> Call:
#>  cdf.bd(y = shrub.data$Width, from = 0, to = 3, correction = "left")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.2031
#> 
#>        x              y         
#>  Min.   :0.00   Min.   :0.0000  
#>  1st Qu.:0.75   1st Qu.:0.6393  
#>  Median :1.50   Median :0.9236  
#>  Mean   :1.50   Mean   :0.7733  
#>  3rd Qu.:2.25   3rd Qu.:0.9943  
#>  Max.   :3.00   Max.   :1.0000
bw_cv <- bw.F.SBCcv(shrub.data$Width, lower = 0.05, upper = 0.2, nh = 250, plot = F)
#> Interval where bandwidth is searched: [0.050000, 0.200000]
cdf.bd(shrub.data$Width, correction = "left", bw = bw_cv)
#> Interval for Estimation: [0.080000, 2.960000]
#> 
#> Call:
#>  cdf.bd(y = shrub.data$Width, bw = bw_cv, correction = "left")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.1199
#> 
#>        x              y         
#>  Min.   :0.08   Min.   :0.0000  
#>  1st Qu.:0.80   1st Qu.:0.6839  
#>  Median :1.52   Median :0.9319  
#>  Mean   :1.52   Mean   :0.7944  
#>  3rd Qu.:2.24   3rd Qu.:0.9961  
#>  Max.   :2.96   Max.   :1.0000
cdf.bd(shrub.data$Width, from = 0, to = 3, correction = "left", bw = "bw.F.SBCpi")
#> Interval for Estimation: [0.000000, 3.000000]
#> Pilot Bandwidth: 0.214089
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-6-1.png" alt="@bose2022 distribution estimator for shrub width using global bandwiths." width="100%" />
<p class="caption">
Bose and Dutta (2022) distribution estimator for shrub width using
global bandwiths.
</p>

</div>

    #> 
    #> Call:
    #>  cdf.bd(y = shrub.data$Width, bw = "bw.F.SBCpi", from = 0, to = 3,     correction = "left")
    #> 
    #> Data: y (89 obs.);   Bandwidth 'bw' = 0.2066
    #> 
    #>        x              y         
    #>  Min.   :0.00   Min.   :0.0000  
    #>  1st Qu.:0.75   1st Qu.:0.6380  
    #>  Median :1.50   Median :0.9232  
    #>  Mean   :1.50   Mean   :0.7728  
    #>  3rd Qu.:2.25   3rd Qu.:0.9942  
    #>  Max.   :3.00   Max.   :1.0000

#### `cdf.efro()`: Efromovich (2008) Series distribution estimator

``` r
par(mfrow = c(2, 4))
cdf.efro(shrub.data$Width)
#> Interval for Estimation: [0.080000, 2.960000]
#> Coefficient of difficulty due to biasing: 1.75791
#> Optimal Cutoff J: 2
#> Cutoff Value for High Frequency Terms U: 0.21102
#> Number of High Frequency Terms Included: 0
#> Interval where c is searched: [-1.100184, 2.442850]
#> Optimal Value for c: -0.99336
#> Value of the Estimator Integral: 1.00267
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" style="display: block; margin: auto;" />

### Quantile estimation

#### `qf.sen()`and `qf.sen()`: Empirical quantile estimators

``` r
par(mfrow = c(1, 2))
plot(qf.sen(shrub.data$Width), xlab = "", ylab = "", col = "blue", main = "Empirical quantile estimator (qf.sen)")
plot(qf.SBC(shrub.data$Width), xlab = "", ylab = "", col = "blue", main = "Empirical quantile estimator (qf.SBC)")
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-8-1.png" alt="Empirical quantile estimators for shrub width." width="80%" />
<p class="caption">
Empirical quantile estimators for shrub width.
</p>

</div>

#### `qf.SBCkernel()`: Kernel quantile estimator

``` r
par(mfrow = c(1, 4))
qf.SBCkernel(shrub.data$Width, kernel = "epanechnikov", bw = 0.01)
#> Interval for Estimation: [0.010000, 0.990000]
#> 
#> Call:
#>  qf.SBCkernel(y = shrub.data$Width, bw = 0.01, kernel = "epanechnikov")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.01
#> 
#>        x               y         
#>  Min.   :0.010   Min.   :0.1146  
#>  1st Qu.:0.255   1st Qu.:0.2326  
#>  Median :0.500   Median :0.4941  
#>  Mean   :0.500   Mean   :0.6016  
#>  3rd Qu.:0.745   3rd Qu.:0.8447  
#>  Max.   :0.990   Max.   :1.9721
qf.SBCkernel(shrub.data$Width, kernel = "epanechnikov", bw = 0.05)
#> Interval for Estimation: [0.010000, 0.990000]
#> 
#> Call:
#>  qf.SBCkernel(y = shrub.data$Width, bw = 0.05, kernel = "epanechnikov")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.05
#> 
#>        x               y          
#>  Min.   :0.010   Min.   :0.07508  
#>  1st Qu.:0.255   1st Qu.:0.23229  
#>  Median :0.500   Median :0.50193  
#>  Mean   :0.500   Mean   :0.58866  
#>  3rd Qu.:0.745   3rd Qu.:0.84423  
#>  Max.   :0.990   Max.   :1.62474
qf.SBCkernel(shrub.data$Width, kernel = "epanechnikov", bw = 0.10)
#> Interval for Estimation: [0.010000, 0.990000]
#> 
#> Call:
#>  qf.SBCkernel(y = shrub.data$Width, bw = 0.1, kernel = "epanechnikov")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.1
#> 
#>        x               y          
#>  Min.   :0.010   Min.   :0.07501  
#>  1st Qu.:0.255   1st Qu.:0.23734  
#>  Median :0.500   Median :0.49917  
#>  Mean   :0.500   Mean   :0.57216  
#>  3rd Qu.:0.745   3rd Qu.:0.84096  
#>  Max.   :0.990   Max.   :1.36825
qf.SBCkernel(shrub.data$Width, kernel = "epanechnikov", bw = 0.25)
#> Interval for Estimation: [0.010000, 0.990000]
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-9-1.png" alt="Kernel quantile estimator for shrub width." width="100%" />
<p class="caption">
Kernel quantile estimator for shrub width.
</p>

</div>

    #> 
    #> Call:
    #>  qf.SBCkernel(y = shrub.data$Width, bw = 0.25, kernel = "epanechnikov")
    #> 
    #> Data: y (89 obs.);   Bandwidth 'bw' = 0.25
    #> 
    #>        x               y          
    #>  Min.   :0.010   Min.   :0.08527  
    #>  1st Qu.:0.255   1st Qu.:0.25643  
    #>  Median :0.500   Median :0.49650  
    #>  Mean   :0.500   Mean   :0.53028  
    #>  3rd Qu.:0.745   3rd Qu.:0.80851  
    #>  Max.   :0.990   Max.   :0.99935

### Sparsity estimation

#### `sp.akbaripi()`: Akbari et al. (2019) plug-in sparsity estimator

``` r
par(mfrow = c(1, 4))
sp.akbaripi(shrub.data$Width, bw = 0.01, kernel = "epanechnikov")
#> Interval for Estimation: [0.120000, 2.540000]
#> Empirical QF 
#> Call: sp.akbaripi(shrub.data$Width, bw = 0.01, kernel = "epanechnikov")
#>  x[1:89] = 0.059043, 0.10333, 0.14062,  ..., 0.99721,      1
sp.akbaripi(shrub.data$Width, bw = 0.05, kernel = "epanechnikov")
#> Interval for Estimation: [0.120000, 2.540000]
#> Empirical QF 
#> Call: sp.akbaripi(shrub.data$Width, bw = 0.05, kernel = "epanechnikov")
#>  x[1:89] = 0.059043, 0.10333, 0.14062,  ..., 0.99721,      1
sp.akbaripi(shrub.data$Width, bw = 0.10, kernel = "epanechnikov")
#> Interval for Estimation: [0.120000, 2.540000]
#> Empirical QF 
#> Call: sp.akbaripi(shrub.data$Width, bw = 0.1, kernel = "epanechnikov")
#>  x[1:89] = 0.059043, 0.10333, 0.14062,  ..., 0.99721,      1
sp.akbaripi(shrub.data$Width, bw = 0.25, kernel = "epanechnikov")
#> Interval for Estimation: [0.120000, 2.540000]
#> Empirical QF 
#> Call: sp.akbaripi(shrub.data$Width, bw = 0.25, kernel = "epanechnikov")
#>  x[1:89] = 0.059043, 0.10333, 0.14062,  ..., 0.99721,      1
```

#### `sp.akbarikernel()`: Kernel sparsity estimator

``` r
par(mfrow = c(1, 4))
sp.akbarikernel(shrub.data$Width, bw = 0.01, kernel = "epanechnikov")
#> Interval for Estimation: [0.010000, 0.990000]
#> 
#> Call:
#>  sp.akbarikernel(y = shrub.data$Width, bw = 0.01, kernel = "epanechnikov")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.01
#> 
#>        x               y          
#>  Min.   :0.010   Min.   : 0.0000  
#>  1st Qu.:0.255   1st Qu.: 0.1128  
#>  Median :0.500   Median : 1.0578  
#>  Mean   :0.500   Mean   : 2.0284  
#>  3rd Qu.:0.745   3rd Qu.: 2.7247  
#>  Max.   :0.990   Max.   :16.3964
sp.akbarikernel(shrub.data$Width, bw = 0.05, kernel = "epanechnikov")
#> Interval for Estimation: [0.010000, 0.990000]
#> 
#> Call:
#>  sp.akbarikernel(y = shrub.data$Width, bw = 0.05, kernel = "epanechnikov")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.05
#> 
#>        x               y           
#>  Min.   :0.010   Min.   : 0.06824  
#>  1st Qu.:0.255   1st Qu.: 0.91352  
#>  Median :0.500   Median : 1.42181  
#>  Mean   :0.500   Mean   : 2.18228  
#>  3rd Qu.:0.745   3rd Qu.: 1.96606  
#>  Max.   :0.990   Max.   :13.38425
sp.akbarikernel(shrub.data$Width, bw = 0.10, kernel = "epanechnikov")
#> Interval for Estimation: [0.010000, 0.990000]
#> 
#> Call:
#>  sp.akbarikernel(y = shrub.data$Width, bw = 0.1, kernel = "epanechnikov")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.1
#> 
#>        x               y         
#>  Min.   :0.010   Min.   :0.4573  
#>  1st Qu.:0.255   1st Qu.:1.0439  
#>  Median :0.500   Median :1.2886  
#>  Mean   :0.500   Mean   :2.0995  
#>  3rd Qu.:0.745   3rd Qu.:1.8851  
#>  Max.   :0.990   Max.   :8.8587
sp.akbarikernel(shrub.data$Width, bw = 0.25, kernel = "epanechnikov")
#> Interval for Estimation: [0.010000, 0.990000]
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-11-1.png" alt="@akbari2019 kernel sparsity estimator for shrub width." width="100%" />
<p class="caption">
Akbari et al. (2019) kernel sparsity estimator for shrub width.
</p>

</div>

    #> 
    #> Call:
    #>  sp.akbarikernel(y = shrub.data$Width, bw = 0.25, kernel = "epanechnikov")
    #> 
    #> Data: y (89 obs.);   Bandwidth 'bw' = 0.25
    #> 
    #>        x               y         
    #>  Min.   :0.010   Min.   :0.5466  
    #>  1st Qu.:0.255   1st Qu.:0.9148  
    #>  Median :0.500   Median :1.2725  
    #>  Mean   :0.500   Mean   :1.9144  
    #>  3rd Qu.:0.745   3rd Qu.:2.5075  
    #>  Max.   :0.990   Max.   :4.7694

#### `sp.SBC()`: Kernel sparsity estimator based on Cox (2005) distribution estimator

``` r
par(mfrow = c(1, 4))
sp.SBC(shrub.data$Width, bw = 0.01, kernel = "epanechnikov")
#> Interval for Estimation: [0.010000, 0.990000]
#> 
#> Call:
#>  sp.SBC(y = shrub.data$Width, bw = 0.01, kernel = "epanechnikov")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.01
#> 
#>        x               y           
#>  Min.   :0.010   Min.   : 0.00000  
#>  1st Qu.:0.255   1st Qu.: 0.09031  
#>  Median :0.500   Median : 1.02242  
#>  Mean   :0.500   Mean   : 1.97514  
#>  3rd Qu.:0.745   3rd Qu.: 2.50028  
#>  Max.   :0.990   Max.   :15.73758
sp.SBC(shrub.data$Width, bw = 0.05, kernel = "epanechnikov")
#> Interval for Estimation: [0.010000, 0.990000]
#> 
#> Call:
#>  sp.SBC(y = shrub.data$Width, bw = 0.05, kernel = "epanechnikov")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.05
#> 
#>        x               y          
#>  Min.   :0.010   Min.   :0.02275  
#>  1st Qu.:0.255   1st Qu.:0.70417  
#>  Median :0.500   Median :1.24700  
#>  Mean   :0.500   Mean   :1.90523  
#>  3rd Qu.:0.745   3rd Qu.:1.92197  
#>  Max.   :0.990   Max.   :8.86825
sp.SBC(shrub.data$Width, bw = 0.10, kernel = "epanechnikov")
#> Interval for Estimation: [0.010000, 0.990000]
#> 
#> Call:
#>  sp.SBC(y = shrub.data$Width, bw = 0.1, kernel = "epanechnikov")
#> 
#> Data: y (89 obs.);   Bandwidth 'bw' = 0.1
#> 
#>        x               y         
#>  Min.   :0.010   Min.   :0.2569  
#>  1st Qu.:0.255   1st Qu.:0.6932  
#>  Median :0.500   Median :1.2355  
#>  Mean   :0.500   Mean   :1.8130  
#>  3rd Qu.:0.745   3rd Qu.:1.8639  
#>  Max.   :0.990   Max.   :6.4130
sp.SBC(shrub.data$Width, bw = 0.25, kernel = "epanechnikov")
#> Interval for Estimation: [0.010000, 0.990000]
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-12-1.png" alt="kernel sparsity estimator based on @cox2005 distribution estimator for shrub width." width="100%" />
<p class="caption">
kernel sparsity estimator based on Cox (2005) distribution estimator for
shrub width.
</p>

</div>

    #> 
    #> Call:
    #>  sp.SBC(y = shrub.data$Width, bw = 0.25, kernel = "epanechnikov")
    #> 
    #> Data: y (89 obs.);   Bandwidth 'bw' = 0.25
    #> 
    #>        x               y         
    #>  Min.   :0.010   Min.   :0.2548  
    #>  1st Qu.:0.255   1st Qu.:0.7867  
    #>  Median :0.500   Median :1.2541  
    #>  Mean   :0.500   Mean   :1.6346  
    #>  3rd Qu.:0.745   3rd Qu.:2.5482  
    #>  Max.   :0.990   Max.   :3.6553

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-akbari2019" class="csl-entry">

Akbari, M., Rezaei, M., Jomhoori, S., and Fakoor, V. (2019),
“Nonparametric estimators for quantile density function under
length-biased sampling,” *Communications in Statistics - Theory and
Methods*, Taylor & Francis, 48, 4918–4935.

</div>

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

<div id="ref-efromovich2008" class="csl-entry">

Efromovich, S. (2008), *Nonparametric curve estimation: Methods, theory,
and applications*, Springer Science & Business Media.

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

</div>
