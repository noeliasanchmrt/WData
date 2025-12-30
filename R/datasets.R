#' Shrub data
#'
#' @name shrub.data
#' @docType data
#' @description Dataset containing the size of the \emph{Cercocarpus montanus} species in an ancient quarry.
#' @details
#' During the fall semester of 1986, students in a graduate course on biological sampling techniques, taught by Lyman L. McDonald at the University of Wyoming, conducted a field study using the linear transect method to measure the size of Cercocarpus montanus shrubs in an old limestone quarry located just east of Laramie, Wyoming. In this area, rock fissures run predominantly north to south, and vegetation is denser within them. To align with this structure, a baseline was established approximately parallel to the fissures, and six transects were placed perpendicular to this baseline, grouped into two independent replicates (I and II), each comprising three equally spaced transects.
#'
#' Students walked along the transects and recorded all Cercocarpus montanus individuals intersected. For each shrub, they measured maximum height, width (defined as the greatest distance between two tangents to the shrub’s contour, parallel to the transect), and the number of stems. Given the rhizomatous nature of the species and possible interconnections between neighboring shrubs, an individual was defined as a group of stems at the base separated by at least fifteen centimeters from the nearest neighbor. Additionally, the length of intersection with the transect line was recorded for each shrub cluster.
#'
#' Due to the nature of the sampling method, wider shrubs had a higher probability of being intersected by a transect. As a result, the sample of shrub widths exhibits length bias. Measurements of height and number of stems are also affected by sampling bias, although the associated bias function is more complex, as it depends on the relationship between width and these other morphological features.
#'
#' For further details on the sampling protocol and data structure, see \insertCite{muttlak1988;textual}{WData}.

#' @format A `data.frame` with the following variables:
#' \tabular{ll}{
#' `Replica` \tab {Replica identifier (I or II).} \cr \tab \cr
#' `Transect` \tab {Transect identifier (1,2 or 3).} \cr \tab \cr
#' `Number` \tab {Shrub/clump identifier.} \cr \tab \cr
#' `Intercept` \tab {Length of the intersection of the clump of overlapping shrubs with the transect.} \cr \tab \cr
#' `Width` \tab {Width between two lines tangent to the shrub and parallel to the transect. } \cr \tab \cr
#' `Height` \tab {Maximum height of the shrub encountered by the transect.} \cr \tab \cr
#' `Stems` \tab {Number of stems on the shrub encountered by the transect.} \cr
#' }
#' @usage data(shrub.data)
#' @source \url{https://www.proquest.com/docview/303721613?%20T&fromopenview=true&pq-origsite=gscholar&sourcetype=Dissertations%20}
#' @references \insertAllCited{}
#' @keywords datasets
#' @examples
#' data(shrub.data)
#' names(shrub.data)
#' str(shrub.data)
#' class(shrub.data)
"shrub.data"

#' \emph{Conger conger} 2024 ICES \emph{Spanish North Coast Bottom Trawl Survey} (DATRAS)
#'
#' @name spnorthdatras.data
#' @docType data
#' @description \emph{Conger conger} 2024 dataset from the ICES \emph{Spanish North Coast Bottom Trawl Survey}, extracted from DATRAS (Database of Trawl Surveys).
#'
#' @details DATRAS is maintained by ICES (International Council for the Exploration of the Sea) to collate, document, and standardise trawl survey data, ensure data quality, and facilitate access for stock assessment and fish community studies. It primarily contains data from bottom trawl surveys coordinated by ICES expert groups, covering the Baltic Sea, Skagerrak, Kattegat, North Sea, English Channel, Celtic Sea, Irish Sea, Bay of Biscay, and the eastern Atlantic from the Shetlands to Gibraltar. A complete list of ICES DATRAS surveys is available \href{https://datras.ices.dk/home/surveycode.aspx}{here}.
#'
#' `spnorthdatras.data` includes _Conger conger_ 2024 length data from the ICES _Spanish North Coast Bottom Trawl Survey_ (DATRAS). Data can be manually downloaded from the \insertCite{ices_datras;textual}{WData} \href{https://datras.ices.dk/Data_products/Download/Download_Data_public.aspx}{web portal} or  obtained using the `icesDatras` R package as follows:
#' ```r
#' icesDatras::getDATRAS(
#'  record = "HL",
#'  survey = "SP-NORTH",
#'  year = 2024,
#'  quarter = 1:4
#' ) %>%
#' dplyr::filter(Valid_Aphia == 126285)
#' ```
#'
#' Fishing gears are selective by species and, within species, by size. The probability that a fish is retained and hauled to the surface depends on its body size because of mesh selectivity. Although retention is more directly determined by girth, length is typically recorded in practice, and the two measures are strongly linearly related. As a result, the observed total length distributions are biased toward larger individuals because of size-dependent trawl selectivity (\insertCite{wileman1996;textual}{WData}).
#'
#' Missing values in the original DATRAS files (\code{-9}) were replaced by  \code{NA}. Some variables were converted to factors, and several fields from the original database were omitted; otherwise, the column names and formats remained aligned with the original specification.
#'
#' @format A `data.frame` with the following variables:
#' \tabular{ll}{
#' `StNo` \tab National station code or number.\cr
#' `HaulNo` \tab Haul number.\cr
#' `TotalNo` \tab Total number of individuals in the catch.\cr
#' `SubFactor` \tab Subsampling factor (\code{SubWgt / CatCatchWgt}).\cr
#' `SubWgt` \tab Subsample weight (grams).\cr
#' `CatCatchWgt` \tab Catch weight per category (grams).\cr
#' `LngtClass` \tab Length (centimetres).\cr
#' `HLNoAtLngt` \tab Number of individuals at length `LngtClass`.\cr
#' }
#' @usage data(spnorthdatras.data)
#' @source
#' \href{https://datras.ices.dk/Data_products/Download/Download_Data_public.aspx}{ICES DATRAS Download web portal}
#'
#' \href{https://ices-library.figshare.com/articles/report/SISP_15_-_Manual_of_the_IBTS_North_Eastern_Atlantic_Surveys/19051037}{Manual of the IBTS North Eastern Atlantic Surveys}
#'
#' \href{https://vocab.ices.dk}{ICES DATRAS Field Descriptions}
#'
#' \href{https://www.marinespecies.org}{World Register of Marine Species (WoRMS)}
#'
#' \href{https://cran.r-project.org/web/packages/icesVocab/index.html}{icesVocab R package: ICES Vocabularies Database Web Services}

#' @references \insertAllCited{}
#' @keywords datasets
#' @examples
#' data(spnorthdatras.data)
#' names(spnorthdatras.data)
#' str(spnorthdatras.data)
#' class(spnorthdatras.data)
"spnorthdatras.data"
