icesVocab::getCodeDetail("SpecWoRMS", 126285)$detail$Description # Conger conger

spnorthdatras.data <- icesDatras::getDATRAS(
  record = "HL",
  survey = "SP-NORTH",
  year = 2024,
  quarter = 1:4
) %>%
  dplyr::filter(Valid_Aphia == 126285) %>%
  dplyr::mutate(dplyr::across(where(is.double), ~ na_if(.x, -9))) |>
  dplyr::mutate(dplyr::across(where(is.character), ~ na_if(.x, "-9"))) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(c(
    "RecordType", "Survey", "Quarter", "Country", "Ship", "Gear", "DoorType",
    "StNo", "Year", "SpecCodeType", "SpecCode", "SpecVal", "Sex",
    "CatIdentifier", "LngtCode", "Valid_Aphia"
  )), factor)) %>%
  dplyr::select(where(~ n_distinct(.) > 1))


duplicated(t(spnorthdatras.data))

spnorthdatras.data <- spnorthdatras.data %>%
  select(which(!duplicated(t(spnorthdatras.data))))


colnames(spnorthdatras.data)

# output for the package
usethis::use_data(spnorthdatras.data, overwrite = TRUE)
