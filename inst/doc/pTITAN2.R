## ----label=setup, include = FALSE---------------------------------------------
library(knitr)
loadNamespace("data.table")
knitr::opts_chunk$set(collapse = TRUE)

## ----label = "workflow", eval = FALSE, fig.cap = "The series of steps needed to impliment pTITAN2", echo = FALSE, results = "asis"----
#  DiagrammeR::grViz(diagram = "digraph flowchart {
#                     node [fontname = arial, shape = oval]
#                     tab1 [label = '@@1']
#                     tab2 [label = '@@2']
#                     tab3 [label = '@@3']
#                     tab4 [label = '@@4']
#                     tab5 [label = '@@5']
#                     tab6 [label = '@@6']
#                     tab7 [label = '@@7']
#                     tab8 [label = '@@8']
#  
#                     tab1 -> tab2 -> tab3 -> tab4 -> tab5 -> tab6 -> tab7 -> tab8;
#  }
#    [1]: 'Prepare and import the environmental gradient dataset into R'
#    [2]: 'Prepare and import the taxonomic dataset into R'
#    [3]: 'Preprocess raw taxonomic data to determine appropriate taxonomic level\\nof resolution (occurrence function)'
#    [4]: 'Select columns for the taxon level retuned by occurrence function'
#    [5]: 'Permut the data across treatment labels to generate list of lists'
#    [6]: 'Set up cluster for parallele processing (optional)'
#    [7]: 'Run TITAN2 series on original and permuted data sets'
#    [8]: 'Analyze probability of exceeding observed diffrence in change-point\\nbetween treatments based on distribution of paired change-point differences'
#    ")

## ----label = "datasets", echo = FALSE, results = "asis"-----------------------
d <-
  data.table::fread(text = "
C_IC_D_06_wID|C_IC_D_06_wID.csv|Environmental Gradient|Chaparral|Dry
C_IC_N_06_wID|C_IC_N_06_wID.csv|Environmental Gradient|Chaparral|Normal
CD_06_Mall_wID|CD_06_Mall_wID.csv|Taxonomic|Chaparral|Dry
CN_06_Mall_wID|CN_06_Mall_wID.csv|Taxonomic|Chaparral|Normal
",
header = FALSE,
col.names = c("R Data", "csv File", "Data Type", "Region", "Treatment")
  )

knitr::kable(d, caption = "Example data sets provided in pTITAN2.")

## ----label="list-example-data-sets"-------------------------------------------
list.files(system.file("extdata", package = "pTITAN2"))

## -----------------------------------------------------------------------------
data(C_IC_D_06_wID, C_IC_N_06_wID, CD_06_Mall_wID, CN_06_Mall_wID,
     package = "pTITAN2")

str(C_IC_D_06_wID)  # Environemntal Gradient, Dry Treatment
str(C_IC_N_06_wID)  # Environemntal Gradient, Normal Treatment
dim(CD_06_Mall_wID) # Taxonomic, Dry Treatment
dim(CN_06_Mall_wID) # Taxonomic, Normal Treatment

## ---- include = FALSE---------------------------------------------------------
set.seed(42) # set the random number generator for reproduciblity
eg_permutation <-
  data.table::data.table(site = c("A", "A", "B", "C", "C", "D", "E"),
                         trt0  = c(1, 2, 1, 1, 2, 2, 2))
eg_permutation[, site_n := .N, keyby = .(site)]
eg_permutation[site_n >  1, trt1 := sample(trt0), by = .(site)]
eg_permutation[site_n == 1, trt1 := sample(trt0)]

eg_permutation[site_n >  1, trt2 := sample(trt0), by = .(site)]
eg_permutation[site_n == 1, trt2 := sample(trt0)]

eg_permutation[site_n >  1, trt3 := sample(trt0), by = .(site)]
eg_permutation[site_n == 1, trt3 := sample(trt0)]

eg_permutation[site_n >  1, trt4 := sample(trt0), by = .(site)]
eg_permutation[site_n == 1, trt4 := sample(trt0)]

eg_permutation[, site_n := NULL]

## ----label = "exampleTreatmentPermutation", echo = FALSE, results = "asis"----
knitr::kable(eg_permutation,
             caption = "Example distribution of sites and perutated treatment labels.",
             align = "c")

## -----------------------------------------------------------------------------
library(pTITAN2)
library(magrittr)
loadNamespace("data.table")
loadNamespace("dplyr")
loadNamespace("tidyr")

## -----------------------------------------------------------------------------
dim(CD_06_Mall_wID)

# top 4 rows and first 10 columns
CD_06_Mall_wID[1:4, 1:10]

## ----label = "occurrences_example"--------------------------------------------
head(occurrences(CN_06_Mall_wID[, -1]))

## -----------------------------------------------------------------------------
# Tidyverse
CN_06_Mall_wID %>%
  dplyr::select(-StationID) %>%
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = 'taxon', values_to = 'count') %>%
  dplyr::mutate(
                Class  = substr(.data$taxon, 1L, 2L),
                Order  = substr(.data$taxon, 3L, 4L),
                Family = substr(.data$taxon, 5L, 6L),
                Genus  = substr(.data$taxon, 7L, 8L)
                ) %>%
  dplyr::group_by(.data$Class, .data$Order, .data$Family, .data$Genus) %>%
  dplyr::summarise(
                   taxon = unique(.data$taxon),
                   count = sum(.data$count > 0),
                   .groups = "keep"
                   ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(.data$Class, .data$Order, .data$Family, .data$Genus)

## -----------------------------------------------------------------------------
# data.table
taxon_count <- data.table::copy(CN_06_Mall_wID)
data.table::setDT(taxon_count)
data.table::set(taxon_count, j = "StationID", value = NULL)
for(j in as.integer(which(sapply(taxon_count, is.integer)))) {
  data.table::set(taxon_count, j = j, value = as.numeric(taxon_count[[j]]))
}

taxon_count <- data.table::melt(taxon_count, variable.factor = FALSE, measure.vars = names(taxon_count), variable.name = "taxon")
taxon_count[, value := sum(value > 0), keyby = .(taxon)]
taxon_count <- unique(taxon_count)

data.table::set(taxon_count, j = "Class",  value = substr(taxon_count[["taxon"]], 1L, 2L))
data.table::set(taxon_count, j = "Order",  value = substr(taxon_count[["taxon"]], 3L, 4L))
data.table::set(taxon_count, j = "Family", value = substr(taxon_count[["taxon"]], 5L, 6L))
data.table::set(taxon_count, j = "Genus",  value = substr(taxon_count[["taxon"]], 7L, 8L))

data.table::setkeyv(taxon_count, c("taxon", "Class", "Order", "Family", "Genus"))

taxon_count[grepl("^(Ar|Bi).+", taxon)]

## -----------------------------------------------------------------------------
eg <-
  permute(taxa = list(CD_06_Mall_wID, CN_06_Mall_wID),
          envs = list(C_IC_D_06_wID, C_IC_N_06_wID),
          sid  = "StationID")
str(eg, max.level = 2)

## -----------------------------------------------------------------------------
system.file("example-scripts/permutation_example.R", package = "pTITAN2")

## -----------------------------------------------------------------------------
permutation_example

## ---- results = "hide"--------------------------------------------------------
# Run TITAN on the observed data
# limited to 2 cores for CRAN policies, end users can use more cores
# the following uses the tidyverse dialect, load and attach the magrittr
# namespace, use dplyr and tidyr namespaces explicitly.
library(magrittr)

CD_obs <-
  TITAN2::titan(
                  env     = C_IC_D_06_wID[["ImpCover"]] # chaparral envgrad dry
                , txa     = subset(CD_06_Mall_wID, select = occurrences(CD_06_Mall_wID[, -1], n = 6L)$taxon)
                , ncpus   = 2
                , numPerm = 50
                , nBoot   = 50
  )

CN_obs <-
  TITAN2::titan(
                  env     = C_IC_N_06_wID[["ImpCover"]] # chaparral envgrad dry
                , txa     = subset(CN_06_Mall_wID, select = occurrences(CN_06_Mall_wID[, -1], n = 6L)$taxon)
                , ncpus   = 2
                , numPerm = 50
                , nBoot   = 50
  )

## -----------------------------------------------------------------------------
# create a table of the median change point from TITAN for calculating the p-values
TITAN_med <-
  rbind(
          cbind(as.data.frame(CD_obs$sumz), "run" = "Tr1_CD")
        , cbind(as.data.frame(CN_obs$sumz), "run" = "Tr2_CN")
  )
TITAN_med[["sumz"]] <- sub("(.+)(\\d)$", "\\1", rownames(TITAN_med))


TITAN_med <-
  TITAN_med %>%
  dplyr::select(run, sumz, `0.50`) %>%
  tidyr::pivot_wider(names_from = run, values_from = "0.50") %>%
  dplyr::mutate(T1T2_abs = abs(Tr1_CD - Tr2_CN)) %>%
  print

## -----------------------------------------------------------------------------
# Create a summary table of the permutation data
tr1_filt <-
  permutation_example %>%
  dplyr::select(permutation, `trt1cpsumz-`,	`trt1cpsumz+`) %>%
  tidyr::pivot_longer(cols      = `trt1cpsumz-`:`trt1cpsumz+`,
                      values_to = "Tr1_CD",
                      names_to  = "sumz")

tr1_filt[which(tr1_filt$sumz == "trt1cpsumz-"), 2] <- "fsumz-"
tr1_filt[which(tr1_filt$sumz == "trt1cpsumz+"), 2] <- "fsumz+"

tr2_filt <-
  permutation_example %>%
  dplyr::select(permutation, `trt2cpsumz-`,	`trt2cpsumz+`) %>%
  tidyr::pivot_longer(cols      = `trt2cpsumz-`:`trt2cpsumz+`,
                      values_to = "Tr2_CN",
                      names_to  = "sumz")
tr2_filt[which(tr2_filt$sumz == "trt2cpsumz-"), 2] <- "fsumz-"
tr2_filt[which(tr2_filt$sumz == "trt2cpsumz+"), 2] <- "fsumz+"


# mutate the data and create a summary table for calculating p-values
out_perm <-
  dplyr::left_join(tr1_filt, tr2_filt) %>%
  dplyr::mutate(T1T2_abs = abs(Tr1_CD - Tr2_CN)) %>%
  dplyr::select(-permutation)

# Calulate the p-values
dplyr::tibble(treatment = "T1T2_CDCN",
              "fsumz-" = (sum((
                               dplyr::filter(out_perm, sumz == "fsumz-")$T1T2_abs >
                               (dplyr::filter(TITAN_med, sumz == "fsumz-")$T1T2_abs)
                               )) + 1) / 1001, "fsumz+" =
              (sum((
                    dplyr::filter(out_perm, sumz == "fsumz+")$T1T2_abs >
                    (dplyr::filter(TITAN_med, sumz == "fsumz+")$T1T2_abs)
                    )) + 1) / 1001)

