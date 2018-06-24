
library(magrittr)

pao_trans_exps <- readr::read_csv("data/PATRIC_PAO_transcriptomics.csv")
myGetGEO <- function(GEO) {
  tryCatch({
    GEOquery::getGEO(GEO)
  }, error = function(err) {
    cat(GEO)
    return(NA)
  })
}
geo_sets <- lapply(pao_trans_exps$Accession, function(x) GEOquery::getGEO(x, destdir="cache"))

geo_sets <- geo_sets[pao_trans_exps$Strain == "[PAO1]"]
esets <- c(geo_sets, recursive=TRUE)

gatherExpression <- function(eset, featureTransform=function(x) x) {
  E <- as.data.frame(Biobase::exprs(eset))
  E$feature <- row.names(E)
  E %>%
    tidyr::gather(experiment, expression, -feature) %>%
    dplyr::mutate(feature=featureTransform(feature)) %>%
    dplyr::group_by(feature, experiment) %>%
    dplyr::summarize(expression=mean(expression))
}

edfs <- vector("list", length=length(esets))
for (i in seq_along(edfs)) {
  edfs[[i]] <- gatherExpression(esets[[i]])
}

# some GSMs are duplicated; remove them
foundGSMs <- c()
for (i in seq_along(edfs)) {
  edfs[[i]] %<>%
    dplyr::filter(!(experiment %in% foundGSMs))
  foundGSMs <- c(foundGSMs, unique(edfs[[i]]$experiment))
}

edf <- dplyr::bind_rows(edfs)
emat <- edf %>%
  dplyr::select(feature, experiment, expression) %>%
  tidyr::spread(experiment, expression) %>%
  dplyr::ungroup() %>%
  dplyr::filter(stringr::str_detect(feature, "^PA\\d{4,4}")) %>%
  dplyr::mutate(feature = substr(feature, 1, 6))

rnames <- emat$feature
emat$feature <- NULL
E <- as.matrix(emat)
row.names(E) <- rnames
for (i in 1:dim(E)[2]) {
  E[is.na(E[ ,i]),i] <- min(E[ ,i], na.rm=TRUE)
  nzmin <- min(E[E[ ,i] > 0,i])
  E[E[ ,i] < nzmin,i] <- nzmin
}
Enorm <- log(E)
for (i in 1:dim(E)[2]) {
  Enorm[ ,i] <- (Enorm[ ,i] - mean(Enorm[ ,i], na.rm=TRUE)) / sd(Enorm[ ,i], na.rm=TRUE)
}

C <- cor(t(Enorm))
