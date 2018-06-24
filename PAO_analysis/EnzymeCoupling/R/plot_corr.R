library(magrittr)
load("data/pao_g0_sets.RData")
load("data/pao_gr0_sets.RData")
load("data/pao_g1_sets.RData")
load("data/pao_gr1_sets.RData")

clean_sets <- function(sets, pattern="PA\\d{4}") {
  genes_only <- function(x) {
    x %>%
      stringr::str_extract(pattern) %>%
      purrr::discard(is.na)
  }
  sets %>%
    purrr::map(genes_only) %>%
    purrr::discard(~length(.) == 0)
}

no_diag <- function(X) {
  diag(X) = NA
  x <- as.numeric(X)
  x[!is.na(x)]
}

plot_in_out <- function(C, sets) {
  plot_title <- dplyr::quo_name(dplyr::enquo(sets))
  all_genes <- dimnames(C)[[1]]
  set_genes <- intersect(unique(purrr::flatten_chr(sets)), all_genes)
  sets <- purrr::map(sets, intersect, set_genes)
  noset_genes <- setdiff(all_genes, set_genes)
  setC = C[set_genes, set_genes]

  nonsing_sets <- sets %>% purrr::discard(~length(.) < 2)
  in_corrs = nonsing_sets %>%
    purrr::map(~ no_diag(setC[.,.])) %>%
    purrr::flatten_dbl()
  out_corrs = nonsing_sets %>%
    purrr::map(~ as.numeric(setC[.,setdiff(set_genes, .)])) %>%
    purrr::flatten_dbl()

  plot(density(no_diag(C)), main=plot_title)
  lines(density(no_diag(setC)), col="red")
  lines(density(in_corrs), col="blue")
  lines(density(out_corrs), col="green")
}

g0 <- clean_sets(g0_sets)
gr0 <- clean_sets(gr0_sets)
g1 <- clean_sets(pao_g1_sets)
gr1 <- clean_sets(gr1_sets)

plot_in_out(C, g0)
plot_in_out(C, gr0)
plot_in_out(C, g1)
plot_in_out(C, gr1)

breaks <- 2:30
set_max <- function(x, m) {x[x > m] <- m; x}
multi_len <- function(x) {
  x %>%
    purrr::map_int(length) %>%
    purrr::discard(~. < 2) %>%
    set_max(max(breaks))
}
hist(multi_len(g0), breaks=breaks)
hist(multi_len(gr0), breaks=breaks)
hist(multi_len(g1), breaks=breaks)
hist(multi_len(gr1), breaks=breaks)


