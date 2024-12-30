suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(nestedcv)
  library(purrr)
  library(readr)
  library(tibble)
})

options(readr.show_col_types = FALSE)

# Random seed for cross-validation
set.seed(42)

# Range for defining genetic variants as cis
CIS_WINDOW <- 2e+4

annot <- readr::read_tsv("../data/albert2018/processed/albert2018_genes.tsv")
expr_irn_rc_train <- readr::read_tsv("../data/albert2018/interim/albert2018_expression_logtpm_irn_regcov_train.tsv")
gen_train <- readr::read_tsv("../data/albert2018/interim/albert2018_genotypes_train.tsv")

res <- list()
res_coef <- list()

for (gene in colnames(expr_irn_rc_train)) {
  gene_expr <- expr_irn_rc_train[[gene]]
  
  gene_annot <- dplyr::filter(annot, gene_id == gene)
  gene_chr <- dplyr::pull(gene_annot, chromosome_name)
  gene_tss <- dplyr::pull(gene_annot, transcription_start_site)

  gene_var <- gen_train |>
    dplyr::filter(chr == gene_chr, dplyr::between(pos, gene_tss - CIS_WINDOW, gene_tss + CIS_WINDOW)) |>
    dplyr::select(-c(chr, pos)) |>
    tibble::column_to_rownames("variant_id") |>
    t()

  if (ncol(gene_var) == 0) {
    res <- append(res, list(data.frame(
      gene_id = gene,
      n_var = 0,
      n_var_sig = NA,
      r2 = NA,
      rmse = NA,
      alpha = NA,
      lambda = NA
    )))
    
    next
  }

  res_glm <- nestedcv::nestcv.glmnet(
    y = gene_expr,
    x = gene_var,
    family = "gaussian",
    alphaSet = seq(0.1, 0.9, 0.2),
    n_inner_folds = 5,
    n_outer_folds = 5,
    cv.cores = 20
  )

  res <- append(res, list(data.frame(
    gene_id = gene,
    n_var = ncol(gene_var),
    n_var_sig = ifelse(is.null(nrow(res_glm[["final_coef"]])), 0, nrow(res_glm[["final_coef"]]) - 1),
    r2 = purrr::pluck(res_glm, "summary", "R.squared"),
    rmse = purrr::pluck(res_glm, "summary", "RMSE"),
    alpha = purrr::pluck(res_glm, "final_param", "alpha"),
    lambda = purrr::pluck(res_glm, "final_param", "lambda")
  )))
  
  if (!is.null(nrow(res_glm[["final_coef"]]))) {
    res_coef <- append(res_coef, list(
      res_glm[["final_coef"]] |>
        tibble::rownames_to_column("variant_id") |>
        dplyr::select(-meanExp) |>
        dplyr::mutate(gene_id = gene, .before = 1)
    ))
  }
}

res <- dplyr::bind_rows(res)
res_coef <- dplyr::bind_rows(res_coef)

readr::write_tsv(res, "../results/albert2018/expr_cisgen_glm_stats.tsv")
readr::write_tsv(res, "../results/albert2018/expr_cisgen_glm_coefs.tsv")
