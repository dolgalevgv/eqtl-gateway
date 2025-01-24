suppressPackageStartupMessages({
  library(FastGGM)
  library(readr)
  library(RcppParallel)
})

options(readr.show_col_types = FALSE)

RcppParallel::setThreadOptions(numThreads = 24)


expr_irn_rc_train <- readr::read_tsv("../../data/albert2018/interim/albert2018_expression_logtpm_irn_regcov_train.tsv")

ggm <- FastGGM::FastGGM_Parallel(as.matrix(expr_irn_rc_train), lambda = 0.9)

ggm[["gene_id"]] <- colnames(expr_irn_rc_train)

saveRDS(ggm, "../../results/albert2018/expr_ggm_lambda08.rds")
