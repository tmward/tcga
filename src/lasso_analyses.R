#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(gtsummary))
suppressPackageStartupMessages(library(survival))

suppressMessages(theme_gtsummary_journal(journal = "jama"))
suppressMessages(theme_gtsummary_compact())


save_gt_table <- function(df, filename) {
    df %>%
        as_gt() %>%
        gt::gtsave(filename, path = "../output")
}


center <- function(x) {
    x - mean(x, na.rm = TRUE)
}


lasso_data <- function(df) {
    x <- df %>%
        select(
            stage_t,
            stage_n,
            lvi_YES,
            pos_margin_YES,
            metastases_YES,
            msi_pos,
            ends_with("Mutated")
        ) %>%
        mutate(across(c(stage_t, stage_n), ~ center(.))) %>%
        as.matrix()
    y <- Surv(df$surv_months_5_yr, df$surv_status_5_yr)
    list(df = df, x = x, y = y)
}


lasso <- function(dat) {
    set.seed(1234)
    lasso_mod <- glmnet(dat$x, dat$y, family = "cox", alpha = 1)
    lasso_mod_cv <- cv.glmnet(dat$x, dat$y, family = "cox", alpha = 1) 
    list(mod = lasso_mod, lambda_min = lasso_mod_cv$lambda.min, dat = dat)
}


standard_lambda <- function(lresults) {
    lresults$lambda_min * nrow(lresults$dat$y)
}


lasso_beta <- function(lresults) {
    as.numeric(coef(lresults$mod, s = lresults$lambda_min, exact = TRUE))
}


sli <- function(lresults) {
    set.seed(1234)
    fli <- selectiveInference::fixedLassoInf(
        lresults$dat$x,
        lresults$dat$df$surv_months_5_yr,
        status = lresults$dat$df$surv_status_5_yr,
        beta = lasso_beta(lresults),
        lambda = standard_lambda(lresults),
        family = "cox",
        bits = 200
    )
    # tbl_regression expects a regression model "output" that has
    # variable "model" to see the original df so add that
    list(model = lresults$dat$df, fli = fli)
}


# custom tidier make sli results work with tbl_regression(); conf.level ignored
tidy_sli <- function(x, exponentiate = TRUE, conf.level = 0.95, ...) {
    # define levels so the terms sort in the table how we want
    sli_var_levels <- c(
        "stage_t",
        "stage_n",
        "metastases_YES",
        "pos_margin_YES",
        "dnah5_Mutated",
        "flg_Mutated",
        "muc16_Mutated",
        "neb_Mutated",
        "obscn_Mutated",
        "smad4_Mutated",
        "syne1_Mutated",
        "ttn_Mutated",
        "ush2a_Mutated"
    )
    tidied <- tibble(
            enframe(x$fli$coef0, name = "term", value = "estimate"),
            conf.low = x$fli$ci[,1],
            conf.high = x$fli$ci[,2],
            p.value = x$fli$pv
        ) %>%
        mutate(
            term = parse_factor(
                str_remove(term, "x\\[, m\\]"),
                levels = sli_var_levels
            )
        ) %>%
        arrange(term)
    if (exponentiate) {
        tidied <- tidied %>%
            mutate(across(c(estimate, conf.low, conf.high), ~ exp(.)))
    }
    tidied
}


sli_table <- function(sli_results) {
    sli_results %>%
        tbl_regression(
            label = list(
                stage_t ~ "T stage",
                stage_n ~ "N stage",
                pos_margin_YES ~ "R1/R2 margin",
                metastases_YES ~ "Metastases",
                ttn_Mutated ~ "TTN mutation",
                syne1_Mutated ~ "SYNE1 mutation",
                obscn_Mutated ~ "OBSCN mutation",
                muc16_Mutated ~ "MUC16 mutation",
                ush2a_Mutated ~ "USH2A mutation",
                dnah5_Mutated ~ "DNAH5 mutation",
                neb_Mutated ~ "NEB mutation",
                flg_Mutated ~ "FLG mutation",
                smad4_Mutated ~ "SMAD4 mutation"
            ),
            exponentiate = TRUE,
            tidy_fun = tidy_sli
        ) %>%
        modify_header(
            update = list(
                label ~ "**Variable**",
                estimate ~ "**HR**"
            )
        ) %>%
        modify_footnote(estimate ~ "HR = Hazard Ratio", abbreviation = TRUE)
}


# lasso regression with TNM, R, mutations
read_rds("../output/model_inputs.rds") %>%
    lasso_data() %>%
    lasso() %>%
    sli() %>%
    sli_table() %>%
    save_gt_table("coxph_lasso.html")
