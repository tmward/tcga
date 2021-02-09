#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(gtsummary))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(magrittr))

suppressMessages(theme_gtsummary_journal(journal = "jama"))
suppressMessages(theme_gtsummary_compact())


save_surv_plot <- function(plt, filename, ...) {
    ggsave(filename = filename, path = "../output", plot = print(plt), ...)
}


save_gt_table <- function(df, filename) {
    df %>%
        as_gt() %>%
        gt::gtsave(filename, path = "../output")
}


table_1 <- function(df) {
    df %>%
        select(
            stage_t,
            stage_n,
            metastases_YES,
            tumor_location_Right.colon,
            pos_margin_YES,
            lvi_YES,
            msi_pos,
            surv_months_5_yr,
            surv_status_5_yr
        ) %>%
        mutate(
            tumor_location = if_else(
                tumor_location_Right.colon == 0,
                "Left colon",
                "Right colon"
            ),
            .keep = "unused"
        ) %>%
        tbl_summary(
            by = "tumor_location",
            label = list(
                stage_t ~ "T Stage",
                stage_n ~ "N Stage",
                metastases_YES ~ "Metastatic disease",
                lvi_YES ~ "Lymphovascular invasion",
                pos_margin_YES ~ "R1/R2 margin",
                msi_pos ~ "Microsatellite unstable",
                surv_status_5_yr ~ "Deaths (within 5 years)",
                surv_months_5_yr ~ "Follow-up time"
            )
        ) %>%
        modify_header(update = list(label ~ "**Variable**"))
}


univariate_data <- function(df) {
    df %>%
        select(
            tumor_location_Right.colon,
            surv_months_5_yr,
            surv_status_5_yr
        ) %>%
        mutate(
            Side = if_else(tumor_location_Right.colon == 0, "Left", "Right"),
            .keep = "unused"
        ) 
}


surv_mod <- function(df) {
    list(
        data = df,
        fit = survfit(
            Surv(surv_months_5_yr, surv_status_5_yr) ~ Side,
            data = df
        )
    )
}


save_km_table <- function(model) {
    tbl_survfit(
        model$fit,
        times = seq.int(12, 60, 12),
        label = "Survival",
        label_header = "{time} months"
    ) %>%
    modify_header(update = list(label ~ "")) %>%
    save_gt_table("survival_probabilities_rl.html")
}


km_plot <- function(model) {
    ggsurvplot(
        model$fit,
        data = model$data,
        palette = c("#B57EDC", "#FF0000"),
        risk.table = TRUE,
        legend.labs = c("Left", "Right"),
        xlab = "Time (months)",
        break.time.by = 12,
        ylim = c(0.50, 1.0),
        pval = TRUE,
        pval.coord = c(0, 0.65),
        pval.method = TRUE,
        pval.method.coord = c(0, 0.7),
        cumevents = TRUE,
        tables.height = 0.2
    )
}


univariate_coxph <- function(df) {
    surv_mod <- coxph(
        Surv(surv_months_5_yr, surv_status_5_yr) ~ Side,
        data = df
    )
    # save Schoenfeld residual plots testing for Coxph assumptions
    ggcoxzph(cox.zph(surv_mod), caption = "Univariate Cox") %>%
        save_surv_plot("coxph_univariate_schoen.png")

    # make table showing regression results
    surv_mod %>%
        tbl_regression(
            exponentiate = TRUE,
            show_single_row = c(Side),
            label = list(Side ~ "Right Cancer")
        ) %>%
        modify_header(update = list(label ~ "**Variable**"))
}


multivariate_data <- function(df) {
    df %>%
        select(
            tumor_location_Right.colon,
            stage_t,
            stage_n,
            metastases_YES,
            pos_margin_YES,
            lvi_YES,
            msi_pos,
            surv_months_5_yr,
            surv_status_5_yr
        )
}


multivariate_coxph <- function(df) {
    surv_mod <- coxph(
        Surv(surv_months_5_yr, surv_status_5_yr) ~ .,
        data = df
    )
    # save Schoenfeld residual plots testing for Coxph assumptions
    ggcoxzph(cox.zph(surv_mod), caption = "Multivariate Cox") %>%
        save_surv_plot(
            "coxph_multivariate_schoen.png",
            width = 25,
            height = 12,
            units = "in"
        )
    # make table showing regression results
    surv_mod %>%
        tbl_regression(
            exponentiate = TRUE,
            label = list(
                stage_t ~ "T stage",
                stage_n ~ "N stage",
                metastases_YES ~ "Metastases",
                tumor_location_Right.colon ~ "Right-sided",
                lvi_YES ~ "Lymphovascular invasion",
                pos_margin_YES ~ "R1/R2 margin",
                msi_pos ~ "Microsatellite instability"
            )
        ) %>%
        modify_header(update = list(label ~ "**Variable**"))
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
        "muc_Mutated",
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
                muc_Mutated ~ "MUC16 mutation",
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


inputs <- read_rds("../output/model_inputs.rds")
# KM plot, and survival probs, R vs L
inputs %>% univariate_data() %>% surv_mod() %T>% save_km_table() %>% km_plot() %>%
    save_surv_plot("km_plot.png", height = 7.5, width = 6)
# Univariate coxph including testing model assumptions
inputs %>% univariate_data() %>% univariate_coxph() %>% save_gt_table("coxph_univariate_table.html")
# multivariate with TNM, R, site, and testing model assumptions
inputs %>% multivariate_data() %>% multivariate_coxph() %>% save_gt_table("coxph_multivariate_table.html")
# lasso regression with TNM, R, mutations
inputs %>% lasso_data() %>% lasso() %>% sli() %>% sli_table() %>% save_gt_table("coxph_lasso.html")
# summarise patients
inputs %>% table_1() %>% save_gt_table("table_1.html")
