#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(gtsummary))

suppressMessages(theme_gtsummary_journal(journal = "jama"))
suppressMessages(theme_gtsummary_compact())


table1 <- function(df) {
    tbl1_vars <- c(
        "stage_t",
        "stage_n",
        "metastases_YES",
        "tumor_location_Right.colon",
        "pos_margin_YES",
        "lvi_YES",
        "msi_pos",
        "surv_months_5_yr",
        "surv_status_5_yr"
    )
    df %>%
        select(all_of(tbl1_vars)) %>%
        mutate(
            tumor_location_Right.colon = if_else(
                tumor_location_Right.colon == 0,
                "Left colon",
                "Right colon"
            ),
            .keep = "unused"
        ) %>%
        tbl_summary(
            by = "tumor_location_Right.colon",
            label = tbl_labels(tbl1_vars)
        ) %>%
        modify_header(update = list(label ~ "**Variable**"))
}



# tbl_regression() only wants "label" to include variables that are
# present in the data. This will take the vars in the model and return
# the appropriate human readable labels
tbl_labels <- function(final_vars) {
    label_map <- list(
        "tumor_location_Right.colon" = tumor_location_Right.colon ~ "Right colon",
        "stage_t" = stage_t ~ "T stage",
        "stage_n" = stage_n ~ "N stage",
        "metastases_YES" = metastases_YES ~ "Metastases",
        "pos_margin_YES" = pos_margin_YES ~ "R1/R2 Margin",
        "lvi_YES" = lvi_YES ~ "Lymphovascular invasion",
        "msi_pos" = msi_pos ~ "Microsatellite instability",
        "apc_Mutated" = apc_Mutated ~ "APC mutation",
        "csmd1_Mutated" = csmd1_Mutated ~ "CSMD1 mutation",
        "csmd3_Mutated" = csmd3_Mutated ~ "CSMD3 mutation",
        "dnah11_Mutated" = dnah11_Mutated ~ "DNAH11 mutation",
        "dnah5_Mutated" = dnah5_Mutated ~ "DNAH5 mutation",
        "fat3_Mutated" = fat3_Mutated ~ "FAT3 mutation",
        "fat4_Mutated" = fat4_Mutated ~ "FAT4 mutation",
        "fbxw7_Mutated" = fbxw7_Mutated ~ "FBXW7 mutation",
        "flg_Mutated" = flg_Mutated ~ "FLG mutation",
        "kras_Mutated" = kras_Mutated ~ "KRAS mutation",
        "lrp_Mutated" = lrp_Mutated ~ "LRP mutation",
        "muc16_Mutated" = muc16_Mutated ~ "MUC16 mutation",
        "neb_Mutated" = neb_Mutated ~ "NEB mutation",
        "obscn_Mutated" = obscn_Mutated ~ "OBSCN mutation",
        "pclo_Mutated" = pclo_Mutated ~ "PCLO mutation",
        "pi3k_Mutated" = pi3k_Mutated ~ "PI3K mutation",
        "ryr_Mutated" = ryr_Mutated ~ "RYR mutation",
        "smad4_Mutated" = smad4_Mutated ~ "SMAD4 mutation",
        "spta1_Mutated" = spta1_Mutated ~ "SPTA1 mutation",
        "syne1_Mutated" = syne1_Mutated ~ "SYNE1 mutation",
        "tp53_Mutated" = tp53_Mutated ~ "TP53 mutation",
        "ttn_Mutated" = ttn_Mutated ~ "TTN mutation",
        "unc13c_Mutated" = unc13c_Mutated ~ "UNC13C mutation",
        "ush2a_Mutated" = ush2a_Mutated ~ "USH2A mutation",
        "zfh_Mutated" = zfh_Mutated ~ "ZFH mutation",
        "surv_status_5_yr" = surv_status_5_yr ~ "Deaths (within 5 years)",
        "surv_months_5_yr" = surv_months_5_yr ~ "Follow-up time"
    )
    enframe(label_map, name = "term", value = "label") %>%
        filter(term %in% final_vars) %>%
        pull(label)
}


# requires df with *only* the Y vars (surv months and status) and X vars
cox_mod <- function(df) {
    coxph(Surv(surv_months_5_yr, surv_status_5_yr) ~ ., data = df)
}


schoen <- function(model, ...) {
    ggcoxzph(cox.zph(model), ...)
}


save_gt_table <- function(model_tbl, filename) {
    model_tbl %>%
        as_gt() %>%
        gt::gtsave(filename, path = "../output")
}


save_surv_plot <- function(plt, filename, ...) {
    ggsave(filename, path = "../output", plot = print(plt), ...)
}


cox_tbl <- function(model, ...) {
    tbl_regression(model, exponentiate = TRUE, ...) %>%
        modify_header(update = list(label ~ "**Variable**"))
}


# song and dance to have the terms print in the table in the order we
# want
order_term_for_table <- function(df, rl_first = TRUE) {
    # rl vs not from chisq.R analyses
    rl_genes <- c(
        "csmd3_Mutated",
        "fat3_Mutated",
        "fat4_Mutated",
        "muc16_Mutated",
        "neb_Mutated",
        "obscn_Mutated",
        "pclo_Mutated",
        "pi3k_Mutated",
        "syne1_Mutated",
        "tp53_Mutated",
        "ush2a_Mutated",
        "zfh_Mutated"
    )
    non_rl_genes <- c(
        "apc_Mutated",
        "csmd1_Mutated",
        "dnah11_Mutated",
        "dnah5_Mutated",
        "fbxw7_Mutated",
        "flg_Mutated",
        "kras_Mutated",
        "lrp_Mutated",
        "ryr_Mutated",
        "smad4_Mutated",
        "spta1_Mutated",
        "ttn_Mutated",
        "unc13c_Mutated"
    )
    other_vars <- c(
        "tumor_location_Right.colon",
        "stage_t",
        "stage_n",
        "metastases_YES",
        "pos_margin_YES",
        "lvi_YES",
        "msi_pos"
    )
    if (rl_first) {
        term_levels <- c(other_vars, rl_genes, non_rl_genes)
    } else {
        term_levels <- c(other_vars, sort(c(rl_genes, non_rl_genes)))
    }
    df %>%
        mutate(term = parse_factor(term, levels = term_levels)) %>%
        arrange(term) %>%
        mutate(term = as.character(term))
}


one_step_select_genes <- function(full_model) {
    full_model %>%
        tidy() %>%
        # genes are what we're using backward selection for so only drop
        # those with p < 0.5
        filter(
                str_detect(term, "Mutated", negate = TRUE) |
                (str_detect(term, "Mutated") & p.value < 0.5)
        ) %>%
        order_term_for_table(rl_first = FALSE) %>%
        pull(term)
}


km_plot <- function(model, dat) {
    ggsurvplot(
        model,
        data = dat,
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


km_table <- function(model) {
    tbl_survfit(
        model,
        times = seq.int(12, 60, 12),
        label = "Survival",
        label_header = "{time} months"
    ) %>%
    modify_header(update = list(label ~ ""))
}


km_data <- function(df) {
    mutate(df, 
        Side = if_else(tumor_location_Right.colon == 0, "Left", "Right"),
        .keep = "unused"
    ) 
}


km_fit <- function(df) {
    survfit(
        Surv(surv_months_5_yr, surv_status_5_yr) ~ Side,
        data = df
    )
}


# causal effect of tumor side and survival
rl_to_surv <- function(df) {
    # cox analyses
    dat <- inputs %>%
        select(ends_with("_5_yr"), tumor_location_Right.colon)
    model <- cox_mod(dat)
    schoen_plot <- schoen(model, caption = "RL on survival") 
    model_tbl <- cox_tbl(model, 
        label = list(tumor_location_Right.colon ~ "Right Cancer")
    )
    save_gt_table(model_tbl, "rl_to_surv.html")
    save_surv_plot(schoen_plot, "rl_to_surv_schoen.png")
    # km analyses
    km_dat <- km_data(dat)
    km_mod <- km_fit(km_dat)
    km_table(km_mod) %>% save_gt_table("survival_probabilities_rl.html")
    km_plot(km_mod, km_dat) %>%
        save_surv_plot("km_plot.png", height = 7.5, width = 6)
}


total_effect_g_on_surv <- function(df) {
    dat <- df %>%
        select(
            ends_with("_5_yr"),
            tumor_location_Right.colon,
            ends_with("Mutated"),
            msi_pos
        )
    full_model <- cox_mod(dat)
    final_vars <- one_step_select_genes(full_model)
    final_model <- select(dat, ends_with("_5_yr"), all_of(final_vars)) %>% cox_mod()
    schoen_plot <- schoen(final_model, caption = "Genes total effect on survival") 
    # only print genes in table because we are looking at the total effect of
    # the gene to avoid Table 2 fallacy
    model_tbl <- cox_tbl(
        final_model,
        label = tbl_labels(final_vars),
        include = contains("mutated")
    )
    save_gt_table(model_tbl, "total_effect_g_on_surv.html")
    save_surv_plot(schoen_plot, "total_effect_g_on_surv.png")
}


inputs <- read_rds("../output/model_inputs.rds")
# Summarise patients
inputs %>% table1() %>% save_gt_table("table_1.html")
# Causal: direct effect RL -> all paths -> survival
rl_to_surv(inputs)
# Causal: total effect G -> survival
total_effect_g_on_surv(inputs)
