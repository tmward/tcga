#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gtsummary))

suppressMessages(theme_gtsummary_journal(journal = "jama"))
suppressMessages(theme_gtsummary_compact())


table1 <- function(df) {
    tbl1_vars <- c(
        "diagnosis_age",
        "sex",
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
            tumor_location = parse_factor(
                if_else(tumor_location_Right.colon == 0,
                    "Left",
                    "Right"
                ),
                levels = c("Right", "Left")
            ),
            .keep = "unused"
        ) %>%
        tbl_summary(
            by = "tumor_location",
            label = list(
                diagnosis_age ~ "Age",
                sex ~ "Sex",
                stage_t ~ "T stage",
                stage_n ~ "N stage",
                metastases_YES ~ "Metastases",
                pos_margin_YES ~ "R1/R2 Margin",
                lvi_YES ~ "Lymphovascular invasion",
                msi_pos ~ "Microsatellite instability",
                surv_status_5_yr ~ "Deaths (within 5 years)",
                surv_months_5_yr ~ "Follow-up time"
            )
        ) %>%
        modify_header(update = list(label ~ "**Variable**")) %>%
        modify_spanning_header(c("stat_1", "stat_2") ~ "**Tumor Site**")
}

inputs <- read_rds("../output/model_inputs.rds")
age_sex <- read_rds("../output/tidied.rds") %>%
    distinct(caseid, diagnosis_age, sex) %>%
    mutate(sex = fct_drop(sex))
# Summarise patients
inputs %>%
    inner_join(age_sex, by = "caseid") %>%
    table1() %>%
    as_gt() %>%
    gt::gtsave("table1.html", path = "../output")
