#!/usr/bin/env Rscript
"
Produces two R objects from the previously tidied patient and msi
data produced by tidy_tcga.R. Notably, it only keeps data for patients
with a right or left colon cancer and completed information for use in
survival analyses (so more than 0 months of follow-up) and completed T
and N staging. The first, analyses_data_tidy.rds, holds the subset of
data for R and L cancers with imputed missing LVI, R0 margin status,
and metastases status. This is tidy (one gene mutation per line). It
also produces a second R object, model_inputs.rds, that holds the data
widened, for use in modeling.

Usage: make_inputs.R [-h] [-t TIDIED_DF] [-m MSI]

Options:
    -h             Print this menu and exit.
    -t TIDIED_DF   Filename that holds pt info/genetics [default: ../output/tidied.rds].
    -m MSI         Filename that holds msi info [default: ../output/msi.rds].
" -> doc

library(docopt)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidymodels))


rl_only <- function(df){
    df %>%
        filter(tumor_location %in% c("Right colon", "Left colon")) %>%
        # reset factor levels to only have right and left
        mutate(tumor_location = fct_drop(tumor_location))
}


rm_t_n_5_yr_nas <- function(df) {
    df %>%
        # whittle down to same population we used for regressions
        # remove na surv_status can't use that and t and n stage (only 2
        # cases removed with the latter)
        filter(across(c(surv_status_5_yr, stage_t, stage_n), ~ !is.na(.)))
}


keep_analyses_cols <- function(df) {
    df %>%
        select(
            caseid,
            gene,
            classification,
            tumor_location,
            stage_t,
            stage_n,
            lvi,
            pos_margin,
            metastases,
            surv_months_5_yr,
            surv_status_5_yr
        )
}


missing_to_na <- function(df) {
    df %>%
        # convert "(Missing)" to NA (for factor columns, can only do this
        # once a character
        mutate(across(where(is.factor), as.character)) %>%
        mutate(across(where(is.character), ~ na_if(., "(Missing)")))
}


group_mut_subtypes <- function(df) {
    df %>%
        mutate(classification = if_else(
                classification == "Wild type", "Wild type", "Mutated"
            )
        )
}


analyses_data <- function(df) {
    df %>%
        rl_only() %>%
        rm_t_n_5_yr_nas() %>%
        keep_analyses_cols() %>%
        missing_to_na() %>%
        group_mut_subtypes()
}


model_inputs <- function(df, msi) {
    wide_df  <- df %>%
        pivot_wider(names_from = gene, values_from = classification)
    # list of genes so we know which ones to relevel to have a baseline
    # of "Wild type"
    genes <- df %>% distinct(gene) %>% pull(gene)
    set.seed(1234)
    recipe(surv_months_5_yr + surv_status_5_yr ~ ., data = wide_df) %>%
        update_role(caseid, new_role = "id") %>%
        step_relevel(pos_margin, metastases, ref_level = "NO") %>%
        step_relevel(all_of(genes), ref_level = "Wild type") %>%
        step_bagimpute(lvi, pos_margin, metastases) %>%
        step_dummy(all_nominal()) %>%
        prep(retain = TRUE) %>%
        juice() %>%
        inner_join(msi, by = "caseid")
}


add_imputed_and_msi <- function(a_df, m_df) {
    m_df <- select(m_df, caseid, lvi_YES, metastases_YES, pos_margin_YES, msi_pos)
    left_join(a_df, m_df, by = "caseid") %>%
        mutate(
            lvi = lvi_YES,
            metastases = metastases_YES,
            pos_margin = pos_margin_YES,
            .keep = "unused"
        )
}


main <- function(opts) {
    a_df <- read_rds(opts$t) %>% analyses_data()
    m_df <- model_inputs(a_df, read_rds(opts$m))
    a_with_imputed_df <- add_imputed_and_msi(a_df, m_df)
    write_rds(m_df, "../output/model_inputs.rds")
    write_rds(a_with_imputed_df, "../output/analyses_data_tidy.rds")
}


main(docopt(doc))
