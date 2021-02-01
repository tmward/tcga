#!/usr/bin/env Rscript
"
Cleans and tidies the CSV holding information on all CRC obtained from
TCGA. Will save a serialized R object (.rds) for a tibble holding
tidied information with patient characteristics and genetics and
a separate tibble holding information on MSI. It only keeps patient data
when the patient has information on a primary CRC tumor and completed
genetic mutation data.

Usage: tidy_tcga.R [-h] [-t TIDYFN] [-m MSI] [-i CSV]

Options:
    -h          Print this menu and exit.
    -t TIDYFN   Filename to save pt info/genetics [default: ../output/tidied.rds].
    -m MSI      Filename to save pt msi [default: ../output/msi.rds].
    -i CSV      File holding TCGA data [default: ../data/gene.csv].
" -> doc

library(docopt)
suppressPackageStartupMessages(library(tidyverse))

get_crc_data <- function(filename) {
  suppressMessages(read_csv(filename)) %>%
    # we only want info on the primary tumors, not the genetics of their
    # mets and recurrences
    filter(
      Sample_Type == "Primary",
      Primary_Site %in% c("Colon", "Rectum", "Rectosigmoid junction")
    ) %>%
    # add an integer rowid col for ease of joining different tables later
    rowid_to_column(var = "caseid") %>%
    janitor::clean_names()
}


add_location <- function(df) {
  df %>%
    mutate(tumor_location = fct_recode(
      # need to parse_factor to make NA be an explicit level
      parse_factor(patient_primary_tumor_site),
      "Right colon" = "Cecum",
      "Right colon" = "Ascending Colon",
      "Right colon" = "Hepatic Flexure",
      "Transverse colon" = "Transverse Colon",
      "Left colon" = "Sigmoid Colon",
      "Left colon" = "Rectosigmoid Jun",
      "Left colon" = "Descending Colon",
      "Left colon" = "Splenic Flexure",
      "Rectum" = "Rectum"
    ))
}


clean_pt_info <- function(df) {
  df %>%
    select(
      caseid,
      # Patient demographics
      diagnosis_age,
      sex,
      patient_height,
      patient_weight,
      # Patient survival info
      patient_s_vital_status,
      overall_survival_months,
      # tumor site
      tumor_location,
      # TNM
      stage_t = diagnosis_pathologict_1,
      stage_n = diagnosis_pathologicn_1,
      stage_m = diagnosis_pathologicm_1,
      # detailed pathology report info
      lns_examined = lymph_node_s_examined_number,
      pni = perineural_invasion,
      lvi = lymphovascular_invasion_indicato,
      vi = vascular_invasion_indicator,
      surgical_margin = surgical_margin_resection_status,
      # detailed genetics info
      fraction_genome_altered
    ) %>%
    # add bmi
    mutate(
      bmi = patient_weight / (patient_height / 100)^2,
      .keep = "unused",
      .after = sex
    ) %>%
    # now clean them up to be usable in a model
    mutate(
      sex = parse_factor(sex),
      patient_s_vital_status = parse_factor(patient_s_vital_status),
      stage_t = case_when(
        stage_t == "T1" ~ 1,
        stage_t == "T2" ~ 2,
        stage_t == "T3" ~ 3,
        stage_t %in% c("T4", "T4a", "T4b") ~ 4,
        stage_t == "Tis" ~ 0,
        is.na(stage_t) ~ NaN,
      ),
      stage_n = case_when(
        stage_n == "N0" ~ 0,
        stage_n %in% c("N1", "N1a", "N1b", "N1c") ~ 1,
        stage_n %in% c("N2", "N2B", "N2a", "N2b") ~ 2,
        stage_n == "NX" ~ NaN,
        is.na(stage_n) ~ NaN
      ),
      stage_m = case_when(
        stage_m %in% c("M0", "Mo") ~ 0,
        stage_m %in% c("M1", "M1a", "M1b") ~ 1,
        stage_m %in% c("MX", "Mx") ~ NaN,
        is.na(stage_m) ~ NaN
      ),
      pni = parse_factor(pni),
      lvi = parse_factor(lvi),
      vi = parse_factor(vi),
      surgical_margin = parse_factor(
        surgical_margin,
        levels = c("R0", "R1", "R2"),
        na = c("RX", NA)
      )
    )
}

fix_info_levels <- function(df) {
  df %>%
    # convert stage_m into a factor column of metastases
    mutate(
      metastases = parse_factor(case_when(
        stage_m == 0 ~ "NO",
        stage_m == 1 ~ "YES",
        is.na(stage_m) ~ NA_character_
      )),
      .keep = "unused"
    ) %>%
    # move NA to the first level to serve as a baseline, but only in
    # factor columns that actually have NA as a level
    mutate(
      across(
        where(
          ~ is.factor(.) && NA %in% levels(.)
        ),
        ~ fct_relevel(fct_explicit_na(.), "(Missing)")
      )
    )
}


add_surv_info <- function(df, yr) {
  m <- yr * 12
  surv_months_x_yr <- str_c("surv_months_", yr, "_yr")
  surv_status_x_yr <- str_c("surv_status_", yr, "_yr")

  df %>%
    mutate(
      {{ surv_months_x_yr }} := if_else(
        overall_survival_months > m, m, overall_survival_months
      ),
      {{ surv_status_x_yr }} := case_when(
        # NA for three possibilites, NA data, NA survival
        # months, or survival is 0 months which is not a valid
        # input for coxnet
        patient_s_vital_status == "(Missing)" |
          is.na(overall_survival_months) |
          overall_survival_months == 0 ~ NA_integer_,
        # Any alive patient is right censored (0) because they never
        # experienced the event
        patient_s_vital_status == "Alive" &
          0 < overall_survival_months ~ 0L,
        # if patient died after m months, they still are
        # technically right-censored since they never had the
        # event during the time period
        patient_s_vital_status == "Dead" &
          m <= overall_survival_months ~ 0L,
        # Failure (1) if died prior to m months
        patient_s_vital_status == "Dead" &
          overall_survival_months < m ~ 1L,
        # test to look for missed items (should be none)
        TRUE ~ -1L
      ),
      # ensure no values fell through otherwise error and stop
      {{ surv_status_x_yr }} := if (-1L %in% .data[[surv_status_x_yr]]) {
        stop("Value not correctly categorized in add_surv_status")
      } else {
        .data[[surv_status_x_yr]]
      }
    )
}

# lump R1 and R2 together
lump_margin <- function(df) {
  df %>%
    mutate(
      pos_margin = fct_recode(surgical_margin,
        "(Missing)" = "(Missing)",
        "NO" = "R0",
        "YES" = "R1",
        "YES" = "R2"
      )
    )
}


patient_data <- function(df) {
  df %>%
    select(!apc_mutation:unc13c_classification, -case) %>%
    add_location() %>%
    clean_pt_info() %>%
    fix_info_levels() %>%
    lump_margin() %>%
    add_surv_info(yr = 1) %>%
    add_surv_info(yr = 2) %>%
    add_surv_info(yr = 3) %>%
    add_surv_info(yr = 5)
}


# Many typos, truncations, etc to standardize classifications across
# genetic mutations
fix_classifications <- function(df) {
  df %>%
    mutate(across(ends_with("_classification"), ~ case_when(
      . == "3_prime_UTR_vari" ~ "3' UTR variant",
      . == "3_prime_UTR_varian" ~ "3' UTR variant",
      . == "3_prime_UTR_variant" ~ "3' UTR variant",
      . == "3_prime_UTR_variant" ~ "3' UTR variant",
      . == "5" ~ "5' UTR variant",
      . == "5_prime_UTR_varian" ~ "5' UTR variant",
      . == "5_prime_UTR_variant" ~ "5' UTR variant",
      . == "i" ~ "Intron variant",
      . == "intron_variant" ~ "Intron variant",
      . == "m" ~ "Missense variant",
      . == "missense_variant" ~ "Missense variant",
      . == "splice_acceptor_" ~ "Splice acceptor variant",
      . == "splice_acceptor_va" ~ "Splice acceptor variant",
      . == "splice_acceptor_vari" ~ "Splice acceptor variant",
      . == "splice_donor_var" ~ "Splice donor variant",
      . == "splice_donor_varia" ~ "Splice donor variant",
      . == "splice_donor_varian" ~ "Splice donor variant",
      . == "splice_donor_variant" ~ "Splice donor variant",
      . == "splice_region_va" ~ "Splice region variant",
      . == "splice_region_vari" ~ "Splice region variant",
      . == "splice_region_varian" ~ "Splice region variant",
      . == "stop_gained" ~ "Stop gained",
      . == "synonymous_varia" ~ "Synonymous variant",
      . == "synonymous_variant" ~ "Synonymous variant",
      . == "s" ~ "Unknown",
      is.na(.) ~ NA_character_,
      TRUE ~ "ERROR"
    )))
}


tidy_genetics <- function(df) {
  df %>%
    pivot_longer(
      -caseid,
      names_to = c("gene", ".value"),
      names_pattern = "([[:alnum:]]+)_+([[:alnum:]]+)"
    ) %>%
    mutate(classification = if_else(
      mutation == 0,
      "Wild type",
      classification
    )) %>%
    # don't need mutation anymore (classification != wild type works)
    select(!mutation)
}


clean_genetics <- function(df) {
  df %>%
    # Remove unknown columns (var 64, 65) and duplicated ryr_class
    select(!var64:ryr_classification_2) %>%
    # fix ttn_cnv which has "-" instead of NA
    mutate(ttn_cnv = as.double(na_if(ttn_cnv, "-"))) %>%
    # fix APC columns that aren't named the same way as the other genes
    rename(apc_classification = apc_class) %>%
    # fix SYNE_Expression which forgets a 1:
    rename(syne1_expression = syne_expression)
}


remove_pts_incomplete_genetics <- function(df) {
  df %>%
    group_by(caseid) %>%
    # add a column that averages the number of na classifications. zero
    # means everything is known. 1 means it's all NA
    mutate(na_class = mean(is.na(classification))) %>%
    ungroup() %>%
    # remove patients with incomplete genetic mutation data
    # could also do na_class != 0 which would gain you one more patient
    # (with an na mucinous gene) but this is most pure
    filter(na_class == 0) %>%
    select(-na_class)
}


# make gene and classification factors
factorize_genes <- function(df) {
  df %>%
    mutate(
      # make "Wild type" be first factor so it will serve as the
      # baseline for dummy variables in the model
      classification = fct_relevel(
        parse_factor(classification), "Wild type"
      ),
      gene = parse_factor(gene)
    )
}


genetics <- function(df) {
  df %>%
    select(caseid, apc_mutation:unc13c_classification) %>%
    clean_genetics() %>%
    fix_classifications() %>%
    tidy_genetics() %>%
    remove_pts_incomplete_genetics() %>%
    factorize_genes()
}

msi <- function(df) {
  df %>%
    select(caseid, pms2_mutations:mlh1_mutations) %>%
    pivot_longer(
      -caseid,
      names_to = "mutation",
      values_to = "mutated"
    ) %>%
    mutate(
      mutated = if_else(
        is.na(mutated),
        0,
        1
      )
    ) %>%
    group_by(caseid) %>%
    summarise(tot_msi = sum(mutated), .groups = "drop") %>%
    mutate(
      msi_pos = as.integer(tot_msi > 0),
      .keep = "unused"
    )
}


main <- function(opts) {
  raw_df <- get_crc_data(opts$i)
  inner_join(genetics(raw_df), patient_data(raw_df), by = "caseid") %>%
    write_rds(opts$t)
  raw_df %>%
    msi() %>%
    write_rds(opts$m)
}

main(docopt(doc))
