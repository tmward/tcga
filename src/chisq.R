#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
library(gtsummary)

suppressMessages(theme_gtsummary_journal(journal = "jama"))
suppressMessages(theme_gtsummary_compact())


get_chisq_input <- function() {
    read_rds("../output/analyses_data_tidy.rds") %>%
        mutate(
            classification = parse_factor(
                classification,
                levels = c("Wild type", "Mutated")
            ),
            tumor_location = parse_factor(
                str_remove(tumor_location, " colon$"),
                levels = c("Right", "Left")
            ),
            nodal_positivity = parse_factor(
                if_else(stage_n == 0, "Negative", "Positive"),
                levels = c("Negative", "Positive")
            )
        )
}


chisq_df <- function(df, variable, by) {
    df %>%
        select(caseid, gene, {{ variable }}, {{ by }}) %>%
        mutate(gene = str_to_upper(gene)) %>%
        # columns will pivot in lexigraphical order now
        arrange(gene) %>%
        pivot_wider(names_from = gene, values_from = {{ variable }}) %>%
        select(-caseid)
}


calc_chisq <- function(data, variable, by, ...) {
    set.seed(1234)
    results <- chisq.test(
        x = data[[variable]],
        y = data[[by]],
        simulate.p.value = TRUE,
        B = 1000000
    )
    list(p = results$p.value, test = results$method)
}


chisq_summary_by_location <- function(df, by_var, span_label) {
    df %>%
        tbl_summary(by = by_var) %>%
        add_p(test = list(all_categorical() ~ calc_chisq)) %>%
        add_q(method = "holm", quiet = TRUE) %>%
        modify_header(update = list(label ~ "**Gene**")) %>%
        modify_spanning_header(c("stat_1", "stat_2") ~ span_label)
}


# Tumor location and mutation table
location_and_mutation <- function() {
    get_chisq_input() %>%
        chisq_df(classification, tumor_location) %>%
        chisq_summary_by_location("tumor_location", "**Tumor Site**")
}


# Tumor location and nodal positivity table
location_and_nodal_positivity <- function() {
    get_chisq_input() %>%
        chisq_df(classification, nodal_positivity) %>%
        chisq_summary_by_location("nodal_positivity", "**Nodal Status**")
}


save_gt <- function(df, filename) {
    df %>%
        as_gt() %>%
        gt::gtsave(filename, path = "../output/")
}


location_and_nodal_positivity() %>%
    save_gt("chisq_location_nodal_positivity.html")

location_and_mutation() %>%
    save_gt("chisq_location_mutated.html")
