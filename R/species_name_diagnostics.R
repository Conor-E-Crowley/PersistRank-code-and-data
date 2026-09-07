# Stage 4 name-resolution and exclusion diagnostics.
#
# Loaded by stage4_workflow.R after candidate selection. It formats and reports
# already-computed matching decisions without changing them. Sourcing has no
# side effects; console messages are emitted only when reporters are called.

format_table_rows <- function(x) {
  if (!nrow(x)) return("  (none)")

  x <- as.data.frame(x, stringsAsFactors = FALSE)
  values <- as.data.frame(
    lapply(x, function(col) {
      out <- trimws(as.character(col))
      out[is.na(out) | !nzchar(out)] <- "-"
      out
    }),
    stringsAsFactors = FALSE
  )
  widths <- pmax(
    nchar(names(values)),
    vapply(values, function(col) as.integer(max(nchar(col))), integer(1))
  )
  format_one <- function(row) {
    paste0("  ", paste(sprintf(paste0("%-", widths, "s"), row), collapse = "  "))
  }

  c(format_one(names(values)), apply(values, 1, format_one))
}

format_detail_rows <- function(x) {
  if (!nrow(x)) return(character())

  x <- as.data.frame(x, stringsAsFactors = FALSE)
  apply(
    x,
    1,
    function(row) {
      values <- trimws(as.character(row))
      keep <- !is.na(values) & nzchar(values)
      if (!any(keep)) return("  -")
      paste0("  - ", paste(paste(names(row)[keep], values[keep], sep = ": "), collapse = " | "))
    }
  )
}

message_stage4_section <- function(title) {
  message("\n", title)
}

message_count_table <- function(title, x) {
  message(title, ":\n", paste(format_table_rows(x), collapse = "\n"))
}
message_detail_table <- function(title, x, cols) {
  if (!nrow(x)) return(invisible(FALSE))
  cols <- cols[cols %in% names(x)]
  message(title, ":\n", paste(format_detail_rows(x[, cols, drop = FALSE]), collapse = "\n"))
  invisible(TRUE)
}

report_match_diagnostics <- function(label, diagnostics, verbose = FALSE) {
  message_stage4_section(paste0(label, " curated input-synonym matches"))
  if (!nrow(diagnostics$manual_used)) {
    message("  none used")
  } else {
    message_detail_table(
      "Selected mappings",
      diagnostics$manual_used,
      c("input_synonym_row", "scientificName", "candidate_name", "selected_target")
    )
  }
  if (!isTRUE(verbose)) return(invisible(TRUE))

  message_detail_table(
    paste0(label, " ambiguous matches excluded"),
    diagnostics$ambiguous,
    c("scientificName", "match_type", "candidate_name", "target_label")
  )
  message_detail_table(
    paste0(label, " matches using input_synonyms.csv"),
    diagnostics$manual_used,
    c("scientificName", "input_synonym_row", "candidate_name", "selected_target", "resolved_ambiguity")
  )
  message_detail_table(
    paste0(label, " ambiguities resolved by input_synonyms.csv"),
    diagnostics$manual_resolved,
    c("scientificName", "prior_match_type", "ambiguous_targets", "input_synonym_row", "candidate_name", "selected_target")
  )
  message_detail_table(
    paste0(label, " input_synonyms.csv fallbacks excluded"),
    diagnostics$manual_excluded,
    c("scientificName", "reason", "input_synonym_row", "candidate_name", "target_label")
  )
}

empty_sdm_matched_exclusion_table <- function() {
  tibble::tibble(
    row_id = integer(),
    scientificName = character(),
    taxon_class = character(),
    sdm_method = character(),
    raster_file = character(),
    name_raster = character(),
    excluded_at = character(),
    reason = character(),
    detail = character()
  )
}

match_diagnostic_detail_lookup <- function(diagnostics) {
  out <- dplyr::bind_rows(
    diagnostics$ambiguous |>
      dplyr::transmute(
        row_id,
        detail = paste0("ambiguous ", match_type, ": ", candidate_name, " -> ", target_label)
      ),
    diagnostics$manual_excluded |>
      dplyr::transmute(
        row_id,
        detail = paste0("input_synonyms.csv ", reason, ": ", candidate_name, " -> ", target_label)
      )
  )

  if (!nrow(out)) {
    return(tibble::tibble(row_id = integer(), detail = character()))
  }

  out |>
    dplyr::group_by(row_id) |>
    dplyr::summarise(detail = collapse_unique(detail), .groups = "drop")
}

sdm_matched_exclusions <- function(excluded, raster_pick, excluded_at,
                                   detail_lookup = NULL) {
  if (!nrow(excluded) || !nrow(raster_pick)) return(empty_sdm_matched_exclusion_table())

  out <- excluded |>
    dplyr::inner_join(
      raster_pick |>
        dplyr::select(row_id, sdm_method, raster_file, name_raster),
      by = "row_id"
    )

  if (!is.null(detail_lookup) && nrow(detail_lookup)) {
    out <- out |>
      dplyr::left_join(detail_lookup, by = "row_id")
  }
  if (!"detail" %in% names(out)) out$detail <- NA_character_

  out |>
    dplyr::transmute(
      row_id,
      scientificName,
      taxon_class,
      sdm_method,
      raster_file,
      name_raster,
      excluded_at = excluded_at,
      reason,
      detail = dplyr::coalesce(detail, "")
    ) |>
    dplyr::arrange(excluded_at, taxon_class, sdm_method, scientificName)
}

message_sdm_matched_exclusions <- function(exclusions, verbose = FALSE) {
  message_stage4_section("Exclusions after SDM raster match")
  counts <- exclusions |>
    dplyr::count(excluded_at, taxon_class, sdm_method, reason, name = "n") |>
    dplyr::arrange(excluded_at, taxon_class, sdm_method, reason)
  message_count_table("By exclusion stage, taxon, SDM source, and reason", counts)
  if (isTRUE(verbose)) {
    message_detail_table(
      "SDM-matched species excluded after raster selection",
      exclusions,
      c("scientificName", "taxon_class", "sdm_method", "raster_file", "name_raster", "excluded_at", "reason", "detail")
    )
  }
}

report_name_resolution <- function(summary_df, raster_pick, species_inputs,
                                   raster_diagnostics, trait_diagnostics,
                                   verbose = FALSE) {
  taxon_lookup <- summary_df |>
    dplyr::filter(taxon_class %in% c("Mammalia", "Aves")) |>
    dplyr::select(row_id, taxon_class)

  input_counts <- taxon_lookup |>
    dplyr::count(taxon_class, name = "n_input") |>
    dplyr::arrange(taxon_class)

  raster_counts <- raster_pick |>
    dplyr::left_join(taxon_lookup, by = "row_id") |>
    dplyr::mutate(name_source = name_resolution_source(match_raster)) |>
    dplyr::count(taxon_class, sdm_method, name_source, name = "n") |>
    dplyr::arrange(taxon_class, sdm_method, name_source)

  trait_counts <- species_inputs |>
    dplyr::mutate(name_source = name_resolution_source(match_trait)) |>
    dplyr::count(taxon_class, sdm_method, name_source, name = "n") |>
    dplyr::arrange(taxon_class, sdm_method, name_source)

  retained_counts <- species_inputs |>
    dplyr::count(taxon_class, sdm_method, name = "n") |>
    dplyr::arrange(taxon_class, sdm_method)

  trait_excluded <- trait_diagnostics$excluded |>
    dplyr::filter(row_id %in% raster_pick$row_id)
  trait_excluded_after_raster <- sdm_matched_exclusions(
    trait_excluded,
    raster_pick,
    excluded_at = "EltonTraits",
    detail_lookup = match_diagnostic_detail_lookup(trait_diagnostics)
  )
  excluded_counts <- dplyr::bind_rows(
    raster_diagnostics$excluded,
    trait_excluded
  ) |>
    dplyr::count(use_site, taxon_class, reason, name = "n") |>
    dplyr::arrange(use_site, taxon_class, reason)

  message("\nStage 4 name-resolution report")
  message_stage4_section("Inputs")
  message_count_table("Input species by taxon", input_counts)
  message_stage4_section("Selected SDM raster matches")
  message_count_table("By taxon, SDM source, and name source", raster_counts)
  message_stage4_section("Selected EltonTraits matches")
  message_count_table("Among raster-matched species, by taxon, SDM source, and name source", trait_counts)
  message_count_table("Rows retained after SDM raster and EltonTraits filtering", retained_counts)
  message_count_table("Species excluded before GenLength by stage and reason", excluded_counts)
  message_detail_table(
    "SDM-matched species excluded before GenLength",
    trait_excluded_after_raster,
    c("scientificName", "taxon_class", "sdm_method", "reason")
  )

  report_match_diagnostics("SDM raster", raster_diagnostics, verbose = verbose)
  report_match_diagnostics("EltonTraits", trait_diagnostics, verbose = verbose)

  invisible(trait_excluded_after_raster)
}
bird_genlength_detail_lookup <- function(match_result) {
  ambiguous <- match_result$ambiguous
  if (nrow(ambiguous)) {
    gen <- ifelse(
      is.finite(ambiguous$GenLength),
      format(signif(ambiguous$GenLength, 6), scientific = FALSE, trim = TRUE),
      "NA"
    )
    ambiguous <- ambiguous |>
      dplyr::transmute(
        query_key,
        reason,
        detail = paste0("ambiguous ", match_source, ": ", candidate_name, " -> ", bird_scientific, " (GenLength=", gen, ")")
      )
  } else {
    ambiguous <- tibble::tibble(query_key = character(), reason = character(), detail = character())
  }

  manual_excluded <- match_result$manual_excluded
  if (nrow(manual_excluded)) {
    gen <- ifelse(
      is.finite(manual_excluded$GenLength),
      format(signif(manual_excluded$GenLength, 6), scientific = FALSE, trim = TRUE),
      "NA"
    )
    manual_excluded <- manual_excluded |>
      dplyr::transmute(
        query_key,
        reason,
        detail = paste0(
          "input_synonyms.csv ", reason,
          ": row ", input_synonym_row,
          ", ", candidate_name, " -> ", bird_scientific,
          " (GenLength=", gen, ")"
        )
      )
  } else {
    manual_excluded <- tibble::tibble(query_key = character(), reason = character(), detail = character())
  }

  unmatched <- match_result$unmatched
  if (nrow(unmatched)) {
    unmatched <- unmatched |>
      dplyr::transmute(
        query_key,
        reason = "no_match",
        detail = "no Bird et al. GenLength target selected"
      )
  } else {
    unmatched <- tibble::tibble(query_key = character(), reason = character(), detail = character())
  }

  out <- dplyr::bind_rows(unmatched, ambiguous, manual_excluded)
  if (!nrow(out)) return(tibble::tibble(query_key = character(), reason = character(), detail = character()))

  out |>
    dplyr::group_by(query_key) |>
    dplyr::summarise(
      reason = collapse_unique(reason),
      detail = collapse_unique(detail),
      .groups = "drop"
    )
}

genlength_exclusions_after_sdm <- function(rows, match_result) {
  if (!nrow(rows)) return(empty_sdm_matched_exclusion_table())

  detail_lookup <- bird_genlength_detail_lookup(match_result)
  out <- rows |>
    dplyr::mutate(query_key = normalize_scientific_name(scientificName))
  if (nrow(detail_lookup)) {
    out <- out |>
      dplyr::left_join(detail_lookup, by = "query_key")
  }
  if (!"reason" %in% names(out)) out$reason <- NA_character_
  if (!"detail" %in% names(out)) out$detail <- NA_character_

  out |>
    dplyr::transmute(
      row_id = NA_integer_,
      scientificName,
      taxon_class,
      sdm_method,
      raster_file,
      name_raster,
      excluded_at = "GenLength",
      reason = dplyr::coalesce(reason, "invalid_or_missing_genlength"),
      detail = dplyr::coalesce(detail, "GenLength missing, non-finite, or non-positive after matching")
    ) |>
    dplyr::arrange(taxon_class, sdm_method, scientificName)
}

report_stage4_final_resolution <- function(species_inputs, exclusions_after_sdm,
                                           genlength_counts, verbose = FALSE) {
  if (is.null(exclusions_after_sdm)) {
    exclusions_after_sdm <- empty_sdm_matched_exclusion_table()
  }

  final_counts <- species_inputs |>
    dplyr::count(taxon_class, sdm_method, name = "n") |>
    dplyr::arrange(taxon_class, sdm_method)

  message_stage4_section("Selected Bird GenLength matches")
  message_count_table("By SDM source and name source", genlength_counts)
  message_stage4_section("Final retained rows")
  message_count_table("By taxon and SDM source", final_counts)
  message_sdm_matched_exclusions(exclusions_after_sdm, verbose = verbose)
}

