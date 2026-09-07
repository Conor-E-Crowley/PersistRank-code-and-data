# Deterministic Stage 4 scientific-name candidate selection.
#
# Loaded by stage4_workflow.R after input contracts. It owns candidate priority,
# target identity, ambiguity handling, and unique resolution. Sourcing defines
# functions only and performs no directory scan or artifact I/O.

name_candidate_priority <- function(match_type) {
  dplyr::case_when(
    match_type == "original" ~ 0L,
    match_type == "synonym" ~ 1L,
    match_type == "input_synonym" ~ 2L,
    TRUE ~ 99L
  )
}

annotate_name_candidates <- function(candidates, traits) {
  trait_index <- unique(paste(traits$taxon_class, traits$trait_key, sep = "||"))

  candidates |>
    dplyr::filter(!is.na(taxon_class), !is.na(candidate_name), nzchar(candidate_name)) |>
    dplyr::mutate(
      candidate_stem = species_to_stem(candidate_name),
      candidate_key = clean_name(candidate_name),
      name_priority = name_candidate_priority(match_type),
      trait_available = paste(taxon_class, candidate_key, sep = "||") %in% trait_index
    )
}

build_name_candidates <- function(summary_df, synonyms, input_synonyms, traits) {
  base <- dplyr::bind_rows(
    summary_df |>
      dplyr::transmute(
        row_id,
        taxon_class,
        summary_scientific = scientificName,
        candidate_name = scientificName,
        match_type = "original",
        synonym_source_row = NA_integer_,
        input_synonym_row = NA_integer_
      ),
    summary_df |>
      dplyr::select(row_id, taxon_class, scientificName) |>
      dplyr::left_join(synonyms, by = "scientificName") |>
      dplyr::transmute(
        row_id,
        taxon_class,
        summary_scientific = scientificName,
        candidate_name = synonym,
        match_type = "synonym",
        synonym_source_row,
        input_synonym_row = NA_integer_
      )
  )

  if (!is.null(input_synonyms)) {
    input <- summary_df |>
      dplyr::select(row_id, taxon_class, scientificName) |>
      dplyr::left_join(input_synonyms, by = "scientificName") |>
      dplyr::transmute(
        row_id,
        taxon_class,
        summary_scientific = scientificName,
        candidate_name = synonym,
        match_type = "input_synonym",
        synonym_source_row = NA_integer_,
        input_synonym_row
      )
    base <- dplyr::bind_rows(base, input)
  }

  annotate_name_candidates(base, traits)
}

name_resolution_source <- function(x) {
  dplyr::case_when(
    x == "original" ~ "direct",
    x == "synonym" ~ "synonyms.csv",
    x == "input_synonym" ~ "input_synonyms.csv",
    TRUE ~ as.character(x)
  )
}

collapse_unique <- function(x) {
  x <- sort(unique(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  paste(x, collapse = "; ")
}
candidate_has_target <- function(x, target_cols) {
  if (!nrow(x)) return(logical())
  Reduce(`&`, lapply(target_cols, function(col) {
    value <- x[[col]]
    !is.na(value) & nzchar(as.character(value))
  }))
}

candidate_target_key <- function(x, target_cols) {
  if (!nrow(x)) return(character())
  do.call(
    paste,
    c(lapply(target_cols, function(col) as.character(x[[col]])), sep = "\r")
  )
}

candidate_target_label <- function(x, target_cols, detail_cols = character()) {
  if (!nrow(x)) return(character())
  cols <- unique(c("candidate_name", detail_cols, target_cols))
  cols <- cols[cols %in% names(x)]
  apply(
    as.data.frame(x[, cols, drop = FALSE]),
    1,
    function(row) {
      row <- row[!is.na(row) & nzchar(as.character(row))]
      paste(paste(names(row), row, sep = "="), collapse = ", ")
    }
  )
}

empty_match_diagnostics <- function() {
  list(
    excluded = tibble::tibble(
      row_id = integer(), taxon_class = character(), scientificName = character(),
      use_site = character(), reason = character()
    ),
    ambiguous = tibble::tibble(
      row_id = integer(), taxon_class = character(), scientificName = character(),
      use_site = character(), match_type = character(), reason = character(),
      candidate_name = character(), target_label = character()
    ),
    manual_used = tibble::tibble(
      row_id = integer(), taxon_class = character(), scientificName = character(),
      use_site = character(), input_synonym_row = integer(), candidate_name = character(),
      selected_target = character(), resolved_ambiguity = logical()
    ),
    manual_resolved = tibble::tibble(
      row_id = integer(), taxon_class = character(), scientificName = character(),
      use_site = character(), prior_match_type = character(), ambiguous_targets = character(),
      input_synonym_row = integer(), candidate_name = character(), selected_target = character()
    ),
    manual_excluded = tibble::tibble(
      row_id = integer(), taxon_class = character(), scientificName = character(),
      use_site = character(), reason = character(), input_synonym_row = integer(),
      candidate_name = character(), target_label = character()
    )
  )
}

diagnostic_identity <- function(rows, use_site, reason) {
  tibble::tibble(
    row_id = rows$row_id[1],
    taxon_class = rows$taxon_class[1],
    scientificName = rows$summary_scientific[1],
    use_site = use_site,
    reason = reason
  )
}

diagnostic_candidate_rows <- function(rows, use_site, reason, target_cols, detail_cols = character()) {
  if (!nrow(rows)) return(empty_match_diagnostics()$ambiguous)
  rows <- rows |>
    dplyr::distinct(
      row_id, taxon_class, summary_scientific, match_type, candidate_name,
      dplyr::across(dplyr::any_of(c("input_synonym_row", target_cols, detail_cols))),
      .keep_all = TRUE
    )
  rows$target_label <- candidate_target_label(rows, target_cols, detail_cols)
  rows |>
    dplyr::transmute(
      row_id,
      taxon_class,
      scientificName = summary_scientific,
      use_site = use_site,
      match_type,
      reason = reason,
      input_synonym_row = if ("input_synonym_row" %in% names(rows)) input_synonym_row else NA_integer_,
      candidate_name,
      target_label
    )
}

select_first_candidate <- function(rows) {
  rows |>
    dplyr::arrange(
      name_priority,
      candidate_name,
      dplyr::across(dplyr::any_of(c("sdm_method", "raster_file", "trait_key_selected")))
    ) |>
    dplyr::slice(1L)
}

resolve_unique_candidates <- function(matches, target_cols, use_site,
                                      detail_cols = character()) {
  diagnostics <- empty_match_diagnostics()
  selected <- list()

  for (id in unique(matches$row_id)) {
    rows <- matches[matches$row_id == id, , drop = FALSE]
    valid_rows <- rows[candidate_has_target(rows, target_cols), , drop = FALSE]
    picked <- NULL
    prior_ambiguous <- NULL
    prior_match_type <- NA_character_

    for (source in c("original", "synonym")) {
      source_rows <- valid_rows[valid_rows$match_type == source, , drop = FALSE]
      if (!nrow(source_rows)) next

      n_targets <- dplyr::n_distinct(candidate_target_key(source_rows, target_cols))
      if (n_targets == 1L) {
        picked <- select_first_candidate(source_rows)
      } else {
        prior_ambiguous <- source_rows
        prior_match_type <- source
      }
      break
    }

    if (is.null(picked)) {
      manual_all <- rows[rows$match_type == "input_synonym" & !is.na(rows$input_synonym_row), , drop = FALSE]
      manual_valid <- manual_all[candidate_has_target(manual_all, target_cols), , drop = FALSE]

      if (nrow(manual_all)) {
        manual_rows <- unique(manual_all$input_synonym_row)
        if (length(manual_rows) > 1L) {
          diagnostics$manual_excluded <- dplyr::bind_rows(
            diagnostics$manual_excluded,
            diagnostic_candidate_rows(manual_all, use_site, "multiple_input_synonym_rows", target_cols, detail_cols) |>
              dplyr::select(row_id, taxon_class, scientificName, use_site, reason, input_synonym_row, candidate_name, target_label)
          )
        } else if (!nrow(manual_valid)) {
          diagnostics$manual_excluded <- dplyr::bind_rows(
            diagnostics$manual_excluded,
            diagnostic_candidate_rows(manual_all, use_site, "no_valid_manual_target", target_cols, detail_cols) |>
              dplyr::select(row_id, taxon_class, scientificName, use_site, reason, input_synonym_row, candidate_name, target_label)
          )
        } else if (dplyr::n_distinct(candidate_target_key(manual_valid, target_cols)) > 1L) {
          diagnostics$manual_excluded <- dplyr::bind_rows(
            diagnostics$manual_excluded,
            diagnostic_candidate_rows(manual_valid, use_site, "multiple_manual_final_targets", target_cols, detail_cols) |>
              dplyr::select(row_id, taxon_class, scientificName, use_site, reason, input_synonym_row, candidate_name, target_label)
          )
        } else {
          picked <- select_first_candidate(manual_valid)
          selected_target <- candidate_target_label(picked, target_cols, detail_cols)
          resolved_ambiguity <- !is.null(prior_ambiguous)
          diagnostics$manual_used <- dplyr::bind_rows(
            diagnostics$manual_used,
            tibble::tibble(
              row_id = picked$row_id,
              taxon_class = picked$taxon_class,
              scientificName = picked$summary_scientific,
              use_site = use_site,
              input_synonym_row = picked$input_synonym_row,
              candidate_name = picked$candidate_name,
              selected_target = selected_target,
              resolved_ambiguity = resolved_ambiguity
            )
          )
          if (resolved_ambiguity) {
            diagnostics$manual_resolved <- dplyr::bind_rows(
              diagnostics$manual_resolved,
              tibble::tibble(
                row_id = picked$row_id,
                taxon_class = picked$taxon_class,
                scientificName = picked$summary_scientific,
                use_site = use_site,
                prior_match_type = prior_match_type,
                ambiguous_targets = collapse_unique(candidate_target_label(prior_ambiguous, target_cols, detail_cols)),
                input_synonym_row = picked$input_synonym_row,
                candidate_name = picked$candidate_name,
                selected_target = selected_target
              )
            )
          }
        }
      }
    }

    if (!is.null(picked)) {
      selected[[length(selected) + 1L]] <- picked
      next
    }

    if (!is.null(prior_ambiguous)) {
      diagnostics$ambiguous <- dplyr::bind_rows(
        diagnostics$ambiguous,
        diagnostic_candidate_rows(
          prior_ambiguous,
          use_site,
          paste0("ambiguous_", prior_match_type, "_match"),
          target_cols,
          detail_cols
        )
      )
      diagnostics$excluded <- dplyr::bind_rows(
        diagnostics$excluded,
        diagnostic_identity(rows, use_site, paste0("ambiguous_", prior_match_type, "_match"))
      )
    } else if (nrow(rows[rows$match_type == "input_synonym" & !is.na(rows$input_synonym_row), , drop = FALSE])) {
      diagnostics$excluded <- dplyr::bind_rows(
        diagnostics$excluded,
        diagnostic_identity(rows, use_site, "input_synonym_fallback_unresolved")
      )
    } else {
      diagnostics$excluded <- dplyr::bind_rows(
        diagnostics$excluded,
        diagnostic_identity(rows, use_site, "no_match")
      )
    }
  }

  out <- if (length(selected)) dplyr::bind_rows(selected) else matches[0, , drop = FALSE]
  list(selected = out, diagnostics = diagnostics)
}

