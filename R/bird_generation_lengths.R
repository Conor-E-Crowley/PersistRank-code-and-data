# Shared Bird et al. generation-length readers and scientific-name matching.
#
# Inputs are the published workbook plus explicit synonym/curation tables.
# Matching is deterministic and preserves query order; ambiguous or incomplete
# matches fail rather than selecting a candidate silently. Generation length is
# returned in years with match provenance for Stage 1 and Stage 4.

# ---- Workbook normalization --------------------------------------------------

bird_genlength_as_num <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(readr::parse_number(x))
}

quiet_unique_name_repair <- function(x) {
  x <- as.character(x)
  blank <- is.na(x) | !nzchar(x)
  x[blank] <- paste0("...", which(blank))
  make.unique(x, sep = "...")
}

normalize_scientific_name <- function(x) {
  x <- stringr::str_squish(as.character(x))
  tolower(x)
}

read_bird_generation_lengths <- function(path) {
  raw <- readxl::read_xlsx(path, sheet = "Table S4", .name_repair = quiet_unique_name_repair)
  need_cols(raw, c("Scientific name", "GenLength"), "Bird et al. generation-length workbook")

  out <- raw |>
    dplyr::transmute(
      bird_scientific = stringr::str_squish(as.character(`Scientific name`)),
      GenLength = bird_genlength_as_num(GenLength),
      bird_key = normalize_scientific_name(bird_scientific)
    ) |>
    dplyr::filter(nzchar(bird_scientific))

  bad <- out |>
    dplyr::filter(!is.finite(GenLength) | GenLength <= 0)
  if (nrow(bad) > 0) {
    stop(
      "Bird et al. generation-length workbook contains non-positive or missing GenLength values for:\n",
      paste(utils::head(bad$bird_scientific, 25), collapse = "\n"),
      if (nrow(bad) > 25) "\n..." else "",
      call. = FALSE
    )
  }

  dup <- out$bird_scientific[duplicated(out$bird_key)]
  if (length(dup) > 0) {
    stop(
      "Bird et al. generation-length workbook contains duplicate normalized scientific names:\n",
      paste(unique(dup), collapse = "\n"),
      call. = FALSE
    )
  }

  out
}

# ---- Query and diagnostic-table contracts -----------------------------------

# Query order is carried explicitly through every match and diagnostic table.
empty_bird_match_table <- function() {
  tibble::tibble(
    query_order = integer(),
    query_scientific = character(),
    query_key = character(),
    bird_key = character(),
    bird_scientific = character(),
    GenLength = numeric(),
    match_source = character(),
    candidate_name = character(),
    input_synonym_row = integer()
  )
}

empty_bird_diagnostic_table <- function() {
  tibble::tibble(
    query_order = integer(),
    query_scientific = character(),
    query_key = character(),
    match_source = character(),
    reason = character(),
    candidate_name = character(),
    bird_scientific = character(),
    GenLength = numeric(),
    input_synonym_row = integer()
  )
}

bird_generation_queries <- function(scientific_names, label) {
  out <- tibble::tibble(
    query_order = seq_along(scientific_names),
    query_scientific = stringr::str_squish(as.character(scientific_names)),
    query_key = normalize_scientific_name(query_scientific)
  ) |>
    dplyr::filter(nzchar(query_scientific), nzchar(query_key))

  if (nrow(out) != length(scientific_names)) {
    stop(label, " contains missing or blank scientific names.", call. = FALSE)
  }

  dup <- out$query_scientific[duplicated(out$query_key)]
  if (length(dup) > 0) {
    stop(label, " contains duplicate normalized scientific names:\n", paste(unique(dup), collapse = "\n"), call. = FALSE)
  }

  out
}

# ---- Synonym and cosynonym candidate lookup ---------------------------------

# Lookup builders enumerate valid published targets; selection and ambiguity
# handling remain centralized in the matching functions below.
build_iucn_cosynonym_lookup <- function(raw_iucn_bird_synonyms, bird_generation_lengths) {
  group_col <- if ("internalTaxonId" %in% names(raw_iucn_bird_synonyms)) {
    "internalTaxonId"
  } else {
    "scientificName"
  }

  rows <- vector("list", nrow(raw_iucn_bird_synonyms) * 2L)
  k <- 0L
  for (i in seq_len(nrow(raw_iucn_bird_synonyms))) {
    group_id <- as.character(raw_iucn_bird_synonyms[[group_col]][i])
    official <- raw_iucn_bird_synonyms$scientificName[i]
    genus <- raw_iucn_bird_synonyms$genusName[i]
    species <- raw_iucn_bird_synonyms$speciesName[i]
    synonym <- if (is.na(genus) || is.na(species)) NA_character_ else paste(genus, species)

    for (nm in c(official, synonym)) {
      key <- normalize_scientific_name(nm)
      if (is.na(group_id) || !nzchar(group_id) || is.na(key) || !nzchar(key)) next
      k <- k + 1L
      rows[[k]] <- tibble::tibble(group_id = group_id, query_key = key)
    }
  }

  if (k == 0L) return(tibble::tibble(query_key = character(), bird_key = character()))

  long <- dplyr::bind_rows(rows[seq_len(k)]) |>
    dplyr::distinct()
  if (nrow(long) == 0) {
    return(tibble::tibble(query_key = character(), bird_key = character()))
  }

  bird_keys <- unique(bird_generation_lengths$bird_key)
  by_group <- split(long$query_key, long$group_id)
  out <- lapply(by_group, function(keys) {
    keys <- unique(keys)
    targets <- intersect(keys, bird_keys)
    if (length(targets) == 0) return(NULL)
    expand.grid(
      query_key = keys,
      bird_key = targets,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  })

  out <- Filter(Negate(is.null), out)
  if (length(out) == 0) {
    return(tibble::tibble(query_key = character(), bird_key = character()))
  }

  dplyr::bind_rows(out) |>
    dplyr::mutate(candidate_name = bird_key) |>
    dplyr::distinct(query_key, bird_key, candidate_name)
}

build_attached_synonym_lookup <- function(raw_synonyms, bird_generation_lengths) {
  synonym_col <- if ("synonym" %in% names(raw_synonyms)) "synonym" else NULL
  if (is.null(synonym_col) && all(c("genusName", "speciesName") %in% names(raw_synonyms))) {
    raw_synonyms <- raw_synonyms |>
      dplyr::mutate(synonym = stringr::str_squish(paste(genusName, speciesName)))
    synonym_col <- "synonym"
  }
  if (is.null(synonym_col)) {
    stop("Attached synonym matching requires a synonym column or genusName/speciesName columns.", call. = FALSE)
  }

  raw_synonyms |>
    dplyr::transmute(
      query_key = normalize_scientific_name(scientificName),
      candidate_name = stringr::str_squish(as.character(.data[[synonym_col]])),
      bird_key = normalize_scientific_name(candidate_name)
    ) |>
    dplyr::filter(nzchar(query_key), nzchar(bird_key)) |>
    dplyr::inner_join(
      dplyr::select(bird_generation_lengths, bird_key),
      by = "bird_key"
    ) |>
    dplyr::distinct(query_key, bird_key, candidate_name)
}

normalize_input_synonym_table <- function(raw_input_synonyms) {
  out <- raw_input_synonyms
  if (!"input_synonym_row" %in% names(out)) {
    out$input_synonym_row <- seq_len(nrow(out))
  }
  if (!"synonym" %in% names(out)) {
    out$synonym <- out$synonym_scientificName
  }
  out |>
    dplyr::transmute(
      input_synonym_row = as.integer(input_synonym_row),
      query_key = normalize_scientific_name(scientificName),
      candidate_name = stringr::str_squish(as.character(synonym)),
      bird_key = normalize_scientific_name(candidate_name)
    ) |>
    dplyr::filter(nzchar(query_key), nzchar(bird_key))
}

# ---- Deterministic match selection -------------------------------------------

direct_bird_generation_matches <- function(queries, bird_generation_lengths) {
  queries |>
    dplyr::inner_join(
      dplyr::select(bird_generation_lengths, bird_key, bird_scientific, GenLength),
      by = c("query_key" = "bird_key")
    ) |>
    dplyr::transmute(
      query_order, query_scientific, query_key,
      bird_key = query_key,
      bird_scientific, GenLength,
      match_source = "direct",
      candidate_name = query_scientific
    )
}

match_from_lookup <- function(unmatched, lookup, bird_generation_lengths, source_label,
                              query_label, strict_ambiguity = TRUE) {
  if (nrow(unmatched) == 0 || nrow(lookup) == 0) {
    return(list(matches = empty_bird_match_table(), ambiguous = empty_bird_diagnostic_table()))
  }

  candidates <- unmatched |>
    dplyr::select(query_order, query_scientific, query_key) |>
    dplyr::inner_join(lookup, by = "query_key") |>
    dplyr::left_join(
      dplyr::select(bird_generation_lengths, bird_key, bird_scientific, GenLength),
      by = "bird_key"
    ) |>
    dplyr::distinct()

  if (nrow(candidates) == 0) {
    return(list(matches = empty_bird_match_table(), ambiguous = empty_bird_diagnostic_table()))
  }

  target_counts <- candidates |>
    dplyr::distinct(query_scientific, bird_scientific) |>
    dplyr::count(query_scientific, name = "n_targets")
  ambig_names <- target_counts$query_scientific[target_counts$n_targets > 1L]

  if (length(ambig_names) && isTRUE(strict_ambiguity)) {
    details <- lapply(ambig_names, function(sp) {
      targets <- candidates$bird_scientific[candidates$query_scientific == sp]
      paste0("  ", sp, " -> ", paste(sort(unique(targets)), collapse = "; "))
    })
    stop(
      "Ambiguous ", query_label, " ", source_label, " Bird et al. generation-length matches detected:\n",
      paste(unlist(details), collapse = "\n"),
      call. = FALSE
    )
  }

  ambiguous <- candidates |>
    dplyr::filter(query_scientific %in% ambig_names) |>
    dplyr::mutate(match_source = source_label, reason = "multiple_final_targets") |>
    dplyr::select(
      query_order, query_scientific, query_key, match_source, reason,
      candidate_name, bird_scientific, GenLength,
      dplyr::any_of("input_synonym_row")
    )
  if (!"input_synonym_row" %in% names(ambiguous)) ambiguous$input_synonym_row <- NA_integer_

  matches <- candidates |>
    dplyr::filter(!query_scientific %in% ambig_names) |>
    dplyr::group_by(query_order, query_scientific, query_key) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(match_source = source_label) |>
    dplyr::select(
      query_order, query_scientific, query_key,
      bird_key, bird_scientific, GenLength, match_source, candidate_name
    )

  list(matches = matches, ambiguous = ambiguous)
}

match_from_input_lookup <- function(unmatched, lookup, bird_generation_lengths) {
  if (nrow(unmatched) == 0 || nrow(lookup) == 0) {
    return(list(
      matches = empty_bird_match_table(),
      manual_excluded = empty_bird_diagnostic_table()
    ))
  }

  candidates <- unmatched |>
    dplyr::select(query_order, query_scientific, query_key) |>
    dplyr::inner_join(lookup, by = "query_key") |>
    dplyr::left_join(
      dplyr::select(bird_generation_lengths, bird_key, bird_scientific, GenLength),
      by = "bird_key"
    ) |>
    dplyr::distinct()

  if (nrow(candidates) == 0) {
    return(list(
      matches = empty_bird_match_table(),
      manual_excluded = empty_bird_diagnostic_table()
    ))
  }

  rows <- split(candidates, candidates$query_key)
  matches <- list()
  excluded <- list()
  for (key in names(rows)) {
    x <- rows[[key]]
    valid <- x[!is.na(x$bird_scientific) & is.finite(x$GenLength) & x$GenLength > 0, , drop = FALSE]
    reason <- NULL
    if (dplyr::n_distinct(x$input_synonym_row) > 1L) {
      reason <- "multiple_input_synonym_rows"
    } else if (!nrow(valid)) {
      reason <- "no_valid_input_synonym_target"
    } else if (dplyr::n_distinct(valid$bird_key) > 1L) {
      reason <- "multiple_final_targets"
    }

    if (is.null(reason)) {
      matches[[length(matches) + 1L]] <- valid |>
        dplyr::slice(1L) |>
        dplyr::mutate(match_source = "input_synonym") |>
        dplyr::select(
          query_order, query_scientific, query_key,
          bird_key, bird_scientific, GenLength, match_source, candidate_name,
          input_synonym_row
        )
    } else {
      excluded[[length(excluded) + 1L]] <- x |>
        dplyr::mutate(match_source = "input_synonym", reason = reason) |>
        dplyr::select(
          query_order, query_scientific, query_key, match_source, reason,
          candidate_name, bird_scientific, GenLength, input_synonym_row
        )
    }
  }

  list(
    matches = if (length(matches)) dplyr::bind_rows(matches) else empty_bird_match_table(),
    manual_excluded = if (length(excluded)) dplyr::bind_rows(excluded) else empty_bird_diagnostic_table()
  )
}

report_curated_input_synonym_matches <- function(matches, process_label, value_label) {
  need_cols(
    matches,
    c("input_synonym_row", "input_name", "curated_synonym", "matched_name", "matched_value"),
    paste0(process_label, " curated input-synonym audit")
  )

  matches <- dplyr::arrange(matches, input_synonym_row, input_name)
  message(
    "input_synonyms.csv selected matches | process=", process_label,
    " | n=", nrow(matches)
  )
  if (nrow(matches) == 0L) {
    message("  none used")
    return(invisible(matches))
  }

  message(paste(
    sprintf(
      "  row=%d | input=%s | curated=%s | matched=%s | %s=%s",
      matches$input_synonym_row,
      matches$input_name,
      matches$curated_synonym,
      matches$matched_name,
      value_label,
      matches$matched_value
    ),
    collapse = "\n"
  ))
  invisible(matches)
}

match_bird_generation_lengths <- function(scientific_names, bird_generation_lengths,
                                          raw_synonyms, raw_input_synonyms,
                                          query_label,
                                          synonym_mode = c("cosynonym", "attached"),
                                          strict_ambiguity = TRUE) {
  # Precedence is direct name, synonym/cosynonym, then curated input synonym.
  # The curated table may resolve an earlier ambiguity; unresolved, ambiguous,
  # and explicitly excluded queries remain separate and preserve query order.
  synonym_mode <- match.arg(synonym_mode)
  queries <- bird_generation_queries(scientific_names, query_label)

  direct <- direct_bird_generation_matches(queries, bird_generation_lengths)
  unresolved <- queries |>
    dplyr::anti_join(direct, by = "query_key")

  synonym_lookup <- if (identical(synonym_mode, "attached")) {
    build_attached_synonym_lookup(raw_synonyms, bird_generation_lengths)
  } else {
    build_iucn_cosynonym_lookup(raw_synonyms, bird_generation_lengths)
  }
  synonym_label <- if (identical(synonym_mode, "attached")) "synonym" else "iucn"
  synonym_result <- match_from_lookup(
    unresolved,
    synonym_lookup,
    bird_generation_lengths,
    synonym_label,
    query_label,
    strict_ambiguity = strict_ambiguity
  )

  unresolved_after_synonyms <- unresolved |>
    dplyr::anti_join(synonym_result$matches, by = "query_key")

  input_lookup <- normalize_input_synonym_table(raw_input_synonyms)
  input_result <- match_from_input_lookup(
    unresolved_after_synonyms,
    input_lookup,
    bird_generation_lengths
  )

  matches <- dplyr::bind_rows(direct, synonym_result$matches, input_result$matches) |>
    dplyr::arrange(query_order)
  unmatched <- unresolved_after_synonyms |>
    dplyr::anti_join(input_result$matches, by = "query_key") |>
    dplyr::anti_join(input_result$manual_excluded, by = "query_key") |>
    dplyr::anti_join(synonym_result$ambiguous, by = "query_key") |>
    dplyr::arrange(query_order)

  manual_resolved <- input_result$matches |>
    dplyr::semi_join(synonym_result$ambiguous, by = "query_key") |>
    dplyr::left_join(
      synonym_result$ambiguous |>
        dplyr::group_by(query_key) |>
        dplyr::summarise(
          ambiguous_targets = paste(
            sort(unique(paste0(bird_scientific, " (GenLength=", signif(GenLength, 6), ")"))),
            collapse = "; "
          ),
          .groups = "drop"
        ),
      by = "query_key"
    )

  list(
    queries = queries,
    matches = matches,
    unmatched = unmatched,
    ambiguous = synonym_result$ambiguous |>
      dplyr::anti_join(input_result$matches, by = "query_key"),
    manual_used = input_result$matches,
    manual_resolved = manual_resolved,
    manual_excluded = input_result$manual_excluded,
    candidates = dplyr::bind_rows(
      synonym_result$ambiguous,
      input_result$manual_excluded
    )
  )
}

# ---- Human-readable matching audit ------------------------------------------

report_bird_generation_matches <- function(match_result, label,
                                           include_all_matches = FALSE,
                                           verbose = FALSE) {
  matches <- match_result$matches
  unmatched <- match_result$unmatched
  ambiguous <- match_result$ambiguous
  manual_used <- match_result$manual_used
  manual_resolved <- match_result$manual_resolved
  manual_excluded <- match_result$manual_excluded

  counts <- matches |>
    dplyr::count(match_source, name = "n") |>
    dplyr::arrange(match_source)
  count_text <- if (nrow(counts) > 0) {
    paste(paste0(counts$match_source, "=", counts$n), collapse = ", ")
  } else {
    "none"
  }

  message(
    "\nBird et al. GenLength matching for ", label, "\n",
    "  selected: ", count_text, "\n",
    "  unmatched: ", nrow(unmatched), "\n",
    "  ambiguous: ", dplyr::n_distinct(ambiguous$query_key), "\n",
    "  input_synonyms.csv excluded: ", dplyr::n_distinct(manual_excluded$query_key)
  )

  curated_audit <- manual_used |>
    dplyr::transmute(
      input_synonym_row,
      input_name = query_scientific,
      curated_synonym = candidate_name,
      matched_name = bird_scientific,
      matched_value = format(signif(GenLength, 6), scientific = FALSE, trim = TRUE)
    )
  report_curated_input_synonym_matches(curated_audit, label, "GenLength")

  detailed_matches <- if (isTRUE(include_all_matches)) {
    matches
  } else {
    dplyr::filter(matches, !match_source %in% c("direct", "input_synonym"))
  }
  if (isTRUE(verbose) && nrow(detailed_matches) > 0L) {
    message(
      paste(
        sprintf(
          "  %s -> %s | GenLength=%s | source=%s",
          detailed_matches$query_scientific,
          detailed_matches$bird_scientific,
          format(signif(detailed_matches$GenLength, 6), scientific = FALSE, trim = TRUE),
          detailed_matches$match_source
        ),
        collapse = "\n"
      )
    )
  }

  if (isTRUE(verbose) && nrow(manual_resolved) > 0) {
    message(
      "Curated input-synonym GenLength matches that resolved ambiguity for ", label, ":\n",
      paste(
        sprintf(
          "  %s | ambiguous: %s | selected: %s",
          manual_resolved$query_scientific,
          manual_resolved$ambiguous_targets,
          manual_resolved$bird_scientific
        ),
        collapse = "\n"
      )
    )
  }

  if (nrow(ambiguous) > 0) {
    message(
      "Ambiguous GenLength matches excluded for ", label, ":\n",
      paste(
        sprintf(
          "  %s -> %s",
          ambiguous$query_scientific,
          ambiguous$bird_scientific
        ),
        collapse = "\n"
      )
    )
  }

  if (isTRUE(verbose) && nrow(manual_excluded) > 0) {
    message(
      "Curated input-synonym GenLength fallbacks excluded for ", label, ":\n",
      paste(
        sprintf(
          "  %s | reason=%s | row=%s | candidate=%s -> %s",
          manual_excluded$query_scientific,
          manual_excluded$reason,
          manual_excluded$input_synonym_row,
          manual_excluded$candidate_name,
          manual_excluded$bird_scientific
        ),
        collapse = "\n"
      )
    )
  }

  if (nrow(unmatched) > 0) {
    message(
      "Unmatched ", label, " for Bird et al. GenLength:\n",
      paste(sprintf("  %s", unmatched$query_scientific), collapse = "\n")
    )
  }

  invisible(match_result)
}
