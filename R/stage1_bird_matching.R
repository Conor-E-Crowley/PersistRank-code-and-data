# Deterministic bird-trait matching for Stage 1 calibration.
#
# Loaded by stage1_workflow.R before stage1_data.R. It requires shared name and
# bird-model contracts. Sourcing has no side effects; input reads occur only
# when matching helpers are called.

as_num <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(readr::parse_number(x))
}

validate_niel_rmax_table <- function(niel_bird_rmax) {
  need_cols(niel_bird_rmax, c("Species", "lambda"), "Niel & Lebreton bird growth-rate table")

  out <- niel_bird_rmax |>
    dplyr::mutate(
      Species = stringr::str_squish(as.character(Species)),
      lambda = as_num(lambda)
    )

  bad_species <- out$Species[
    is.na(out$Species) | !nzchar(out$Species) |
      !is.finite(out$lambda) | out$lambda <= 1
  ]

  if (length(bad_species) > 0) {
    stop(
      "Niel & Lebreton bird growth records need non-empty Species values and finite lambda values greater than 1. ",
      "Check:\n",
      paste(unique(bad_species), collapse = "\n"),
      call. = FALSE
    )
  }

  if (anyDuplicated(out$Species) > 0) {
    stop("Niel & Lebreton bird growth-rate table contains duplicate Species values.", call. = FALSE)
  }

  out
}

scientific_binomial_key <- function(x) {
  key <- normalize_scientific_name(x)
  parts <- strsplit(key, "\\s+")
  vapply(
    parts,
    function(z) {
      z <- z[nzchar(z)]
      if (length(z) >= 2L) paste(z[1:2], collapse = " ") else paste(z, collapse = " ")
    },
    character(1)
  )
}

read_bird_elton_traits <- function(path) {
  raw <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA", "NaN", "NULL"),
    show_col_types = FALSE,
    progress = FALSE
  )
  need_cols(raw, c("Scientific", "Diet-5Cat"), "EltonTraits bird_data.txt")

  out <- raw |>
    dplyr::transmute(
      elton_scientific = stringr::str_squish(as.character(Scientific)),
      elton_key = normalize_scientific_name(elton_scientific),
      Diet_5Cat = stringr::str_squish(as.character(`Diet-5Cat`))
    ) |>
    dplyr::filter(nzchar(elton_scientific), nzchar(elton_key))

  duplicated_keys <- out$elton_scientific[duplicated(out$elton_key)]
  if (length(duplicated_keys) > 0L) {
    stop(
      "EltonTraits bird_data.txt contains duplicate normalized Scientific names:\n",
      paste(unique(duplicated_keys), collapse = "\n"),
      call. = FALSE
    )
  }

  missing_diet <- out |>
    dplyr::filter(is.na(Diet_5Cat) | !nzchar(Diet_5Cat))
  if (nrow(missing_diet) > 0L) {
    stop(
      "EltonTraits bird_data.txt contains missing Diet-5Cat values for:\n",
      paste(utils::head(missing_diet$elton_scientific, 25), collapse = "\n"),
      if (nrow(missing_diet) > 25L) "\n..." else "",
      call. = FALSE
    )
  }

  out
}

build_elton_iucn_lookup <- function(raw_iucn_bird_synonyms, bird_elton_traits) {
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

  if (k == 0L) {
    return(tibble::tibble(query_key = character(), elton_key = character(), candidate_name = character()))
  }

  long <- dplyr::bind_rows(rows[seq_len(k)]) |>
    dplyr::distinct()
  elton_keys <- unique(bird_elton_traits$elton_key)
  by_group <- split(long$query_key, long$group_id)
  out <- lapply(by_group, function(keys) {
    keys <- unique(keys)
    targets <- intersect(keys, elton_keys)
    if (length(targets) == 0L) return(NULL)
    expand.grid(
      query_key = keys,
      elton_key = targets,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  })

  out <- Filter(Negate(is.null), out)
  if (length(out) == 0L) {
    return(tibble::tibble(query_key = character(), elton_key = character(), candidate_name = character()))
  }

  dplyr::bind_rows(out) |>
    dplyr::left_join(
      dplyr::select(bird_elton_traits, elton_key, elton_scientific),
      by = "elton_key"
    ) |>
    dplyr::transmute(query_key, elton_key, candidate_name = elton_scientific) |>
    dplyr::distinct(query_key, elton_key, candidate_name)
}

build_elton_input_synonym_lookup <- function(raw_input_synonyms) {
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
      candidate_key = normalize_scientific_name(candidate_name),
      candidate_binomial_key = scientific_binomial_key(candidate_name)
    ) |>
    dplyr::filter(nzchar(query_key), nzchar(candidate_key))
}

# Match each Brook et al. bird record through a deterministic precedence:
# exact EltonTraits name, binomial, IUCN synonym group, then curated input
# synonym. Ambiguous or unmatched records stop the fit rather than being
# silently assigned a diet.
select_unique_elton_matches <- function(candidates, source_label) {
  if (nrow(candidates) == 0L) return(candidates)

  target_counts <- candidates |>
    dplyr::distinct(query_scientific, elton_key) |>
    dplyr::count(query_scientific, name = "n_targets")
  ambiguous <- target_counts$query_scientific[target_counts$n_targets > 1L]

  if (length(ambiguous) > 0L) {
    details <- lapply(ambiguous, function(sp) {
      targets <- candidates$elton_scientific[candidates$query_scientific == sp]
      paste0("  ", sp, " -> ", paste(sort(unique(targets)), collapse = "; "))
    })
    stop(
      "Ambiguous EltonTraits ", source_label, " matches for bird sigma calibration:\n",
      paste(unlist(details), collapse = "\n"),
      call. = FALSE
    )
  }

  candidates |>
    dplyr::group_by(query_order, query_scientific, query_key) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(elton_match_source = source_label) |>
    dplyr::select(
      query_order, query_scientific, query_key,
      elton_scientific, Diet_5Cat, elton_match_source,
      dplyr::any_of(c("input_synonym_row", "candidate_name"))
    )
}

# ---- Calibration-table assembly ----------------------------------------------

match_bird_sigma_elton_traits <- function(bird_sigma_data, bird_elton_traits,
                                          raw_iucn_bird_synonyms,
                                          raw_input_synonyms,
                                          verbose = FALSE) {
  queries <- bird_sigma_data |>
    dplyr::transmute(
      query_order = dplyr::row_number(),
      query_scientific = sigma_scientific,
      query_key = normalize_scientific_name(sigma_scientific),
      query_binomial_key = scientific_binomial_key(sigma_scientific)
    )

  direct <- queries |>
    dplyr::inner_join(
      bird_elton_traits,
      by = c("query_key" = "elton_key")
    ) |>
    dplyr::mutate(elton_key = query_key) |>
    select_unique_elton_matches("direct")

  unresolved <- queries |>
    dplyr::anti_join(direct, by = "query_key")

  binomial <- unresolved |>
    dplyr::inner_join(
      bird_elton_traits,
      by = c("query_binomial_key" = "elton_key")
    ) |>
    dplyr::mutate(elton_key = query_binomial_key) |>
    select_unique_elton_matches("binomial")

  unresolved <- unresolved |>
    dplyr::anti_join(binomial, by = "query_key")

  iucn_lookup <- build_elton_iucn_lookup(raw_iucn_bird_synonyms, bird_elton_traits)
  iucn <- unresolved |>
    dplyr::inner_join(iucn_lookup, by = "query_key") |>
    dplyr::left_join(
      bird_elton_traits,
      by = "elton_key"
    ) |>
    select_unique_elton_matches("iucn")

  unresolved <- unresolved |>
    dplyr::anti_join(iucn, by = "query_key")

  input_lookup <- build_elton_input_synonym_lookup(raw_input_synonyms)
  input_direct <- unresolved |>
    dplyr::inner_join(input_lookup, by = "query_key") |>
    dplyr::inner_join(
      bird_elton_traits,
      by = c("candidate_key" = "elton_key")
    ) |>
    dplyr::mutate(elton_key = candidate_key) |>
    select_unique_elton_matches("input_synonym")

  unresolved <- unresolved |>
    dplyr::anti_join(input_direct, by = "query_key")

  input_binomial <- unresolved |>
    dplyr::inner_join(input_lookup, by = "query_key") |>
    dplyr::inner_join(
      bird_elton_traits,
      by = c("candidate_binomial_key" = "elton_key")
    ) |>
    dplyr::mutate(elton_key = candidate_binomial_key) |>
    select_unique_elton_matches("input_synonym_binomial")

  matches <- dplyr::bind_rows(direct, binomial, iucn, input_direct, input_binomial) |>
    dplyr::arrange(query_order)
  unmatched <- queries |>
    dplyr::anti_join(matches, by = "query_key")

  if (nrow(unmatched) > 0L) {
    stop(
      "Every Brook bird sigma calibration record must match EltonTraits. Unmatched records:\n",
      paste(sprintf("  %s", unmatched$query_scientific), collapse = "\n"),
      call. = FALSE
    )
  }

  counts <- matches |>
    dplyr::count(elton_match_source, name = "n") |>
    dplyr::arrange(elton_match_source)
  message(
    "\nEltonTraits Diet-5Cat matching for Brook bird sigma records\n",
    "  selected: ",
    paste(paste0(counts$elton_match_source, "=", counts$n), collapse = ", ")
  )

  input_synonym_used <- matches |>
    dplyr::filter(.data$elton_match_source %in% c("input_synonym", "input_synonym_binomial"))
  curated_audit <- input_synonym_used |>
    dplyr::transmute(
      input_synonym_row,
      input_name = query_scientific,
      curated_synonym = candidate_name,
      matched_name = elton_scientific,
      matched_value = Diet_5Cat
    )
  report_curated_input_synonym_matches(
    curated_audit,
    "Brook bird sigma records to EltonTraits Diet-5Cat",
    "Diet-5Cat"
  )

  if (isTRUE(verbose)) {
    additional_matches <- matches |>
      dplyr::filter(!elton_match_source %in% c("direct", "input_synonym", "input_synonym_binomial"))
    if (nrow(additional_matches) > 0L) {
      message(
        "Additional non-direct EltonTraits matches for Brook bird sigma records:\n",
        paste(
          sprintf(
            "  %s -> %s | Diet-5Cat=%s | source=%s",
            additional_matches$query_scientific,
            additional_matches$elton_scientific,
            additional_matches$Diet_5Cat,
            additional_matches$elton_match_source
          ),
          collapse = "\n"
        )
      )
    }
  }

  matches
}

attach_bird_sigma_model_groups <- function(
  bird_sigma_data,
  bird_elton_traits,
  raw_iucn_bird_synonyms,
  raw_input_synonyms,
  bird_model,
  verbose = FALSE
) {
  matches <- match_bird_sigma_elton_traits(
    bird_sigma_data,
    bird_elton_traits,
    raw_iucn_bird_synonyms,
    raw_input_synonyms,
    verbose = verbose
  )

  out <- bird_sigma_data |>
    dplyr::left_join(
      matches |>
        dplyr::transmute(
          sigma_scientific = query_scientific,
          elton_scientific,
          Diet_5Cat,
          elton_match_source
        ),
      by = "sigma_scientific"
    ) |>
    dplyr::mutate(
      diet5_group = validate_bird_diet5(Diet_5Cat),
      bird_sigma_model_group =
        map_bird_diet_to_model_group(diet5_group, bird_model)
    )

  validate_bird_sigma_diet_matches(out, bird_model)
  out
}

validate_bird_sigma_diet_matches <- function(bird_sigma_data, bird_model) {
  need_cols(
    bird_sigma_data,
    c(
      "sigma_scientific", "elton_scientific", "Diet_5Cat",
      "diet5_group", "bird_sigma_model_group"
    ),
    "bird sigma calibration data"
  )

  missing_match <- is.na(bird_sigma_data$elton_scientific) | !nzchar(bird_sigma_data$elton_scientific)
  if (any(missing_match)) {
    stop(
      "Bird sigma calibration records missing EltonTraits matches:\n",
      paste(unique(bird_sigma_data$sigma_scientific[missing_match]), collapse = "\n"),
      call. = FALSE
    )
  }

  validate_bird_diet5(
    bird_sigma_data$diet5_group,
    "bird sigma calibration Diet-5Cat"
  )
  validate_bird_model_group(
    bird_sigma_data$bird_sigma_model_group,
    bird_model,
    "bird sigma calibration model group"
  )
  assert(
    identical(
      bird_sigma_data$bird_sigma_model_group,
      map_bird_diet_to_model_group(bird_sigma_data$diet5_group, bird_model)
    ),
    "Bird sigma calibration model groups are inconsistent with Diet-5Cat."
  )

  invisible(TRUE)
}

