# Stage 4 IUCN suitable-habitat queries.
#
# Build mode resolves each species from the shared cache, querying only the
# rows required by the selected cache policy. Suitable habitat records are
# normalized to one summary row per species, and cache updates are committed
# transactionally. Figure mode never sources this module.

iucn_level1 <- c(
  `1` = "Forest", `2` = "Savanna", `3` = "Shrubland", `4` = "Grassland",
  `5` = "Wetlands (inland)", `6` = "Rocky Areas", `7` = "Caves & Subterranean (non-aquatic)",
  `8` = "Desert", `9` = "Marine Neritic", `10` = "Marine Oceanic", `11` = "Marine Deep Ocean Floor",
  `12` = "Marine Intertidal", `13` = "Marine Coastal/Supratidal", `14` = "Artificial - Terrestrial",
  `15` = "Artificial - Aquatic", `16` = "Introduced Vegetation", `17` = "Other", `18` = "Unknown"
)

normalize_habitat_code <- function(x) gsub("_", ".", trimws(as.character(x)), fixed = TRUE)

summarize_artificial_terrestrial <- function(codes) {
  codes_14 <- codes[stringr::str_detect(codes, "^14(\\.|$)")]
  if (!length(codes_14)) return(character(0))
  level2 <- suppressWarnings(as.integer(stringr::str_match(codes_14, "^14\\.(\\d+)")[, 2]))
  level2 <- unique(level2[!is.na(level2)])

  out <- character(0)
  if (any(level2 %in% c(1, 2))) out <- c(out, "Arable & Pastureland")
  if (any(level2 %in% c(3, 6))) out <- c(out, "Plantations & Heavily Degraded Former Forest")
  if (any(level2 %in% c(4, 5))) out <- c(out, "Urban & Rural Gardens")
  unique(out)
}

empty_habitat_row <- function() {
  tibble::tibble(
    habitat_codes_suitable = NA_character_,
    habitats_level1 = NA_character_,
    habitats_mixed = NA_character_
  )
}

iucn_api_error_message <- function(genus, species, error) {
  msg <- conditionMessage(error)
  msg_lower <- tolower(msg)
  species_label <- paste(genus, species)

  if (grepl("could not resolve host|could not resolve hostname|name resolution", msg_lower)) {
    return(paste0(
      "IUCN API call failed for ", species_label, " because api.iucnredlist.org could not be resolved.\n",
      "This is a network, DNS, VPN, firewall, or proxy problem rather than a Stage 4 name-matching problem. ",
      "Reconnect to the internet or a network that can resolve api.iucnredlist.org, then rerun Stage 4.\n",
      "Original error: ", msg
    ))
  }

  if (grepl("timeout|timed out|connection.*failed|couldn't connect|ssl|certificate", msg_lower)) {
    return(paste0(
      "IUCN API call failed for ", species_label, " because the API connection failed.\n",
      "Check internet/VPN/proxy/firewall access to api.iucnredlist.org, then rerun Stage 4.\n",
      "Original error: ", msg
    ))
  }

  if (grepl("unauthori[sz]ed|forbidden|api key|401|403", msg_lower)) {
    return(paste0(
      "IUCN API call failed for ", species_label, " because the request was not authorized.\n",
      "Check that IUCN_REDLIST_KEY is set to a valid IUCN Red List API key, then rerun Stage 4.\n",
      "Original error: ", msg
    ))
  }

  paste0("IUCN API call failed for ", species_label, ": ", msg)
}

normalize_iucn_habitats <- function(habitats) {
  if (is.null(habitats) || !nrow(habitats)) return(empty_habitat_row())
  need_cols(habitats, c("suitability", "code"), "IUCN habitat response")

  suitable <- habitats |>
    dplyr::mutate(suitability = stringr::str_to_lower(stringr::str_squish(as.character(suitability)))) |>
    dplyr::filter(suitability == "suitable")

  if (!nrow(suitable)) return(empty_habitat_row())

  codes <- sort(unique(normalize_habitat_code(suitable$code)))
  codes <- codes[!is.na(codes) & nzchar(codes)]
  level1_id <- suppressWarnings(as.integer(stringr::str_extract(codes, "^\\d+")))
  level1_labels <- unique(unname(iucn_level1[as.character(level1_id)]))
  level1_labels <- level1_labels[!is.na(level1_labels)]

  # Specific 14.x records select only their subgroups; the broad parent mask
  # is included only when Level 1 code 14 is explicitly suitable.
  mask_level1_labels <- if ("14" %in% codes) {
    level1_labels
  } else {
    setdiff(level1_labels, unname(iucn_level1[["14"]]))
  }
  mixed_labels <- unique(c(mask_level1_labels, summarize_artificial_terrestrial(codes)))

  tibble::tibble(
    habitat_codes_suitable = if (length(codes)) paste(codes, collapse = ", ") else NA_character_,
    habitats_level1 = if (length(level1_labels)) paste(level1_labels, collapse = ", ") else NA_character_,
    habitats_mixed = if (length(mixed_labels)) paste(mixed_labels, collapse = ", ") else NA_character_
  )
}

query_iucn_habitats <- function(genus, species) {
  res <- tryCatch(
    rredlist::rl_species_latest(genus = genus, species = species, scope = "1", parse = TRUE),
    error = function(e) species_abort(iucn_api_error_message(genus, species, e))
  )
  normalize_iucn_habitats(res$habitats)
}

query_species_habitats <- function(species_inputs, pause_seconds, verbose = FALSE) {
  need_cols(species_inputs, c("scientificName", "genusName", "speciesName"), "Stage 4 IUCN query inputs")
  invalid <- is.na(species_inputs$genusName) |
    !nzchar(stringr::str_squish(as.character(species_inputs$genusName))) |
    is.na(species_inputs$speciesName) |
    !nzchar(stringr::str_squish(as.character(species_inputs$speciesName)))
  if (any(invalid)) {
    labels <- unique(stringr::str_squish(as.character(species_inputs$scientificName[invalid])))
    labels[is.na(labels) | !nzchar(labels)] <- "<missing scientificName>"
    species_abort(
      "Stage 4 cannot query IUCN because retained species have blank genusName or speciesName values:\n",
      paste(utils::head(labels, 25L), collapse = "\n"),
      if (length(labels) > 25L) "\n..." else ""
    )
  }

  query_species <- species_inputs |>
    dplyr::transmute(
      genusName = stringr::str_squish(as.character(genusName)),
      speciesName = stringr::str_squish(as.character(speciesName))
    ) |>
    dplyr::distinct(genusName, speciesName)

  assert(nrow(query_species) > 0L, "No valid genus/species pairs are available for IUCN queries.")
  query_started <- Sys.time()
  report_every <- max(1L, floor(nrow(query_species) / 10L))
  log_msg("Stage 4 | IUCN queries start | species=", nrow(query_species))

  out <- dplyr::bind_rows(Map(
    f = function(genus, species, idx) {
      species_label <- paste(genus, species)

      t0 <- Sys.time()

      out <- query_iucn_habitats(genus, species)
      elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)

      if (isTRUE(verbose) || idx == 1L || idx == nrow(query_species) || idx %% report_every == 0L) {
        log_msg(
          "Stage 4 | IUCN", paste0(idx, "/", nrow(query_species)), "|", species_label,
          "| elapsed=", elapsed, "s",
          if (isTRUE(verbose)) paste0("| suitable_codes=", dplyr::coalesce(out$habitat_codes_suitable[[1]], "none")) else ""
        )
      }

      if (pause_seconds > 0 && idx < nrow(query_species)) {
        Sys.sleep(pause_seconds)
      }

      dplyr::mutate(out, genusName = genus, speciesName = species)
    },
    query_species$genusName,
    query_species$speciesName,
    seq_len(nrow(query_species))
  ))
  assert(
    nrow(out) == nrow(query_species) && !anyDuplicated(out[c("genusName", "speciesName")]),
    "IUCN queries did not return exactly one row per requested genus/species pair."
  )
  log_msg(
    "Stage 4 | IUCN queries complete | species=", nrow(out),
    "| with_suitable_habitat=", sum(!is.na(out$habitats_mixed) & nzchar(out$habitats_mixed)),
    "| elapsed=", round(as.numeric(difftime(Sys.time(), query_started, units = "secs")), 2), "s"
  )
  out
}

iucn_cache_schema <- function() 1L

read_iucn_habitat_cache <- function(path) {
  if (!file.exists(path)) {
    return(list(
      schema = iucn_cache_schema(),
      rows = tibble::tibble(
        genusName = character(), speciesName = character(),
        habitat_codes_suitable = character(), habitats_level1 = character(),
        habitats_mixed = character()
      )
    ))
  }
  cache <- readRDS(path)
  assert(identical(cache$schema, iucn_cache_schema()), "Unsupported IUCN habitat-cache schema.")
  assert(is.data.frame(cache$rows), "IUCN habitat cache rows must be tabular.")
  need_cols(
    cache$rows,
    c("genusName", "speciesName", "habitat_codes_suitable", "habitats_level1", "habitats_mixed"),
    "IUCN habitat cache"
  )
  assert(
    !anyDuplicated(cache$rows[c("genusName", "speciesName")]),
    "IUCN habitat cache contains duplicate species."
  )
  cache
}

write_iucn_habitat_cache <- function(cache, path) {
  ensure_writable_dir(dirname(path), "IUCN habitat-cache directory")
  staged <- tempfile("iucn_cache_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(staged, force = TRUE), add = TRUE)
  saveRDS(cache, staged, version = 3L)
  validated <- read_iucn_habitat_cache(staged)
  assert(identical(validated, cache), "Staged IUCN habitat cache validation failed.")
  project_file_set_transaction(
    staged, path, overwrite = TRUE, label = "IUCN habitat cache"
  )
  invisible(path)
}

resolve_species_habitats <- function(species_inputs,
                                     cache_path,
                                     mode = c("cache_only", "cache_or_query", "refresh"),
                                     pause_seconds = 2,
                                     verbose = FALSE) {
  mode <- match.arg(mode)
  need_cols(species_inputs, c("genusName", "speciesName"), "Stage 4 IUCN inputs")
  requested <- species_inputs |>
    dplyr::transmute(
      genusName = stringr::str_squish(as.character(genusName)),
      speciesName = stringr::str_squish(as.character(speciesName))
    ) |>
    dplyr::distinct(genusName, speciesName)
  cache <- read_iucn_habitat_cache(cache_path)
  keys <- paste(requested$genusName, requested$speciesName, sep = "\r")
  cached_keys <- paste(cache$rows$genusName, cache$rows$speciesName, sep = "\r")
  query_rows <- if (identical(mode, "refresh")) {
    requested
  } else {
    requested[!keys %in% cached_keys, , drop = FALSE]
  }
  if (identical(mode, "cache_only") && nrow(query_rows)) {
    missing <- paste(query_rows$genusName, query_rows$speciesName)
    species_abort(
      "IUCN cache_only mode is missing requested species:\n",
      paste(utils::head(missing, 25L), collapse = "\n"),
      if (length(missing) > 25L) "\n..." else ""
    )
  }
  if (nrow(query_rows)) {
    if (!nzchar(Sys.getenv("IUCN_REDLIST_KEY"))) {
      species_abort(
        "IUCN_REDLIST_KEY is not set and ", nrow(query_rows),
        " requested species are absent from the shared IUCN cache."
      )
    }
    load_packages("rredlist")
    queried <- query_species_habitats(
      dplyr::mutate(
        query_rows,
        scientificName = paste(genusName, speciesName),
        .before = genusName
      ),
      pause_seconds = pause_seconds, verbose = verbose
    )
    replaced_keys <- paste(queried$genusName, queried$speciesName, sep = "\r")
    retained <- cache$rows[!cached_keys %in% replaced_keys, , drop = FALSE]
    cache <- list(
      schema = iucn_cache_schema(),
      rows = dplyr::bind_rows(retained, queried) |>
        dplyr::arrange(genusName, speciesName)
    )
    write_iucn_habitat_cache(cache, cache_path)
  }
  cache_keys <- paste(cache$rows$genusName, cache$rows$speciesName, sep = "\r")
  used <- cache$rows[match(keys, cache_keys), , drop = FALSE]
  assert(nrow(used) == nrow(requested) && !any(is.na(used$genusName)),
         "IUCN habitat cache could not supply every requested species.")
  list(rows = used, queried = nrow(query_rows))
}
