# High-level Stage 4 species-name resolution coordination.
#
# Loaded by stage4_workflow.R after candidate and diagnostic definitions. It
# coordinates raster, trait, and bird-generation-length selection. Sourcing is
# definition-only; matching work and reporting occur only on call.

pick_rasters <- function(name_candidates, raster_manifest) {
  raster_matches <- name_candidates |>
    dplyr::left_join(
      dplyr::select(
        raster_manifest,
        taxon_class, raster_stem, sdm_method, raster_dir, raster_path, raster_file
      ),
      by = c("taxon_class", "candidate_stem" = "raster_stem")
    )

  resolved <- resolve_unique_candidates(
    raster_matches,
    target_cols = c("sdm_method", "raster_path"),
    use_site = "SDM raster",
    detail_cols = c("sdm_method", "raster_file")
  )

  out <- resolved$selected |>
    dplyr::transmute(
      row_id,
      name_raster = candidate_name,
      match_raster = match_type,
      raster_stem = candidate_stem,
      sdm_method,
      raster_dir,
      raster_file,
      raster_path
    )
  attr(out, "diagnostics") <- resolved$diagnostics
  out
}

pick_traits <- function(name_candidates) {
  trait_candidates <- name_candidates |>
    dplyr::mutate(
      trait_key_selected = dplyr::if_else(trait_available, candidate_key, NA_character_)
    )

  resolved <- resolve_unique_candidates(
    trait_candidates,
    target_cols = "trait_key_selected",
    use_site = "EltonTraits",
    detail_cols = "candidate_name"
  )

  out <- resolved$selected |>
    dplyr::transmute(
      row_id,
      name_trait = candidate_name,
      match_trait = match_type,
      trait_key = trait_key_selected
    )
  attr(out, "diagnostics") <- resolved$diagnostics
  out
}

resolve_species_inputs <- function(core, verbose = FALSE) {
  name_candidates <- build_name_candidates(
    summary_df = core$summary_df,
    synonyms = core$synonyms,
    input_synonyms = core$input_synonyms,
    traits = core$traits
  )

  raster_pick <- pick_rasters(name_candidates, core$raster_manifest)
  raster_diagnostics <- attr(raster_pick, "diagnostics")
  trait_pick <- pick_traits(name_candidates)
  trait_diagnostics <- attr(trait_pick, "diagnostics")

  species_inputs <- core$summary_df |>
    dplyr::left_join(raster_pick, by = "row_id") |>
    dplyr::left_join(trait_pick, by = "row_id") |>
    dplyr::filter(!is.na(raster_path), !is.na(trait_key)) |>
    dplyr::left_join(core$traits, by = c("taxon_class", "trait_key")) |>
    dplyr::select(-row_id)

  exclusions_after_sdm <- report_name_resolution(
    core$summary_df,
    raster_pick,
    species_inputs,
    raster_diagnostics,
    trait_diagnostics,
    verbose = verbose
  )
  attr(species_inputs, "stage4_exclusions_after_sdm") <- exclusions_after_sdm

  list(
    species_inputs = species_inputs,
    raster_pick = raster_pick,
    trait_pick = trait_pick,
    raster_diagnostics = raster_diagnostics,
    trait_diagnostics = trait_diagnostics
  )
}
attach_bird_generation_lengths <- function(species_inputs, bird_generation_lengths,
                                           synonyms, input_synonyms,
                                           verbose = FALSE) {
  out <- species_inputs
  out$GenLength <- NA_real_
  out$bird_genlength_scientific <- NA_character_
  out$genlength_match_source <- NA_character_
  exclusions_after_sdm <- attr(out, "stage4_exclusions_after_sdm")

  bird_rows <- which(out$taxon_class == "Aves")
  if (!length(bird_rows)) {
    report_stage4_final_resolution(
      out,
      exclusions_after_sdm,
      tibble::tibble(sdm_method = character(), genlength_match_source = character(), n = integer()),
      verbose = verbose
    )
    attr(out, "stage4_exclusions_after_sdm") <- NULL
    return(out)
  }

  query_names <- stringr::str_squish(as.character(out$scientificName[bird_rows]))
  if (any(is.na(query_names) | !nzchar(query_names))) {
    bad <- unique(out$scientificName[bird_rows][is.na(query_names) | !nzchar(query_names)])
    species_abort(
      "Retained bird rows have missing scientificName values, so Bird et al. GenLength cannot be matched:\n",
      paste(bad, collapse = "\n")
    )
  }

  unique_queries <- unique(query_names)
  match_result <- match_bird_generation_lengths(
    unique_queries,
    bird_generation_lengths = bird_generation_lengths,
    raw_synonyms = synonyms,
    raw_input_synonyms = input_synonyms,
    query_label = "Stage 4 retained bird species",
    synonym_mode = "attached",
    strict_ambiguity = FALSE
  )
  report_bird_generation_matches(
    match_result,
    "Stage 4 retained bird species",
    verbose = verbose
  )

  lookup <- match_result$matches |>
    dplyr::select(query_key, bird_scientific, GenLength, match_source)
  idx <- match(normalize_scientific_name(query_names), lookup$query_key)
  matched <- !is.na(idx)

  out$GenLength[bird_rows[matched]] <- lookup$GenLength[idx[matched]]
  out$bird_genlength_scientific[bird_rows[matched]] <- lookup$bird_scientific[idx[matched]]
  out$genlength_match_source[bird_rows[matched]] <- lookup$match_source[idx[matched]]

  bad_rows <- bird_rows[
    !matched |
      !is.finite(out$GenLength[bird_rows]) |
      out$GenLength[bird_rows] <= 0
  ]
  genlength_exclusions <- empty_sdm_matched_exclusion_table()
  if (length(bad_rows) > 0) {
    genlength_exclusions <- genlength_exclusions_after_sdm(out[bad_rows, , drop = FALSE], match_result)
    out <- out[-bad_rows, , drop = FALSE]
  }
  exclusions_after_sdm <- dplyr::bind_rows(exclusions_after_sdm, genlength_exclusions)
  attr(out, "stage4_exclusions_after_sdm") <- exclusions_after_sdm

  genlength_counts <- out |>
    dplyr::filter(taxon_class == "Aves") |>
    dplyr::count(sdm_method, genlength_match_source, name = "n") |>
    dplyr::arrange(sdm_method, genlength_match_source)
  report_stage4_final_resolution(out, exclusions_after_sdm, genlength_counts, verbose = verbose)
  attr(out, "stage4_exclusions_after_sdm") <- NULL

  out
}

