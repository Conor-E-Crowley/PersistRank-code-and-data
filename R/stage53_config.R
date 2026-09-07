# Stage 5.3 initial-persistence report configuration.
#
# Configuration selects taxa, one persistence curve, and canonical inputs and
# figure output. It reads no scientific data and writes nothing; the workflow
# owns the calculation and figure transaction.

stage53_config <- function(params, paths = project_paths(),
                           contract = canonical_analysis_contract()) {
  contract <- validate_analysis_contract(contract, "Stage 5.3 abundance contract")
  selection <- validate_taxa_selector(params$taxa, "params$taxa")
  list(
    taxa = selection$value,
    selected_mammals = selection$selected_mammals,
    selected_birds = selection$selected_birds,
    curve = validate_persistence_curve(params$curve, "params$curve"),
    contract = contract,
    paths = list(
      species_csv = file.path(paths$clean, "species_table.csv"),
      patch_lookup_rds = paths$patch_lookup %||%
        file.path(paths$clean, "all_patch_lookup.rds"),
      figure_dir = paths$spatial_figures %||% paths$figures,
      figure = file.path(
        paths$spatial_figures %||% paths$si_figures, "initial_persistence.png"
      )
    )
  )
}

stage53_application_config <- function(params, context) {
  stage53_config(
    list(taxa = context$taxa, curve = params$curve),
    context$paths, context$contract
  )
}

describe_stage53_config <- function(config) {
  log_msg("Stage 5.3 | curve=", config$curve, "| taxa=", config$taxa)
  invisible(config)
}

check_stage53_preflight <- function(config) {
  need_file(config$paths$species_csv, "Stage 4 species table")
  need_file(config$paths$patch_lookup_rds, "Stage 5 patch lookup")
  ensure_writable_dir(dirname(config$paths$figure), "Stage 5.3 SI figure directory")
  invisible(TRUE)
}
