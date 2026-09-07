# Dynamic demographic-model contracts shared by Stages 1–4.
#
# EltonTraits supplies five mutually exclusive bird Diet-5Cat values. Stage 1
# may fit a separate environmental-variation intercept for any proper subset;
# every unselected category shares the reference branch named "Other". The
# returned specification is serialized with Stage 2/3 artifacts so downstream
# stages never need to guess how raw diets map to fitted branches.


bird_diet5_categories <- function() {
  c("FruiNect", "Invertebrate", "Omnivore", "PlantSeed", "VertFishScav")
}

normalize_bird_separate_intercepts <- function(
  x,
  label = "bird_sigma_separate_intercepts"
) {
  values <- unlist(x, recursive = TRUE, use.names = FALSE)
  if (!length(values)) return(character())
  assert(is.character(values), paste0(label, " must be a character vector."))
  values <- trimws(values)
  assert(
    all(!is.na(values) & nzchar(values)),
    paste0(label, " must not contain blank or missing values.")
  )
  assert(!anyDuplicated(values), paste0(label, " must not contain duplicates."))
  invalid <- setdiff(values, bird_diet5_categories())
  assert(
    !length(invalid),
    paste0(
      label, " contains unsupported Diet-5Cat value(s): ",
      paste(invalid, collapse = ", "), ". Valid values are: ",
      paste(bird_diet5_categories(), collapse = ", "), "."
    )
  )
  ordered <- bird_diet5_categories()[bird_diet5_categories() %in% values]
  assert(
    length(ordered) < length(bird_diet5_categories()),
    paste0(label, " must leave at least one Diet-5Cat category in Other.")
  )
  ordered
}

bird_diet_coefficient_name <- function(group) {
  paste0("beta_diet_", group)
}

bird_demographic_model_spec <- function(
  separate_intercepts = "VertFishScav",
  observed_diet5 = NULL,
  require_observed = FALSE
) {
  separate <- normalize_bird_separate_intercepts(separate_intercepts)
  reference <- setdiff(bird_diet5_categories(), separate)

  if (!is.null(observed_diet5)) {
    observed <- validate_bird_diet5(observed_diet5, "observed bird Diet-5Cat")
    counts <- table(factor(observed, levels = bird_diet5_categories()))
    if (isTRUE(require_observed)) {
      absent <- separate[counts[separate] == 0L]
      assert(
        !length(absent),
        paste0(
          "Selected separate-intercept Diet-5Cat value(s) are absent from the ",
          "calibration data: ", paste(absent, collapse = ", "), "."
        )
      )
      assert(
        sum(counts[reference]) > 0L,
        "No calibration observations remain in the Other reference branch."
      )
    }
  } else {
    counts <- NULL
  }

  raw_to_model <- stats::setNames(
    ifelse(bird_diet5_categories() %in% separate, bird_diet5_categories(), "Other"),
    bird_diet5_categories()
  )
  model_groups <- c("Other", separate)

  list(
    schema_version = 1L,
    raw_categories = bird_diet5_categories(),
    separate_intercepts = separate,
    reference_categories = reference,
    model_groups = model_groups,
    coefficient_names = bird_diet_coefficient_name(separate),
    raw_to_model = raw_to_model,
    raw_counts = counts
  )
}

validate_bird_diet5 <- function(x, label = "bird Diet-5Cat") {
  values <- trimws(as.character(x))
  assert(
    length(values) > 0L &&
      all(!is.na(values) & nzchar(values)) &&
      all(values %in% bird_diet5_categories()),
    paste0(
      label, " must contain only: ",
      paste(bird_diet5_categories(), collapse = ", "), "."
    )
  )
  values
}

validate_bird_model_group <- function(x, spec, label = "bird model group") {
  assert(is.list(spec) && length(spec$model_groups), "A valid bird model specification is required.")
  values <- trimws(as.character(x))
  assert(
    length(values) > 0L &&
      all(!is.na(values) & nzchar(values)) &&
      all(values %in% spec$model_groups),
    paste0(label, " must contain only: ", paste(spec$model_groups, collapse = ", "), ".")
  )
  values
}

map_bird_diet_to_model_group <- function(diet5, spec) {
  diet5 <- validate_bird_diet5(diet5)
  unname(spec$raw_to_model[diet5])
}

bird_model_spec_from_posterior <- function(columns) {
  prefix <- "beta_diet_"
  diet_columns <- columns[startsWith(columns, prefix)]
  separate <- sub(paste0("^", prefix), "", diet_columns)
  spec <- bird_demographic_model_spec(separate)
  assert(
    identical(diet_columns, spec$coefficient_names),
    paste0(
      "Bird environmental-variation posterior diet coefficients must be in ",
      "canonical Diet-5Cat order: ",
      paste(spec$coefficient_names, collapse = ", "), "."
    )
  )
  spec
}

bird_model_design_matrix <- function(generation_length, diet5, spec) {
  generation_length <- suppressWarnings(as.numeric(generation_length))
  assert(
    length(generation_length) == length(diet5) &&
      all(is.finite(generation_length) & generation_length > 0),
    "Bird environmental-variation generation lengths must be positive and finite."
  )
  diet5 <- validate_bird_diet5(diet5)
  indicators <- vapply(
    spec$separate_intercepts,
    function(group) as.numeric(diet5 == group),
    numeric(length(diet5))
  )
  if (!length(spec$separate_intercepts)) {
    indicators <- matrix(numeric(), nrow = length(diet5), ncol = 0L)
  } else if (is.null(dim(indicators))) {
    indicators <- matrix(indicators, ncol = 1L)
  }
  colnames(indicators) <- sub("^beta_", "", spec$coefficient_names)
  cbind(logGenLength = log10(generation_length), indicators)
}

validate_bird_model_spec <- function(spec) {
  expected <- bird_demographic_model_spec(spec$separate_intercepts)
  fields <- c(
    "schema_version", "raw_categories", "separate_intercepts",
    "reference_categories", "model_groups", "coefficient_names", "raw_to_model"
  )
  assert(all(fields %in% names(spec)), "Bird demographic model specification is incomplete.")
  for (field in fields) {
    assert(
      identical(spec[[field]], expected[[field]]),
      paste0("Bird demographic model specification has invalid ", field, ".")
    )
  }
  invisible(TRUE)
}

bird_model_signature <- function(spec) {
  validate_bird_model_spec(spec)
  paste(names(spec$raw_to_model), spec$raw_to_model, sep = ":", collapse = ",")
}

bird_model_spec_from_signature <- function(signature) {
  signature <- validate_scalar_string(signature, "bird model signature")
  pieces <- strsplit(signature, ",", fixed = TRUE)[[1L]]
  pairs <- strsplit(pieces, ":", fixed = TRUE)
  assert(
    all(lengths(pairs) == 2L),
    "Bird model signature is malformed."
  )
  raw <- vapply(pairs, `[[`, character(1), 1L)
  modeled <- vapply(pairs, `[[`, character(1), 2L)
  assert(
    identical(raw, bird_diet5_categories()),
    "Bird model signature has invalid Diet-5Cat order."
  )
  separate <- raw[modeled == raw]
  spec <- bird_demographic_model_spec(separate)
  assert(
    identical(bird_model_signature(spec), signature),
    "Bird model signature has an invalid raw-to-model mapping."
  )
  spec
}
