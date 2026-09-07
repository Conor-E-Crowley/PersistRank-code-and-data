# Foundational project helpers loaded by every numbered workflow.
#
# This definition-only module owns assertions, scalar validation, optional path
# handling, package checks, stable species identifiers, and console logging.
# Filesystem artifacts, Rmd adaptation, runtime logs, transactions, and spatial
# runtime initialization are owned by focused modules loaded after this file.
# Sourcing has no side effects and loads no expensive dependency.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}

validate_parameter_names <- function(params, allowed, label) {
  assert(is.list(params), paste0(label, " parameters must be supplied as a list."))
  unsupported <- setdiff(names(params), allowed)
  assert(!length(unsupported), paste0(
    label, " does not accept: ", paste(unsupported, collapse = ", "), "."
  ))
  invisible(params)
}

validate_scalar_string <- function(x, label, allow_empty = FALSE) {
  if (is.character(x)) x <- trimws(x)
  assert(
    is.character(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      (isTRUE(allow_empty) || nzchar(x)),
    paste0(label, " must be a single", if (isTRUE(allow_empty)) "" else " non-empty", " string.")
  )
  x
}

validate_scalar_choice <- function(x, choices, label) {
  value <- validate_scalar_string(x, label)
  assert(
    value %in% choices,
    paste0(label, " must be one of: ", paste(choices, collapse = ", "), ".")
  )
  value
}

validate_scalar_logical <- function(x, label) {
  assert(
    is.logical(x) && length(x) == 1L && !is.na(x),
    paste0(label, " must be TRUE/FALSE.")
  )
  isTRUE(x)
}

validate_scalar_number <- function(
  x,
  label,
  minimum = -Inf,
  maximum = Inf,
  minimum_open = FALSE,
  maximum_open = FALSE
) {
  assert(
    length(x) == 1L && !is.logical(x),
    paste0(label, " must be one numeric value.")
  )
  value <- suppressWarnings(as.numeric(x))
  lower_ok <- if (minimum_open) value > minimum else value >= minimum
  upper_ok <- if (maximum_open) value < maximum else value <= maximum
  assert(
    is.finite(value) && lower_ok && upper_ok,
    paste0(
      label, " must be one finite number in ",
      if (minimum_open) "(" else "[",
      minimum, ", ", maximum,
      if (maximum_open) ")" else "]."
    )
  )
  value
}

validate_scalar_integer <- function(
  x,
  label,
  minimum = 1L,
  maximum = .Machine$integer.max
) {
  value <- validate_scalar_number(
    x,
    label,
    minimum = minimum,
    maximum = maximum
  )
  assert(value == floor(value), paste0(label, " must be an integer."))
  as.integer(value)
}

validate_taxa_selector <- function(x, label = "taxa") {
  value <- validate_scalar_choice(x, c("mammals", "birds", "both"), label)
  list(
    value = value,
    selected_mammals = value %in% c("mammals", "both"),
    selected_birds = value %in% c("birds", "both")
  )
}

validate_path_param <- function(path, label = "path") {
  path.expand(validate_scalar_string(path, label))
}

need_cols <- function(x, cols, label) {
  missing_cols <- setdiff(cols, names(x))
  assert(
    length(missing_cols) == 0L,
    paste0(label, " is missing required column(s): ", paste(missing_cols, collapse = ", "))
  )
}

optional_path <- function(path, label = "path") {
  if (is.null(path) || !length(path)) return(NA_character_)
  assert(
    is.character(path) && length(path) == 1L,
    paste0(label, " must be NULL, empty, or a single path string.")
  )
  path <- trimws(path)
  if (is.na(path) || !nzchar(path) || toupper(path) %in% c("NA", "NULL")) return(NA_character_)
  path.expand(path)
}

check_packages <- function(packages) {
  packages <- unique(as.character(packages))
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  assert(
    length(missing) == 0L,
    paste0(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      "\nInstall them before running this script."
    )
  )
  invisible(packages)
}

load_packages <- function(packages) {
  packages <- check_packages(packages)
  suppressPackageStartupMessages(
    invisible(lapply(packages, library, character.only = TRUE))
  )
  invisible(packages)
}

species_id <- function(scientific_name) {
  x <- trimws(as.character(scientific_name))
  x <- gsub("\\s+", "_", x)
  x <- gsub("[/\\\\:<>\"|?*]+", "_", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

patch_filename_from_scientific <- function(scientific_name) {
  paste0(species_id(scientific_name), ".tif")
}

log_msg <- function(..., flush = FALSE) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " ")))
  if (isTRUE(flush)) flush.console()
  invisible(NULL)
}

log_space <- function(x_min, x_max, n) {
  assert(
    all(is.finite(c(x_min, x_max, n))) && x_min > 0 && x_max > x_min && n >= 2,
    "log_space() requires 0 < x_min < x_max and n >= 2."
  )
  10^seq(log10(x_min), log10(x_max), length.out = as.integer(n))
}
