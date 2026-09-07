# Standalone numbered-Rmd adaptation loaded by each numbered workflow.
#
# This definition-only module resolves Knit versus interactive Run All inputs
# and validates in-memory chunk ordering. It requires foundational and artifact
# helpers. Sourcing performs no render, package attachment, or filesystem write.

# Fail with an actionable chunk-order message when an interactive Rmd user runs
# a presentation chunk before its in-memory prerequisite.
require_chunk_object <- function(object_name, chunk_name,
                                 environment = parent.frame()) {
  assert(
    exists(object_name, envir = environment, inherits = FALSE),
    paste0("Run the '", chunk_name, "' chunk before this chunk.")
  )
  invisible(get(object_name, envir = environment, inherits = FALSE))
}

# TRUE only while an R Markdown document is being rendered. An explicit knitr
# flag is authoritative; interactive Run All always rereads the selected YAML.
project_knitr_active <- function() {
  knitr_state <- getOption("knitr.in.progress", NULL)
  if (!is.null(knitr_state)) return(isTRUE(knitr_state))
  if (interactive()) return(FALSE)
  if (!requireNamespace("knitr", quietly = TRUE)) return(FALSE)
  !is.null(knitr::opts_knit$get("rmarkdown.pandoc.to"))
}

# Return the Rmd currently owned by knitr, or NULL outside an active render.
project_knitr_input <- function() {
  if (!requireNamespace("knitr", quietly = TRUE)) return(NULL)
  input <- tryCatch(
    knitr::current_input(dir = TRUE),
    error = function(e) character(0)
  )
  if (length(input) != 1L || is.na(input) || !nzchar(input)) return(NULL)
  normalizePath(input, winslash = "/", mustWork = FALSE)
}

# Resolve unchanged YAML parameters for Knit and RStudio Run All. Render-time
# overrides are accepted only for the Rmd currently owned by knitr; otherwise
# the selected document is read so stale parameters cannot leak between Rmds.
project_rmd_params <- function(rmd_path, envir = parent.frame()) {
  assert(is.environment(envir), "envir must be an environment.")
  rmd_path <- validate_path_param(rmd_path, "rmd_path")
  need_file(rmd_path, "R Markdown driver")
  assert(!dir.exists(rmd_path), paste0("R Markdown driver is a directory: ", rmd_path))

  requested_rmd <- normalizePath(rmd_path, winslash = "/", mustWork = TRUE)
  rendering_this_rmd <- project_knitr_active() &&
    identical(project_knitr_input(), requested_rmd)
  has_render_params <- rendering_this_rmd &&
    exists("params", envir = envir, inherits = FALSE)

  if (has_render_params) {
    resolved <- get("params", envir = envir, inherits = FALSE)
  } else {
    assert(
      requireNamespace("rmarkdown", quietly = TRUE),
      "Package 'rmarkdown' is required to read Rmd parameters during Run All."
    )
    resolved <- rmarkdown::yaml_front_matter(rmd_path)$params
  }

  assert(is.list(resolved), "Rmd YAML must define a named params list.")
  param_names <- names(resolved)
  assert(
    length(resolved) > 0L &&
      !is.null(param_names) &&
      all(!is.na(param_names) & nzchar(param_names)) &&
      !anyDuplicated(param_names),
    "Rmd YAML must define a nonempty list of uniquely named parameters."
  )
  resolved
}
