# Optional GRASS/fasterRaster runtime initialization for Stages 5 and 5.1.
#
# Stage 5 workflows load these definitions explicitly; Stage 8 inspection does
# not. Sourcing performs no discovery, package load, environment mutation, or
# initialization. Those effects occur only when the active workflow calls the
# documented resolver/initializer; terra fallback behavior is unchanged.

grass_config_dir <- function() {
  grass_bin <- Sys.which(c("grass", "grass84", "grass83", "grass82"))
  grass_bin <- grass_bin[nzchar(grass_bin)]
  if (!length(grass_bin)) return(NA_character_)

  for (bin in grass_bin) {
    out <- tryCatch(
      system2(bin, c("--config", "path"), stdout = TRUE, stderr = FALSE),
      error = function(e) character(0)
    )
    out <- optional_path(out[1L], "GRASS --config path output")
    if (!is.na(out) && dir.exists(out)) return(out)
  }

  NA_character_
}

resolve_grass_dir <- function(grass_dir = NULL, label = "grass_dir") {
  explicit <- optional_path(grass_dir, label)
  if (!is.na(explicit)) {
    if (!dir.exists(explicit)) {
      stop(
        paste0(label, " does not point to an existing GRASS GIS installation directory.\nPath: ", explicit),
        call. = FALSE
      )
    }
    return(normalizePath(explicit, winslash = "/", mustWork = TRUE))
  }

  env_candidates <- vapply(
    c("FASTER_RASTER_GRASS_DIR", "GRASS_DIR", "GISBASE"),
    function(x) optional_path(Sys.getenv(x, unset = NA_character_), paste0("environment variable ", x)),
    character(1)
  )
  path_candidate <- grass_config_dir()

  conda_prefix <- optional_path(
    Sys.getenv("CONDA_PREFIX", unset = NA_character_),
    "environment variable CONDA_PREFIX"
  )
  conda_candidates <- character(0)
  if (!is.na(conda_prefix)) {
    conda_candidates <- c(
      Sys.glob(file.path(conda_prefix, "grass*")),
      Sys.glob(file.path(conda_prefix, "lib", "grass*"))
    )
  }

  os_candidates <- c(
    Sys.glob("C:/Program Files/GRASS GIS *"),
    Sys.glob("/usr/local/grass*"),
    Sys.glob("/usr/lib/grass*"),
    Sys.glob("/opt/grass*")
  )

  candidates <- unique(c(explicit, env_candidates, path_candidate, conda_candidates, os_candidates))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  existing <- candidates[dir.exists(candidates)]

  if (!length(existing)) {
    stop(
      paste(
        "Could not find a GRASS GIS installation directory for fasterRaster.",
        "Set params$grass_dir, FASTER_RASTER_GRASS_DIR, GRASS_DIR, or GISBASE.",
        "For conda, activate the environment before rendering; if auto-detection fails,",
        "set FASTER_RASTER_GRASS_DIR to the GRASS directory, often $CONDA_PREFIX/lib/grass84.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  normalizePath(existing[1L], winslash = "/", mustWork = TRUE)
}

init_faster_raster <- function(grass_dir = NULL, ...) {
  if (!requireNamespace("fasterRaster", quietly = TRUE)) {
    warning(
      "fasterRaster is not installed; raster clumping will use terra::patches() fallback.",
      call. = FALSE
    )
    return(invisible(NA_character_))
  }

  resolved_grass_dir <- tryCatch(
    resolve_grass_dir(grass_dir),
    error = function(e) {
      warning(
        paste0(
          "Could not initialize fasterRaster because GRASS GIS was not found; ",
          "raster clumping will use terra::patches() fallback. ",
          "Original error: ", conditionMessage(e)
        ),
        call. = FALSE
      )
      NA_character_
    }
  )

  if (is.na(resolved_grass_dir)) {
    return(invisible(NA_character_))
  }

  configure_grass_path(resolved_grass_dir)
  fasterRaster::faster(grassDir = resolved_grass_dir, ...)
  invisible(resolved_grass_dir)
}

configure_grass_path <- function(grass_dir) {
  grass_dir <- normalizePath(grass_dir, winslash = "/", mustWork = TRUE)
  grass_path_dirs <- file.path(grass_dir, c("bin", "extrabin", "scripts", "lib"))
  grass_path_dirs <- normalizePath(
    grass_path_dirs[dir.exists(grass_path_dirs)],
    winslash = "/",
    mustWork = TRUE
  )

  current_path <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  current_path_norm <- normalizePath(
    current_path[nzchar(current_path) & dir.exists(current_path)],
    winslash = "/",
    mustWork = FALSE
  )
  add_path_dirs <- grass_path_dirs[!grass_path_dirs %in% current_path_norm]

  if (length(add_path_dirs)) {
    Sys.setenv(PATH = paste(c(add_path_dirs, Sys.getenv("PATH")), collapse = .Platform$path.sep))
  }
  Sys.setenv(GISBASE = grass_dir)

  invisible(grass_dir)
}
