# Shared Stage 5-7 patch artifact and CSR contracts.
#
# This module contains no raster-processing workflow. It validates the stable
# lookup/connectivity schemas and provides rollback-capable file transactions.
# It deliberately depends only on the general project utilities: callers do
# not need to initialize Stage 5 paths, scientific settings, or a raster
# backend merely to read and validate an existing patch artifact.

patch_abort <- function(...) stop(paste0(...), call. = FALSE)

patch_lookup_required_cols <- function() {
  c("scientificName", "patch_id", "pu_id", "patch_area_km2")
}

patch_squish <- function(x) {
  gsub("\\s+", " ", trimws(as.character(x)))
}

validate_patch_positive_integer_column <- function(x, column, label) {
  values <- suppressWarnings(as.numeric(x[[column]]))
  invalid <- !is.finite(values) | values <= 0 | values != floor(values)
  if (any(invalid)) {
    patch_abort(
      label, " contains invalid ", column,
      " values; all ", column, " values must be positive integers."
    )
  }

  invisible(values)
}

validate_patch_positive_integer_values <- function(values, label, field) {
  values_num <- suppressWarnings(as.numeric(values))
  invalid <- !is.finite(values_num) | values_num <= 0 | values_num != floor(values_num)
  if (any(invalid)) {
    patch_abort(
      label, " contains invalid ", field,
      " values; all ", field, " values must be positive integers."
    )
  }

  invisible(values_num)
}

validate_patch_nonnegative_integer_values <- function(values, label, field) {
  values_num <- suppressWarnings(as.numeric(values))
  invalid <- !is.finite(values_num) | values_num < 0 | values_num != floor(values_num)
  if (any(invalid)) {
    patch_abort(
      label, " contains invalid ", field,
      " values; all ", field, " values must be non-negative integers."
    )
  }

  invisible(values_num)
}

set_patch_connectivity_attrs <- function(x) {
  attr(x, "csr_version") <- 2L
  attr(x, "row_ptr_base") <- "0-based"
  attr(x, "col_idx_space") <- "pu_local_index"
  attr(x, "col_idx_base") <- "0-based"
  x
}

validate_patch_connectivity_attrs <- function(x, label = "all_connectivity.rds") {
  expected <- list(
    csr_version = 2L,
    row_ptr_base = "0-based",
    col_idx_space = "pu_local_index",
    col_idx_base = "0-based"
  )

  for (nm in names(expected)) {
    found <- attr(x, nm, exact = TRUE)
    if (!identical(found, expected[[nm]])) {
      patch_abort(
        label, " has missing or invalid ", nm, " attribute. Expected ",
        expected[[nm]], "; found ",
        if (is.null(found)) "missing" else as.character(found),
        "."
      )
    }
  }

  invisible(TRUE)
}

validate_patch_lookup_object <- function(x, label = "all_patch_lookup.rds") {
  expected_columns <- patch_lookup_required_cols()
  if (!identical(names(x), expected_columns)) {
    patch_abort(
      label, " must contain exactly these ordered columns: ",
      paste(expected_columns, collapse = ", "), "."
    )
  }

  species <- patch_squish(x$scientificName)
  if (any(!nzchar(species))) {
    patch_abort(label, " contains blank scientificName values.")
  }
  if (!identical(as.character(x$scientificName), species)) {
    patch_abort(label, " contains scientificName values with non-normalized whitespace.")
  }
  patch_id <- validate_patch_positive_integer_column(x, "patch_id", label)
  pu_id <- validate_patch_positive_integer_column(x, "pu_id", label)
  if (anyDuplicated(paste0(species, "|", as.integer(patch_id)))) {
    patch_abort(label, " contains duplicate species-level patch_id values.")
  }
  if (any(!is.finite(as.numeric(x$patch_area_km2)) | as.numeric(x$patch_area_km2) <= 0)) {
    patch_abort(label, " contains missing, non-finite, or non-positive patch_area_km2 values.")
  }

  species_rows <- split(seq_along(species), species)
  for (species_name in names(species_rows)) {
    idx <- species_rows[[species_name]]
    found_patch_ids <- sort(unique(as.integer(patch_id[idx])))
    found_pu_ids <- sort(unique(as.integer(pu_id[idx])))
    if (!identical(found_patch_ids, seq_along(found_patch_ids))) {
      patch_abort(label, " patch_id values must equal 1..n within species: ", species_name, ".")
    }
    if (!identical(found_pu_ids, seq_along(found_pu_ids))) {
      patch_abort(label, " pu_id values must equal 1..n within species: ", species_name, ".")
    }
  }

  invisible(TRUE)
}

validate_patch_connectivity_object <- function(x, label = "all_connectivity.rds") {
  if (!is.list(x)) patch_abort(label, " must be a list.")
  validate_patch_connectivity_attrs(x, label)

  entry_names <- names(x)
  if (is.null(entry_names) || any(!nzchar(entry_names))) {
    patch_abort(label, " entries must be named as species|pu_id.")
  }
  dup_keys <- unique(entry_names[duplicated(entry_names)])
  if (length(dup_keys)) {
    patch_abort(label, " contains duplicate connectivity key(s): ", paste(dup_keys, collapse = ", "))
  }

  required_entry_cols <- c("species", "pu_id", "patch_ids", "row_ptr", "col_idx")
  for (key in entry_names) {
    entry <- x[[key]]
    if (!is.list(entry)) patch_abort(label, " entry is not a list: ", key)
    missing <- setdiff(required_entry_cols, names(entry))
    if (length(missing)) {
      patch_abort(label, " entry ", key, " is missing field(s): ", paste(missing, collapse = ", "))
    }

    species <- patch_squish(entry$species)
    if (!identical(length(species), 1L) || !nzchar(species)) {
      patch_abort(label, " entry ", key, " has a missing or blank species value.")
    }

    pu_id <- validate_patch_positive_integer_values(entry$pu_id, label, paste0(key, "$pu_id"))
    if (!identical(length(pu_id), 1L)) {
      patch_abort(label, " entry ", key, " must contain exactly one pu_id value.")
    }

    patch_ids <- validate_patch_positive_integer_values(entry$patch_ids, label, paste0(key, "$patch_ids"))
    if (!length(patch_ids)) {
      patch_abort(label, " entry ", key, " contains no patch_ids.")
    }
    if (anyDuplicated(patch_ids) || is.unsorted(patch_ids, strictly = TRUE)) {
      patch_abort(label, " entry ", key, " patch_ids must be strictly increasing and unique.")
    }

    row_ptr <- validate_patch_nonnegative_integer_values(entry$row_ptr, label, paste0(key, "$row_ptr"))
    col_idx <- validate_patch_nonnegative_integer_values(entry$col_idx, label, paste0(key, "$col_idx"))
    if (length(row_ptr) != length(patch_ids) + 1L) {
      patch_abort(label, " entry ", key, " has row_ptr length incompatible with patch_ids.")
    }
    if (!identical(row_ptr[[1L]], 0)) {
      patch_abort(label, " entry ", key, " row_ptr must start at 0.")
    }
    if (any(diff(row_ptr) < 0)) {
      patch_abort(label, " entry ", key, " row_ptr must be non-decreasing.")
    }
    if (!identical(row_ptr[[length(row_ptr)]], as.numeric(length(col_idx)))) {
      patch_abort(label, " entry ", key, " row_ptr terminus must equal the number of col_idx values.")
    }
    if (length(col_idx) && any(col_idx >= length(patch_ids))) {
      patch_abort(label, " entry ", key, " col_idx contains an out-of-range PU-local index.")
    }

    neighbors <- lapply(seq_along(patch_ids), function(i) {
      start <- row_ptr[[i]] + 1L
      end <- row_ptr[[i + 1L]]
      if (start > end) integer(0) else as.integer(col_idx[start:end])
    })
    for (i in seq_along(neighbors)) {
      nbrs <- neighbors[[i]]
      if (length(nbrs) && (anyDuplicated(nbrs) || is.unsorted(nbrs, strictly = TRUE))) {
        patch_abort(label, " entry ", key, " has unsorted or duplicate neighbor indices.")
      }
      if ((i - 1L) %in% nbrs) {
        patch_abort(label, " entry ", key, " contains a self-loop.")
      }
      if (length(nbrs)) {
        asymmetric <- vapply(nbrs, function(j) !(i - 1L) %in% neighbors[[j + 1L]], logical(1))
        if (any(asymmetric)) {
          patch_abort(label, " entry ", key, " contains asymmetric adjacency.")
        }
      }
    }
    if (length(patch_ids) == 1L && (length(col_idx) || !identical(as.integer(row_ptr), c(0L, 0L)))) {
      patch_abort(label, " entry ", key, " singleton PU must use an edge-free CSR representation.")
    }
    if (length(patch_ids) > 1L) {
      visited <- rep(FALSE, length(patch_ids))
      queue <- 1L
      visited[[1L]] <- TRUE
      while (length(queue)) {
        current <- queue[[1L]]
        queue <- queue[-1L]
        next_nodes <- neighbors[[current]] + 1L
        next_nodes <- next_nodes[!visited[next_nodes]]
        if (length(next_nodes)) {
          visited[next_nodes] <- TRUE
          queue <- c(queue, next_nodes)
        }
      }
      if (!all(visited)) patch_abort(label, " entry ", key, " is not connected.")
    }

    expected_key <- paste0(species, "|", as.integer(pu_id))
    if (!identical(key, expected_key)) {
      patch_abort(label, " entry key does not match species|pu_id. Found ", key, "; expected ", expected_key, ".")
    }
  }

  invisible(TRUE)
}

stage5_file_set_transaction <- function(staged_paths,
                                        target_paths,
                                        overwrite = FALSE,
                                        rename_file = file.rename) {
  project_file_set_transaction(
    staged_paths = staged_paths,
    target_paths = target_paths,
    overwrite = overwrite,
    rename_file = rename_file,
    defer_finalize = TRUE,
    label = "Stage 5"
  )
}

validate_patch_output_consistency <- function(patch_lookup,
                                              all_connectivity,
                                              label = "Stage 5 outputs") {
  validate_patch_lookup_object(patch_lookup, paste0(label, " patch lookup"))
  validate_patch_connectivity_object(all_connectivity, paste0(label, " connectivity"))

  lookup_species <- patch_squish(patch_lookup$scientificName)
  lookup_patch_id <- as.integer(suppressWarnings(as.numeric(patch_lookup$patch_id)))
  lookup_pu_id <- as.integer(suppressWarnings(as.numeric(patch_lookup$pu_id)))
  lookup_key <- paste0(lookup_species, "|", lookup_pu_id)

  duplicate_species_patch <- duplicated(paste0(lookup_species, "|", lookup_patch_id))
  if (any(duplicate_species_patch)) {
    patch_abort(
      label,
      " has duplicate species-level patch_id values in the patch lookup. First affected species:\n",
      paste(utils::head(unique(lookup_species[duplicate_species_patch]), 20L), collapse = "\n")
    )
  }

  lookup_keys <- sort(unique(lookup_key))
  connectivity_keys <- sort(names(all_connectivity))

  missing_connectivity <- setdiff(lookup_keys, connectivity_keys)
  if (length(missing_connectivity)) {
    patch_abort(
      label,
      " is incomplete. Patch lookup species|pu_id keys without connectivity entries:\n",
      paste(utils::head(missing_connectivity, 25L), collapse = "\n"),
      if (length(missing_connectivity) > 25L) "\n..." else ""
    )
  }

  extra_connectivity <- setdiff(connectivity_keys, lookup_keys)
  if (length(extra_connectivity)) {
    patch_abort(
      label,
      " is inconsistent. Connectivity species|pu_id keys without patch lookup rows:\n",
      paste(utils::head(extra_connectivity, 25L), collapse = "\n"),
      if (length(extra_connectivity) > 25L) "\n..." else ""
    )
  }

  lookup_patch_ids_by_key <- split(lookup_patch_id, lookup_key)
  bad_patch_sets <- character(0)
  for (key in connectivity_keys) {
    lookup_ids <- sort(unique(as.integer(lookup_patch_ids_by_key[[key]])))
    connectivity_ids <- sort(unique(as.integer(all_connectivity[[key]]$patch_ids)))

    if (!identical(lookup_ids, connectivity_ids)) {
      bad_patch_sets <- c(bad_patch_sets, key)
      if (length(bad_patch_sets) >= 25L) break
    }
  }

  if (length(bad_patch_sets)) {
    patch_abort(
      label,
      " has connectivity entries whose patch_ids do not match patch lookup rows. First affected keys:\n",
      paste(bad_patch_sets, collapse = "\n")
    )
  }

  invisible(TRUE)
}
