# Hashed live storage for Stage 6 population-unit graphs.
#
# Bundles, checkpoints, and returned results continue using the established
# named-list contract. During prioritization, this small wrapper provides
# constant-time graph lookup and mutation without repeatedly copying a large
# named list.

is_priority_graph_store <- function(x) {
  is.environment(x) && inherits(x, "priority_graph_store")
}

new_priority_graph_store <- function(graphs) {
  if (is_priority_graph_store(graphs)) {
    return(graphs)
  }
  if (!is.list(graphs) || is.null(names(graphs))) {
    stop("PU graphs must be supplied as a named list.", call. = FALSE)
  }

  graph_keys <- names(graphs)
  if (anyNA(graph_keys) || any(!nzchar(graph_keys)) || anyDuplicated(graph_keys)) {
    stop("PU graph keys must be unique, nonblank strings.", call. = FALSE)
  }

  store <- new.env(parent = emptyenv())
  class(store) <- c("priority_graph_store", "environment")
  store$graphs <- new.env(hash = TRUE, parent = emptyenv(), size = max(29L, length(graphs)))
  store$positions <- new.env(hash = TRUE, parent = emptyenv(), size = max(29L, length(graphs)))
  store$order <- as.character(graph_keys)
  store$count <- as.integer(length(graphs))

  for (i in seq_along(graph_keys)) {
    graph_key <- graph_keys[[i]]
    assign(graph_key, graphs[[i]], envir = store$graphs)
    assign(graph_key, as.integer(i), envir = store$positions)
  }

  store
}

priority_graph_keys <- function(graphs) {
  if (is_priority_graph_store(graphs)) {
    return(graphs$order[!is.na(graphs$order)])
  }
  names(graphs)
}

priority_graph_has <- function(graphs, graph_key) {
  graph_key <- as.character(graph_key)
  if (is_priority_graph_store(graphs)) {
    return(vapply(
      graph_key,
      exists,
      logical(1L),
      envir = graphs$graphs,
      inherits = FALSE
    ))
  }
  !is.na(match(graph_key, names(graphs)))
}

priority_graph_get <- function(graphs, graph_key) {
  graph_key <- as.character(graph_key)
  if (length(graph_key) != 1L || is.na(graph_key) || !nzchar(graph_key)) {
    stop("graph_key must be one nonblank string.", call. = FALSE)
  }
  if (is_priority_graph_store(graphs)) {
    return(get0(graph_key, envir = graphs$graphs, inherits = FALSE, ifnotfound = NULL))
  }
  graphs[[graph_key]]
}

priority_graph_snapshot <- function(graphs, graph_keys, require_all = TRUE) {
  graph_keys <- as.character(graph_keys)
  if (!length(graph_keys)) {
    return(list())
  }
  if (anyNA(graph_keys) || any(!nzchar(graph_keys)) || anyDuplicated(graph_keys)) {
    stop("Graph snapshot keys must be unique, nonblank strings.", call. = FALSE)
  }

  present <- priority_graph_has(graphs, graph_keys)
  if (isTRUE(require_all) && any(!present)) {
    stop(
      "Missing PU graph key(s): ",
      paste(utils::head(graph_keys[!present], 5L), collapse = ", "),
      call. = FALSE
    )
  }

  selected_keys <- graph_keys[present]
  if (!length(selected_keys)) {
    return(list())
  }
  if (is_priority_graph_store(graphs)) {
    return(mget(selected_keys, envir = graphs$graphs, inherits = FALSE))
  }
  graphs[selected_keys]
}

normalize_priority_graph_actions <- function(actions) {
  if (!length(actions)) {
    return(list())
  }
  normalized <- vector("list", length(actions))
  for (i in seq_along(actions)) {
    action <- actions[[i]]
    if (!is.list(action) || is.null(action$key) || is.null(action$kind)) {
      stop("Every PU graph action requires kind and key fields.", call. = FALSE)
    }
    graph_key <- as.character(action$key)
    if (length(graph_key) != 1L || is.na(graph_key) || !nzchar(graph_key)) {
      stop("A PU graph action contains an invalid key.", call. = FALSE)
    }
    action_kind <- as.character(action$kind)
    remove_action <- action_kind %in% c("remove", "remove_original")
    set_action <- action_kind %in% c("set", "replace_original", "add_new")
    if (length(action_kind) != 1L || (!remove_action && !set_action)) {
      stop("A PU graph action contains an unsupported kind.", call. = FALSE)
    }
    if (set_action && is.null(action$graph)) {
      stop("A graph-setting action requires a graph object.", call. = FALSE)
    }
    normalized[[i]] <- list(
      kind = if (remove_action) "remove" else "set",
      key = graph_key,
      graph = if (set_action) action$graph else NULL
    )
  }
  normalized
}

priority_graph_apply_actions <- function(
  graphs,
  actions,
  assign_fn = base::assign,
  remove_fn = base::rm
) {
  actions <- normalize_priority_graph_actions(actions)
  if (!length(actions)) {
    return(graphs)
  }

  # Preserve the lower-level named-list interface used by compact tests and
  # standalone helpers. Production prioritization converts to the hashed store
  # once, before entering any iterative phase.
  if (!is_priority_graph_store(graphs)) {
    for (action in actions) {
      if (identical(action$kind, "remove")) {
        graphs[[action$key]] <- NULL
      } else {
        graphs[[action$key]] <- action$graph
      }
    }
    return(graphs)
  }

  affected_keys <- unique(vapply(actions, `[[`, character(1L), "key"))
  existed_before <- priority_graph_has(graphs, affected_keys)
  old_graphs <- priority_graph_snapshot(
    graphs,
    affected_keys[existed_before],
    require_all = TRUE
  )
  old_positions <- vapply(
    affected_keys[existed_before],
    function(key) get(key, envir = graphs$positions, inherits = FALSE),
    integer(1L)
  )
  old_order <- graphs$order
  old_count <- graphs$count

  apply_actions <- function() {
    updated_order <- old_order
    updated_count <- old_count

    for (action in actions) {
      graph_key <- action$key
      currently_present <- exists(graph_key, envir = graphs$graphs, inherits = FALSE)

      if (identical(action$kind, "remove")) {
        if (currently_present) {
          position <- get(graph_key, envir = graphs$positions, inherits = FALSE)
          updated_order[[position]] <- NA_character_
          remove_fn(list = graph_key, envir = graphs$graphs)
          remove_fn(list = graph_key, envir = graphs$positions)
          updated_count <- updated_count - 1L
        }
      } else if (currently_present) {
        # Ordinary list replacement retains the key's position.
        assign_fn(graph_key, action$graph, envir = graphs$graphs)
      } else {
        # Adding a new key, including a removed-and-readded key, appends it.
        updated_order <- c(updated_order, graph_key)
        assign_fn(graph_key, action$graph, envir = graphs$graphs)
        assign_fn(
          graph_key,
          as.integer(length(updated_order)),
          envir = graphs$positions
        )
        updated_count <- updated_count + 1L
      }
    }

    graphs$order <- updated_order
    graphs$count <- as.integer(updated_count)
    graphs
  }

  tryCatch(
    apply_actions(),
    error = function(error) {
      # Restore only keys touched by this batch. Unaffected graph bindings are
      # never copied or traversed.
      for (graph_key in affected_keys) {
        if (exists(graph_key, envir = graphs$graphs, inherits = FALSE)) {
          rm(list = graph_key, envir = graphs$graphs)
        }
        if (exists(graph_key, envir = graphs$positions, inherits = FALSE)) {
          rm(list = graph_key, envir = graphs$positions)
        }
      }
      if (length(old_graphs)) {
        for (graph_key in names(old_graphs)) {
          assign(graph_key, old_graphs[[graph_key]], envir = graphs$graphs)
          assign(graph_key, old_positions[[graph_key]], envir = graphs$positions)
        }
      }
      graphs$order <- old_order
      graphs$count <- old_count
      stop(error)
    }
  )
}

priority_graph_reposition_keys_to_end <- function(graphs, graph_keys) {
  graph_keys <- as.character(graph_keys)
  if (!length(graph_keys)) return(graphs)
  if (anyNA(graph_keys) || any(!nzchar(graph_keys)) || anyDuplicated(graph_keys) ||
      any(!priority_graph_has(graphs, graph_keys))) {
    stop("Graph keys repositioned after repair must be unique existing keys.", call. = FALSE)
  }

  if (!is_priority_graph_store(graphs)) {
    retained_keys <- names(graphs)[!names(graphs) %in% graph_keys]
    return(c(graphs[retained_keys], graphs[graph_keys]))
  }

  old_order <- graphs$order
  old_positions <- vapply(
    graph_keys,
    function(key) get(key, envir = graphs$positions, inherits = FALSE),
    integer(1L)
  )
  tryCatch(
    {
      updated_order <- old_order
      updated_order[old_positions] <- NA_character_
      updated_order <- c(updated_order, graph_keys)
      new_positions <- seq.int(
        length(updated_order) - length(graph_keys) + 1L,
        length(updated_order)
      )
      for (i in seq_along(graph_keys)) {
        assign(graph_keys[[i]], as.integer(new_positions[[i]]), envir = graphs$positions)
      }
      graphs$order <- updated_order
      graphs
    },
    error = function(error) {
      graphs$order <- old_order
      for (i in seq_along(graph_keys)) {
        assign(graph_keys[[i]], old_positions[[i]], envir = graphs$positions)
      }
      stop(error)
    }
  )
}

priority_graph_as_list <- function(graphs) {
  if (!is_priority_graph_store(graphs)) {
    return(graphs)
  }
  graph_keys <- priority_graph_keys(graphs)
  if (!length(graph_keys)) {
    return(stats::setNames(list(), character()))
  }
  mget(graph_keys, envir = graphs$graphs, inherits = FALSE)
}
