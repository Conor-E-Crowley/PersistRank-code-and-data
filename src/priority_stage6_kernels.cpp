#include <Rcpp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Stage 6 compiled-kernel contract:
// - R-facing cell, patch, PU, and row IDs are one-based; C++ vector offsets are
//   zero-based. CSR row pointers are zero-based, while neighbor-row encoding is
//   validated according to each named graph boundary.
// - Compact patch indexes are sorted by patch ID with unique cells. Set-like
//   inputs are either required unique or canonicalized explicitly by the kernel.
// - Output ordering and ties are deterministic. Large outputs are sized or
//   reserved before filling, and no kernel retains state after returning to R.
// - Interrupt checks sit at bounded outer-loop boundaries, outside hot inner
//   loops where checking would alter the performance contract.

// [[Rcpp::plugins(cpp11)]]

namespace {

bool contains_sorted(const std::vector<int>& values, int value) {
  return std::binary_search(values.begin(), values.end(), value);
}

std::vector<int> sorted_unique(const Rcpp::IntegerVector& values,
                               const char* label) {
  std::vector<int> output;
  output.reserve(values.size());
  for (int value : values) {
    if (value == NA_INTEGER || value < 1) {
      Rcpp::stop("%s must contain positive integers.", label);
    }
    output.push_back(value);
  }
  std::sort(output.begin(), output.end());
  output.erase(std::unique(output.begin(), output.end()), output.end());
  return output;
}

}  // namespace


// ---- Runtime contract and sparse cell-species ownership --------------------

// Version tag consumed by R-side compatibility checks.
// [[Rcpp::export]]
Rcpp::String stage6_runtime_kernel_contract() {
  return "stage6_runtime_kernels_v11";
}


// Build one cell-major membership index. Raw species IDs are zero-based bytes;
// integer species IDs remain one-based for R consumers.
// [[Rcpp::export]]
Rcpp::List stage6_build_cell_species_index_cpp(
    const Rcpp::List& cells_by_species,
    const int n_cells,
    const bool use_raw_species_ids = true) {
  if (n_cells < 1) {
    Rcpp::stop("Cell-species index requires a positive cell count.");
  }
  const R_xlen_t n_species = cells_by_species.size();
  if (use_raw_species_ids && n_species > 255) {
    Rcpp::stop("Raw cell-species IDs support at most 255 species.");
  }

  Rcpp::IntegerVector counts(n_cells);
  double total_entries_double = 0;
  for (R_xlen_t species_index = 0; species_index < n_species; ++species_index) {
    if ((species_index & 31) == 0) Rcpp::checkUserInterrupt();
    if (Rf_isNull(cells_by_species[species_index])) continue;
    const Rcpp::IntegerVector cells = cells_by_species[species_index];
    std::unordered_set<int> seen;
    seen.reserve(static_cast<std::size_t>(cells.size()));
    for (R_xlen_t i = 0; i < cells.size(); ++i) {
      const int cell = cells[i];
      if (cell == NA_INTEGER || cell < 1 || cell > n_cells) {
        Rcpp::stop("Cell-species index received an invalid raster cell.");
      }
      if (!seen.insert(cell).second) {
        Rcpp::stop("One species contains a duplicated raster cell.");
      }
      if (counts[cell - 1] == std::numeric_limits<int>::max()) {
        Rcpp::stop("Cell-species membership count exceeds integer limits.");
      }
      counts[cell - 1] += 1;
      total_entries_double += 1;
    }
  }
  if (total_entries_double > static_cast<double>(std::numeric_limits<int>::max())) {
    Rcpp::stop("Cell-species index exceeds the supported integer offset range.");
  }
  const int total_entries = static_cast<int>(total_entries_double);
  Rcpp::IntegerVector offsets(n_cells + 1);
  for (int cell = 0; cell < n_cells; ++cell) {
    offsets[cell + 1] = offsets[cell] + counts[cell];
  }
  Rcpp::IntegerVector cursor = Rcpp::clone(offsets);

  SEXP species_ids;
  if (use_raw_species_ids) {
    Rcpp::RawVector output(total_entries);
    for (R_xlen_t species_index = 0; species_index < n_species; ++species_index) {
      if (Rf_isNull(cells_by_species[species_index])) continue;
      const Rcpp::IntegerVector cells = cells_by_species[species_index];
      for (R_xlen_t i = 0; i < cells.size(); ++i) {
        const int cell_offset = cells[i] - 1;
        output[cursor[cell_offset]++] = static_cast<Rbyte>(species_index);
      }
    }
    species_ids = output;
  } else {
    Rcpp::IntegerVector output(total_entries);
    for (R_xlen_t species_index = 0; species_index < n_species; ++species_index) {
      if (Rf_isNull(cells_by_species[species_index])) continue;
      const Rcpp::IntegerVector cells = cells_by_species[species_index];
      for (R_xlen_t i = 0; i < cells.size(); ++i) {
        const int cell_offset = cells[i] - 1;
        output[cursor[cell_offset]++] = static_cast<int>(species_index + 1);
      }
    }
    species_ids = output;
  }

  return Rcpp::List::create(
    Rcpp::Named("offsets") = offsets,
    Rcpp::Named("species_ids") = species_ids,
    Rcpp::Named("membership_count_by_cell") = counts,
    Rcpp::Named("entry_count") = static_cast<double>(total_entries)
  );
}


// ---- Population-unit graph repair ------------------------------------------

// Rebuild affected PUs from their zero-based CSR graphs. Components and new PU
// IDs are emitted in deterministic graph, node, and component order.
// [[Rcpp::export]]
Rcpp::List stage6_rebuild_species_pus_cpp(
    const Rcpp::List& pu_graphs,
    const Rcpp::List& dead_patch_ids_by_pu,
    const Rcpp::List& patch_area_by_pu,
    const double pu_area_threshold,
    int next_available_pu_id) {
  if (!R_finite(pu_area_threshold) || pu_area_threshold < 0) {
    Rcpp::stop("PU-area threshold must be finite and non-negative.");
  }
  if (next_available_pu_id == NA_INTEGER || next_available_pu_id < 0) {
    Rcpp::stop("Next available PU ID must be a non-negative integer.");
  }
  if (pu_graphs.size() != dead_patch_ids_by_pu.size() ||
      pu_graphs.size() != patch_area_by_pu.size()) {
    Rcpp::stop("Compiled PU-repair inputs must have equal lengths.");
  }

  Rcpp::List updates(pu_graphs.size());
  double graph_nodes_examined = 0;
  double adjacency_entries_examined = 0;
  double dead_patch_nodes = 0;
  double surviving_components = 0;
  double split_pus = 0;
  double component_dropped_patches = 0;

  for (R_xlen_t graph_index = 0; graph_index < pu_graphs.size(); ++graph_index) {
    if ((graph_index & 31) == 0) Rcpp::checkUserInterrupt();
    const Rcpp::List graph = pu_graphs[graph_index];
    const int pu_id = Rcpp::as<int>(graph["pu_id"]);
    const Rcpp::IntegerVector id2patch = graph["id2patch"];
    const Rcpp::IntegerVector row_ptr = graph["row_ptr"];
    const Rcpp::IntegerVector col_idx = graph["col_idx"];
    const Rcpp::IntegerVector dead_input = dead_patch_ids_by_pu[graph_index];
    const Rcpp::NumericVector patch_area = patch_area_by_pu[graph_index];
    const int node_count = id2patch.size();

    if (pu_id == NA_INTEGER || pu_id < 1 ||
        row_ptr.size() != static_cast<R_xlen_t>(node_count) + 1 ||
        row_ptr[0] != 0 || row_ptr[node_count] != col_idx.size() ||
        patch_area.size() != node_count) {
      Rcpp::stop("Compiled PU repair received an invalid graph.");
    }
    int previous_offset = 0;
    for (R_xlen_t i = 0; i < row_ptr.size(); ++i) {
      const int offset = row_ptr[i];
      if (offset == NA_INTEGER || offset < previous_offset ||
          offset < 0 || offset > col_idx.size()) {
        Rcpp::stop("Compiled PU repair received invalid CSR offsets.");
      }
      previous_offset = offset;
    }
    std::unordered_set<int> seen_patch_ids;
    seen_patch_ids.reserve(static_cast<std::size_t>(node_count));
    for (int node = 0; node < node_count; ++node) {
      if (id2patch[node] == NA_INTEGER || id2patch[node] < 1 ||
          !seen_patch_ids.insert(id2patch[node]).second ||
          !R_finite(patch_area[node]) || patch_area[node] < 0) {
        Rcpp::stop("Compiled PU repair received invalid node data.");
      }
    }
    for (R_xlen_t i = 0; i < col_idx.size(); ++i) {
      if (col_idx[i] == NA_INTEGER || col_idx[i] < 1 ||
          col_idx[i] > node_count) {
        Rcpp::stop("Compiled PU repair received invalid CSR adjacency.");
      }
    }

    std::unordered_set<int> dead_ids;
    dead_ids.reserve(static_cast<std::size_t>(dead_input.size()));
    for (int patch_id : dead_input) {
      if (patch_id == NA_INTEGER || patch_id < 1) {
        Rcpp::stop("Compiled PU repair received an invalid dead patch ID.");
      }
      dead_ids.insert(patch_id);
    }
    std::vector<unsigned char> alive(static_cast<std::size_t>(node_count), 1);
    int alive_count = 0;
    for (int node = 0; node < node_count; ++node) {
      if (dead_ids.find(id2patch[node]) != dead_ids.end()) {
        alive[static_cast<std::size_t>(node)] = 0;
      } else {
        ++alive_count;
      }
    }

    graph_nodes_examined += node_count;
    adjacency_entries_examined += col_idx.size();
    dead_patch_nodes += node_count - alive_count;
    Rcpp::List surviving_graphs;
    std::vector<int> dropped_patch_ids;

    if (alive_count == 1) {
      int surviving_node = -1;
      for (int node = 0; node < node_count; ++node) {
        if (alive[static_cast<std::size_t>(node)]) {
          surviving_node = node;
          break;
        }
      }
      if (patch_area[surviving_node] > pu_area_threshold) {
        surviving_graphs = Rcpp::List::create(Rcpp::List::create(
          Rcpp::Named("pu_id") = pu_id,
          Rcpp::Named("id2patch") =
            Rcpp::IntegerVector::create(id2patch[surviving_node]),
          Rcpp::Named("row_ptr") = Rcpp::IntegerVector::create(0, 0),
          Rcpp::Named("col_idx") = Rcpp::IntegerVector(0)
        ));
      } else {
        dropped_patch_ids.push_back(id2patch[surviving_node]);
      }
    } else if (alive_count > 1) {
      std::vector<int> component_by_node(static_cast<std::size_t>(node_count), 0);
      std::vector<double> component_area;
      std::vector<std::vector<int> > component_nodes;
      int component_count = 0;

      for (int seed = 0; seed < node_count; ++seed) {
        if (!alive[static_cast<std::size_t>(seed)] ||
            component_by_node[static_cast<std::size_t>(seed)] != 0) continue;
        ++component_count;
        std::vector<int> queue;
        queue.reserve(static_cast<std::size_t>(alive_count));
        queue.push_back(seed);
        component_by_node[static_cast<std::size_t>(seed)] = component_count;
        std::size_t head = 0;
        while (head < queue.size()) {
          const int node = queue[head++];
          for (int offset = row_ptr[node]; offset < row_ptr[node + 1]; ++offset) {
            const int neighbor = col_idx[offset] - 1;
            if (alive[static_cast<std::size_t>(neighbor)] &&
                component_by_node[static_cast<std::size_t>(neighbor)] == 0) {
              component_by_node[static_cast<std::size_t>(neighbor)] =
                component_count;
              queue.push_back(neighbor);
            }
          }
        }
      }

      component_area.assign(static_cast<std::size_t>(component_count), 0.0);
      component_nodes.resize(static_cast<std::size_t>(component_count));
      for (int node = 0; node < node_count; ++node) {
        const int component = component_by_node[static_cast<std::size_t>(node)];
        if (component > 0) {
          component_area[static_cast<std::size_t>(component - 1)] +=
            patch_area[node];
          component_nodes[static_cast<std::size_t>(component - 1)].push_back(node);
        }
      }

      std::vector<int> surviving_labels;
      for (int component = 0; component < component_count; ++component) {
        if (component_area[static_cast<std::size_t>(component)] >
            pu_area_threshold) {
          surviving_labels.push_back(component);
        }
      }
      surviving_graphs = Rcpp::List(surviving_labels.size());

      for (std::size_t surviving_index = 0;
           surviving_index < surviving_labels.size(); ++surviving_index) {
        const int component = surviving_labels[surviving_index];
        const std::vector<int>& nodes =
          component_nodes[static_cast<std::size_t>(component)];
        const int output_pu_id = surviving_index == 0
          ? pu_id
          : ++next_available_pu_id;
        std::vector<int> old_to_new(static_cast<std::size_t>(node_count), -1);
        Rcpp::IntegerVector output_patch_ids(nodes.size());
        for (std::size_t i = 0; i < nodes.size(); ++i) {
          old_to_new[static_cast<std::size_t>(nodes[i])] =
            static_cast<int>(i);
          output_patch_ids[i] = id2patch[nodes[i]];
        }
        Rcpp::IntegerVector output_row_ptr(nodes.size() + 1);
        std::vector<int> output_col_idx;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
          const int old_node = nodes[i];
          for (int offset = row_ptr[old_node];
               offset < row_ptr[old_node + 1]; ++offset) {
            const int new_neighbor =
              old_to_new[static_cast<std::size_t>(col_idx[offset] - 1)];
            if (new_neighbor >= 0) output_col_idx.push_back(new_neighbor + 1);
          }
          output_row_ptr[i + 1] = output_col_idx.size();
        }
        surviving_graphs[surviving_index] = Rcpp::List::create(
          Rcpp::Named("pu_id") = output_pu_id,
          Rcpp::Named("id2patch") = output_patch_ids,
          Rcpp::Named("row_ptr") = output_row_ptr,
          Rcpp::Named("col_idx") = Rcpp::wrap(output_col_idx)
        );
      }

      for (int node = 0; node < node_count; ++node) {
        const int component = component_by_node[static_cast<std::size_t>(node)];
        if (component > 0 &&
            !(component_area[static_cast<std::size_t>(component - 1)] >
              pu_area_threshold)) {
          dropped_patch_ids.push_back(id2patch[node]);
        }
      }
    }

    surviving_components += surviving_graphs.size();
    split_pus += surviving_graphs.size() > 1 ? 1 : 0;
    component_dropped_patches += dropped_patch_ids.size();
    updates[graph_index] = Rcpp::List::create(
      Rcpp::Named("surviving_pu_graphs") = surviving_graphs,
      Rcpp::Named("dropped_patch_ids") = Rcpp::wrap(dropped_patch_ids)
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("updates") = updates,
    Rcpp::Named("next_available_pu_id") = next_available_pu_id,
    Rcpp::Named("profiling_counts") = Rcpp::NumericVector::create(
      Rcpp::Named("graph_nodes_examined") = graph_nodes_examined,
      Rcpp::Named("adjacency_entries_examined") = adjacency_entries_examined,
      Rcpp::Named("dead_patch_nodes") = dead_patch_nodes,
      Rcpp::Named("surviving_components") = surviving_components,
      Rcpp::Named("split_pus") = split_pus,
      Rcpp::Named("component_dropped_patches") =
        component_dropped_patches,
      Rcpp::Named("compiled_pu_workloads") =
        static_cast<double>(pu_graphs.size())
    )
  );
}


// ---- Distance-edge reconstruction and canonical patch lookup ---------------

// Collect unique undirected distance-recheck edges and sorted candidate patch
// IDs without materializing a dense patch-by-patch matrix.
// [[Rcpp::export]]
Rcpp::List stage6_extract_distance_recheck_edges_cpp(
    const Rcpp::List& affected_graphs,
    const Rcpp::IntegerVector& recheck_patch_ids_input) {
  const std::vector<int> recheck_patch_ids = sorted_unique(
    recheck_patch_ids_input,
    "Distance-recheck patch IDs"
  );
  std::unordered_set<int> recheck_patch_set(
    recheck_patch_ids.begin(), recheck_patch_ids.end()
  );

  std::vector<int> edge_u;
  std::vector<int> edge_v;
  std::unordered_set<std::uint64_t> seen_edges;
  std::vector<int> candidate_patch_ids = recheck_patch_ids;

  double graph_nodes_scanned = 0;
  double recheck_rows_total = 0;
  double csr_entries_visited = 0;
  double overlay_edges_scanned = 0;

  auto validate_patch_ids = [](const Rcpp::IntegerVector& patch_ids,
                               int graph_position,
                               const char* label) {
    std::unordered_set<int> seen;
    seen.reserve(static_cast<std::size_t>(patch_ids.size()));
    for (R_xlen_t i = 0; i < patch_ids.size(); ++i) {
      const int patch_id = patch_ids[i];
      if (patch_id == NA_INTEGER || patch_id < 1) {
        Rcpp::stop("%s contains an invalid patch ID in affected graph %d.",
                   label, graph_position);
      }
      if (!seen.insert(patch_id).second) {
        Rcpp::stop("%s contains duplicate patch IDs in affected graph %d.",
                   label, graph_position);
      }
    }
  };

  auto validate_csr = [](const Rcpp::IntegerVector& row_ptr,
                         const Rcpp::IntegerVector& col_idx,
                         int node_count,
                         int graph_position,
                         const char* label) {
    if (row_ptr.size() != static_cast<R_xlen_t>(node_count) + 1) {
      Rcpp::stop("%s row_ptr length is invalid in affected graph %d.",
                 label, graph_position);
    }
    if (row_ptr[0] != 0) {
      Rcpp::stop("%s row_ptr must begin at zero in affected graph %d.",
                 label, graph_position);
    }
    int previous = 0;
    for (R_xlen_t i = 0; i < row_ptr.size(); ++i) {
      const int offset = row_ptr[i];
      if (offset == NA_INTEGER || offset < previous || offset < 0) {
        Rcpp::stop("%s row_ptr is invalid in affected graph %d.",
                   label, graph_position);
      }
      previous = offset;
    }
    if (row_ptr[row_ptr.size() - 1] != col_idx.size()) {
      Rcpp::stop("%s final CSR offset is invalid in affected graph %d.",
                 label, graph_position);
    }
    for (R_xlen_t i = 0; i < col_idx.size(); ++i) {
      const int neighbor = col_idx[i];
      if (neighbor == NA_INTEGER || neighbor < 1 || neighbor > node_count) {
        Rcpp::stop("%s contains an invalid neighbor row in affected graph %d.",
                   label, graph_position);
      }
    }
  };

  auto append_edge = [&](int patch_from, int patch_to) {
    const int low = std::min(patch_from, patch_to);
    const int high = std::max(patch_from, patch_to);
    const std::uint64_t packed =
      (static_cast<std::uint64_t>(static_cast<std::uint32_t>(low)) << 32) |
      static_cast<std::uint32_t>(high);
    if (seen_edges.insert(packed).second) {
      edge_u.push_back(low);
      edge_v.push_back(high);
      candidate_patch_ids.push_back(low);
      candidate_patch_ids.push_back(high);
    }
  };

  auto reserve_edge_capacity = [&](std::size_t additional_edges) {
    if (additional_edges > edge_u.max_size() - edge_u.size() ||
        additional_edges > edge_v.max_size() - edge_v.size() ||
        additional_edges > (candidate_patch_ids.max_size() -
                            candidate_patch_ids.size()) / 2) {
      Rcpp::stop("Distance recheck-edge output exceeds vector allocation limits.");
    }
    const std::size_t requested_edges = edge_u.size() + additional_edges;
    edge_u.reserve(requested_edges);
    edge_v.reserve(requested_edges);
    seen_edges.reserve(requested_edges);
    candidate_patch_ids.reserve(
      candidate_patch_ids.size() + 2 * additional_edges
    );
  };

  for (R_xlen_t graph_index = 0; graph_index < affected_graphs.size(); ++graph_index) {
    if ((graph_index & 255) == 0) Rcpp::checkUserInterrupt();
    const int graph_position = static_cast<int>(graph_index + 1);
    if (Rf_isNull(affected_graphs[graph_index]) ||
        TYPEOF(affected_graphs[graph_index]) != VECSXP) {
      Rcpp::stop("Affected graph %d is not a graph list.", graph_position);
    }
    const Rcpp::List graph = Rcpp::as<Rcpp::List>(affected_graphs[graph_index]);
    if (!graph.containsElementNamed("id2patch")) {
      Rcpp::stop("Affected graph %d is missing id2patch.", graph_position);
    }
    const Rcpp::IntegerVector id2patch = graph["id2patch"];
    validate_patch_ids(id2patch, graph_position, "Effective id2patch");
    const int node_count = static_cast<int>(id2patch.size());
    graph_nodes_scanned += static_cast<double>(node_count);

    std::vector<int> recheck_rows;
    for (int row = 1; row <= node_count; ++row) {
      if (recheck_patch_set.count(id2patch[row - 1])) recheck_rows.push_back(row);
    }
    recheck_rows_total += static_cast<double>(recheck_rows.size());

    bool lazy = false;
    if (graph.containsElementNamed("graph_type")) {
      SEXP graph_type = graph["graph_type"];
      if (!Rf_isNull(graph_type)) {
        lazy = Rcpp::as<std::string>(graph_type) == "lazy_added_node_overlay";
      }
    }

    if (!lazy) {
      if (!graph.containsElementNamed("row_ptr") ||
          !graph.containsElementNamed("col_idx")) {
        Rcpp::stop("Standard affected graph %d is missing CSR fields.", graph_position);
      }
      const Rcpp::IntegerVector row_ptr = graph["row_ptr"];
      const Rcpp::IntegerVector col_idx = graph["col_idx"];
      validate_csr(row_ptr, col_idx, node_count, graph_position, "Standard CSR");
      if (recheck_rows.empty()) continue;
      std::size_t candidate_entries = 0;
      for (int row : recheck_rows) {
        candidate_entries += static_cast<std::size_t>(
          row_ptr[row] - row_ptr[row - 1]
        );
      }
      reserve_edge_capacity(candidate_entries);
      for (int row : recheck_rows) {
        const int first = row_ptr[row - 1];
        const int last = row_ptr[row];
        csr_entries_visited += static_cast<double>(last - first);
        for (int edge_index = first; edge_index < last; ++edge_index) {
          append_edge(id2patch[row - 1], id2patch[col_idx[edge_index] - 1]);
        }
      }
      continue;
    }

    if (!graph.containsElementNamed("base_graph") ||
        !graph.containsElementNamed("base_node_count") ||
        !graph.containsElementNamed("overlay_edges")) {
      Rcpp::stop("Lazy affected graph %d is missing overlay fields.", graph_position);
    }
    const int base_node_count = Rcpp::as<int>(graph["base_node_count"]);
    if (base_node_count < 0 || base_node_count > node_count) {
      Rcpp::stop("Lazy base-node count is invalid in affected graph %d.", graph_position);
    }
    const Rcpp::List base_graph = graph["base_graph"];
    if (!base_graph.containsElementNamed("id2patch") ||
        !base_graph.containsElementNamed("row_ptr") ||
        !base_graph.containsElementNamed("col_idx")) {
      Rcpp::stop("Lazy base graph %d is missing CSR fields.", graph_position);
    }
    const Rcpp::IntegerVector base_patch_ids = base_graph["id2patch"];
    const Rcpp::IntegerVector base_row_ptr = base_graph["row_ptr"];
    const Rcpp::IntegerVector base_col_idx = base_graph["col_idx"];
    validate_patch_ids(base_patch_ids, graph_position, "Lazy base id2patch");
    if (base_patch_ids.size() != base_node_count) {
      Rcpp::stop("Lazy base graph size is invalid in affected graph %d.", graph_position);
    }
    for (int i = 0; i < base_node_count; ++i) {
      if (base_patch_ids[i] != id2patch[i]) {
        Rcpp::stop("Lazy base graph does not match the effective prefix in affected graph %d.",
                   graph_position);
      }
    }
    validate_csr(
      base_row_ptr, base_col_idx, base_node_count, graph_position, "Lazy base CSR"
    );

    const Rcpp::List overlay = graph["overlay_edges"];
    if (!overlay.containsElementNamed("u") || !overlay.containsElementNamed("v")) {
      Rcpp::stop("Lazy overlay table is malformed in affected graph %d.", graph_position);
    }
    const Rcpp::IntegerVector overlay_u = overlay["u"];
    const Rcpp::IntegerVector overlay_v = overlay["v"];
    if (overlay_u.size() != overlay_v.size()) {
      Rcpp::stop("Lazy overlay columns have unequal lengths in affected graph %d.",
                 graph_position);
    }
    overlay_edges_scanned += static_cast<double>(overlay_u.size());

    std::unordered_set<int> recheck_row_set(recheck_rows.begin(), recheck_rows.end());
    std::unordered_map<int, std::vector<int> > overlay_neighbors;
    overlay_neighbors.reserve(recheck_rows.size());
    for (R_xlen_t i = 0; i < overlay_u.size(); ++i) {
      const int u = overlay_u[i];
      const int v = overlay_v[i];
      if (u == NA_INTEGER || v == NA_INTEGER || u < 1 || v < 1 ||
          u > node_count || v > node_count) {
        Rcpp::stop("Lazy overlay contains an invalid endpoint in affected graph %d.",
                   graph_position);
      }
      if (recheck_row_set.count(u)) overlay_neighbors[u].push_back(v);
      if (recheck_row_set.count(v)) overlay_neighbors[v].push_back(u);
    }

    std::size_t candidate_entries = 0;
    for (int row : recheck_rows) {
      if (row <= base_node_count) {
        candidate_entries += static_cast<std::size_t>(
          base_row_ptr[row] - base_row_ptr[row - 1]
        );
      }
      const auto overlay_it = overlay_neighbors.find(row);
      if (overlay_it != overlay_neighbors.end()) {
        candidate_entries += overlay_it->second.size();
      }
    }
    reserve_edge_capacity(candidate_entries);

    for (int row : recheck_rows) {
      std::vector<int> neighbors;
      if (row <= base_node_count) {
        const int first = base_row_ptr[row - 1];
        const int last = base_row_ptr[row];
        csr_entries_visited += static_cast<double>(last - first);
        neighbors.insert(
          neighbors.end(), base_col_idx.begin() + first, base_col_idx.begin() + last
        );
      }
      const auto overlay_it = overlay_neighbors.find(row);
      if (overlay_it != overlay_neighbors.end()) {
        neighbors.insert(
          neighbors.end(), overlay_it->second.begin(), overlay_it->second.end()
        );
      }
      std::sort(neighbors.begin(), neighbors.end());
      neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
      for (int neighbor : neighbors) {
        append_edge(id2patch[row - 1], id2patch[neighbor - 1]);
      }
    }
  }

  std::sort(candidate_patch_ids.begin(), candidate_patch_ids.end());
  candidate_patch_ids.erase(
    std::unique(candidate_patch_ids.begin(), candidate_patch_ids.end()),
    candidate_patch_ids.end()
  );

  return Rcpp::List::create(
    Rcpp::Named("patch_u") = Rcpp::wrap(edge_u),
    Rcpp::Named("patch_v") = Rcpp::wrap(edge_v),
    Rcpp::Named("candidate_patch_ids") = Rcpp::wrap(candidate_patch_ids),
    Rcpp::Named("graphs_inspected") = static_cast<double>(affected_graphs.size()),
    Rcpp::Named("graph_nodes_scanned") = graph_nodes_scanned,
    Rcpp::Named("recheck_rows") = recheck_rows_total,
    Rcpp::Named("csr_adjacency_entries_visited") = csr_entries_visited,
    Rcpp::Named("overlay_edges_scanned") = overlay_edges_scanned,
    Rcpp::Named("unique_recheck_edges") = static_cast<double>(edge_u.size())
  );
}


// Extract contiguous candidate segments from a patch-sorted compact index;
// output vectors are allocated once after their exact combined size is known.
// [[Rcpp::export]]
Rcpp::List stage6_extract_candidate_patch_cells_cpp(
    const Rcpp::IntegerVector& index_pid,
    const Rcpp::IntegerVector& index_cell,
    const Rcpp::IntegerVector& candidate_patch_ids_input) {
  if (index_pid.size() != index_cell.size()) {
    Rcpp::stop("Candidate-cell compact-index vectors must have equal lengths.");
  }

  const std::vector<int> candidate_patch_ids = sorted_unique(
    candidate_patch_ids_input,
    "Candidate distance patch IDs"
  );
  if (candidate_patch_ids.empty()) {
    return Rcpp::List::create(
      Rcpp::Named("patch_id") = Rcpp::IntegerVector(),
      Rcpp::Named("cell") = Rcpp::IntegerVector()
    );
  }
  if (index_pid.size() == 0) {
    Rcpp::stop("Candidate distance patch cells are missing from the compact index.");
  }

  int previous_patch_id = 0;
  for (R_xlen_t i = 0; i < index_pid.size(); ++i) {
    if ((i & 1048575) == 0) Rcpp::checkUserInterrupt();
    const int patch_id = index_pid[i];
    if (patch_id == NA_INTEGER || patch_id < 1 ||
        (i > 0 && patch_id < previous_patch_id)) {
      Rcpp::stop("Candidate-cell compact-index patch IDs must be sorted positive integers.");
    }
    previous_patch_id = patch_id;
  }

  const int* pid_begin = index_pid.begin();
  const int* pid_end = index_pid.end();
  std::vector<std::pair<R_xlen_t, R_xlen_t> > segments;
  segments.reserve(candidate_patch_ids.size());
  R_xlen_t output_size = 0;
  for (std::size_t i = 0; i < candidate_patch_ids.size(); ++i) {
    const int patch_id = candidate_patch_ids[i];
    const int* lower = std::lower_bound(pid_begin, pid_end, patch_id);
    const int* upper = std::upper_bound(lower, pid_end, patch_id);
    if (lower == upper) {
      Rcpp::stop("Candidate distance patch ID %d is missing from the compact index.", patch_id);
    }
    const R_xlen_t first = static_cast<R_xlen_t>(lower - pid_begin);
    const R_xlen_t last = static_cast<R_xlen_t>(upper - pid_begin);
    const R_xlen_t segment_size = last - first;
    if (segment_size > R_XLEN_T_MAX - output_size) {
      Rcpp::stop("Candidate-cell output length exceeds R vector storage.");
    }
    output_size += segment_size;
    segments.emplace_back(first, last);
  }

  Rcpp::IntegerVector output_pid = Rcpp::no_init(output_size);
  Rcpp::IntegerVector output_cell = Rcpp::no_init(output_size);
  std::unordered_set<int> selected_cells;
  selected_cells.reserve(static_cast<std::size_t>(output_size));
  R_xlen_t output_position = 0;
  for (std::size_t segment = 0; segment < segments.size(); ++segment) {
    if ((segment & 1023) == 0) Rcpp::checkUserInterrupt();
    for (R_xlen_t position = segments[segment].first;
         position < segments[segment].second;
         ++position) {
      const int cell = index_cell[position];
      if (cell == NA_INTEGER || cell < 1 || !selected_cells.insert(cell).second) {
        Rcpp::stop("Candidate cells must be unique positive integers.");
      }
      output_pid[output_position] = index_pid[position];
      output_cell[output_position] = cell;
      ++output_position;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("patch_id") = output_pid,
    Rcpp::Named("cell") = output_cell
  );
}


// Canonicalize PU and patch rows using stable spatial keys, returning an R-side
// one-based row permutation and dense one-based identifiers.
// [[Rcpp::export]]
Rcpp::List stage6_canonicalize_rank_lut_cpp(
    const Rcpp::IntegerVector& patch_ids,
    const Rcpp::IntegerVector& pu_ids,
    const Rcpp::IntegerVector& index_pid,
    const Rcpp::IntegerVector& index_cell) {
  const R_xlen_t n_patches = patch_ids.size();
  if (pu_ids.size() != n_patches) {
    Rcpp::stop("Canonical LUT patch and PU vectors must have equal lengths.");
  }
  if (index_pid.size() != index_cell.size()) {
    Rcpp::stop("Canonical LUT compact-index vectors must have equal lengths.");
  }
  if (n_patches == 0) {
    if (index_pid.size() != 0) {
      Rcpp::stop("Canonical LUT patch rows and compact index disagree.");
    }
    return Rcpp::List::create(
      Rcpp::Named("row_order") = Rcpp::IntegerVector(),
      Rcpp::Named("patch_id") = Rcpp::IntegerVector(),
      Rcpp::Named("pu_id") = Rcpp::IntegerVector(),
      Rcpp::Named("patch_n_cells") = Rcpp::IntegerVector()
    );
  }
  if (index_pid.size() == 0) {
    Rcpp::stop("Canonical LUT patch rows and compact index disagree.");
  }

  std::unordered_map<int, int> row_by_patch;
  row_by_patch.reserve(static_cast<std::size_t>(n_patches));
  for (R_xlen_t row = 0; row < n_patches; ++row) {
    const int patch_id = patch_ids[row];
    const int pu_id = pu_ids[row];
    if (patch_id == NA_INTEGER || patch_id < 1 ||
        !row_by_patch.emplace(patch_id, static_cast<int>(row)).second) {
      Rcpp::stop("Canonical LUT patch IDs must be unique positive integers.");
    }
    if (pu_id == NA_INTEGER || pu_id < 1) {
      Rcpp::stop("Canonical LUT PU IDs must be positive integers.");
    }
  }

  std::unordered_set<int> cells_seen;
  cells_seen.reserve(static_cast<std::size_t>(index_cell.size()));
  std::vector<int> patch_n_cells(static_cast<std::size_t>(n_patches), 0);
  std::vector<int> first_cell(static_cast<std::size_t>(n_patches), 0);
  std::vector<unsigned char> represented(static_cast<std::size_t>(n_patches), 0);
  R_xlen_t position = 0;
  int previous_patch_id = 0;
  while (position < index_pid.size()) {
    if ((position & 1048575) == 0) Rcpp::checkUserInterrupt();
    const int patch_id = index_pid[position];
    if (patch_id == NA_INTEGER || patch_id < 1 ||
        (position > 0 && patch_id < previous_patch_id)) {
      Rcpp::stop("Canonical LUT compact-index patch IDs must be sorted positive integers.");
    }
    const auto row_found = row_by_patch.find(patch_id);
    if (row_found == row_by_patch.end()) {
      Rcpp::stop("Canonical LUT compact index contains an unexpected patch ID.");
    }
    const int row = row_found->second;
    if (represented[static_cast<std::size_t>(row)]) {
      Rcpp::stop("Canonical LUT compact-index patch segments must be contiguous.");
    }
    represented[static_cast<std::size_t>(row)] = 1;
    const R_xlen_t start = position;
    const int segment_first_cell = index_cell[start];
    if (segment_first_cell == NA_INTEGER || segment_first_cell < 1) {
      Rcpp::stop("Canonical LUT compact-index cells must be positive integers.");
    }
    int minimum_cell = segment_first_cell;
    while (position < index_pid.size() && index_pid[position] == patch_id) {
      const int cell = index_cell[position];
      if (cell == NA_INTEGER || cell < 1 || !cells_seen.insert(cell).second) {
        Rcpp::stop("Canonical LUT compact-index cells must be unique positive integers.");
      }
      minimum_cell = std::min(minimum_cell, cell);
      ++position;
    }
    if (segment_first_cell != minimum_cell) {
      Rcpp::stop("Canonical LUT compact patch cells are not in canonical cell order.");
    }
    const R_xlen_t segment_length = position - start;
    if (segment_length > static_cast<R_xlen_t>(std::numeric_limits<int>::max())) {
      Rcpp::stop("Canonical LUT patch cell count exceeds integer storage.");
    }
    patch_n_cells[static_cast<std::size_t>(row)] =
      static_cast<int>(segment_length);
    first_cell[static_cast<std::size_t>(row)] = segment_first_cell;
    previous_patch_id = patch_id;
  }
  if (std::find(represented.begin(), represented.end(), 0) != represented.end()) {
    Rcpp::stop("Canonical LUT patch rows and compact index disagree.");
  }

  std::vector<int> scan_rows(static_cast<std::size_t>(n_patches));
  for (R_xlen_t row = 0; row < n_patches; ++row) {
    scan_rows[static_cast<std::size_t>(row)] = static_cast<int>(row);
  }
  std::sort(scan_rows.begin(), scan_rows.end(), [&first_cell](int a, int b) {
    return first_cell[static_cast<std::size_t>(a)] <
      first_cell[static_cast<std::size_t>(b)];
  });

  std::vector<int> scan_patch_id(static_cast<std::size_t>(n_patches), 0);
  std::unordered_map<int, int> first_scan_patch_by_pu;
  first_scan_patch_by_pu.reserve(static_cast<std::size_t>(n_patches));
  for (R_xlen_t scan = 0; scan < n_patches; ++scan) {
    const int row = scan_rows[static_cast<std::size_t>(scan)];
    const int scan_id = static_cast<int>(scan) + 1;
    scan_patch_id[static_cast<std::size_t>(row)] = scan_id;
    const int pu_id = pu_ids[row];
    const auto inserted = first_scan_patch_by_pu.emplace(pu_id, scan_id);
    if (!inserted.second && scan_id < inserted.first->second) {
      inserted.first->second = scan_id;
    }
  }
  std::vector<std::pair<int, int> > pu_order;
  pu_order.reserve(first_scan_patch_by_pu.size());
  for (const auto& entry : first_scan_patch_by_pu) {
    pu_order.emplace_back(entry.second, entry.first);
  }
  std::sort(pu_order.begin(), pu_order.end());
  std::unordered_map<int, int> canonical_pu_by_internal;
  canonical_pu_by_internal.reserve(pu_order.size());
  for (std::size_t i = 0; i < pu_order.size(); ++i) {
    canonical_pu_by_internal.emplace(pu_order[i].second, static_cast<int>(i) + 1);
  }

  std::vector<int> final_rows(scan_rows);
  std::sort(final_rows.begin(), final_rows.end(), [
      &pu_ids, &canonical_pu_by_internal, &scan_patch_id](int a, int b) {
    const int pu_a = canonical_pu_by_internal.at(pu_ids[a]);
    const int pu_b = canonical_pu_by_internal.at(pu_ids[b]);
    if (pu_a != pu_b) return pu_a < pu_b;
    return scan_patch_id[static_cast<std::size_t>(a)] <
      scan_patch_id[static_cast<std::size_t>(b)];
  });

  Rcpp::IntegerVector row_order = Rcpp::no_init(n_patches);
  Rcpp::IntegerVector canonical_patch_id = Rcpp::no_init(n_patches);
  Rcpp::IntegerVector canonical_pu_id = Rcpp::no_init(n_patches);
  Rcpp::IntegerVector canonical_patch_n_cells = Rcpp::no_init(n_patches);
  for (R_xlen_t output = 0; output < n_patches; ++output) {
    const int row = final_rows[static_cast<std::size_t>(output)];
    row_order[output] = row + 1;
    canonical_patch_id[output] = static_cast<int>(output) + 1;
    canonical_pu_id[output] = canonical_pu_by_internal.at(pu_ids[row]);
    canonical_patch_n_cells[output] =
      patch_n_cells[static_cast<std::size_t>(row)];
  }
  return Rcpp::List::create(
    Rcpp::Named("row_order") = row_order,
    Rcpp::Named("patch_id") = canonical_patch_id,
    Rcpp::Named("pu_id") = canonical_pu_id,
    Rcpp::Named("patch_n_cells") = canonical_patch_n_cells
  );
}


// ---- Incremental distance and fragmentation graphs -------------------------

// Remove invalid undirected distance edges and merge the validated overlay;
// adjacency rows are sorted before the one-based CSR neighbor export.
// [[Rcpp::export]]
Rcpp::List stage6_filter_distance_graph_cpp(
    const Rcpp::IntegerVector& id2patch,
    int base_node_count,
    const Rcpp::IntegerVector& base_row_ptr,
    const Rcpp::IntegerVector& base_col_idx,
    const Rcpp::IntegerVector& overlay_u,
    const Rcpp::IntegerVector& overlay_v,
    const Rcpp::IntegerVector& invalid_patch_u,
    const Rcpp::IntegerVector& invalid_patch_v,
    int pu_id) {
  const R_xlen_t n_nodes = id2patch.size();
  if (pu_id == NA_INTEGER || pu_id < 1) {
    Rcpp::stop("The distance-filter PU ID must be a positive integer.");
  }
  if (n_nodes > static_cast<R_xlen_t>(std::numeric_limits<int>::max()) ||
      base_node_count < 0 || base_node_count > n_nodes) {
    Rcpp::stop("The distance-filter graph has an invalid node count.");
  }
  if (base_row_ptr.size() != static_cast<R_xlen_t>(base_node_count) + 1 ||
      base_row_ptr[0] != 0 || base_row_ptr[base_node_count] != base_col_idx.size()) {
    Rcpp::stop("The distance-filter graph has an invalid CSR row pointer.");
  }
  if (overlay_u.size() != overlay_v.size() ||
      invalid_patch_u.size() != invalid_patch_v.size()) {
    Rcpp::stop("Distance-filter edge endpoint vectors must have equal lengths.");
  }

  std::unordered_set<int> patch_seen;
  patch_seen.reserve(static_cast<std::size_t>(n_nodes));
  for (R_xlen_t i = 0; i < n_nodes; ++i) {
    const int patch_id = id2patch[i];
    if (patch_id == NA_INTEGER || patch_id < 1 ||
        !patch_seen.insert(patch_id).second) {
      Rcpp::stop("Distance-filter patch IDs must be unique positive integers.");
    }
  }
  int previous_offset = 0;
  for (int i = 0; i < base_node_count; ++i) {
    const int offset = base_row_ptr[i + 1];
    if (offset == NA_INTEGER || offset < previous_offset ||
        offset > base_col_idx.size()) {
      Rcpp::stop("The distance-filter graph has malformed CSR offsets.");
    }
    previous_offset = offset;
  }
  for (R_xlen_t i = 0; i < base_col_idx.size(); ++i) {
    const int neighbor = base_col_idx[i];
    if (neighbor == NA_INTEGER || neighbor < 1 || neighbor > base_node_count) {
      Rcpp::stop("The distance-filter graph contains an invalid base neighbor row.");
    }
  }
  for (R_xlen_t i = 0; i < overlay_u.size(); ++i) {
    if (overlay_u[i] == NA_INTEGER || overlay_v[i] == NA_INTEGER ||
        overlay_u[i] < 1 || overlay_u[i] > n_nodes ||
        overlay_v[i] < 1 || overlay_v[i] > n_nodes) {
      Rcpp::stop("The distance-filter graph contains an invalid overlay edge.");
    }
  }

  const auto pack_pair = [](std::uint32_t a, std::uint32_t b) {
    const std::uint32_t low = std::min(a, b);
    const std::uint32_t high = std::max(a, b);
    return (static_cast<std::uint64_t>(low) << 32) |
      static_cast<std::uint64_t>(high);
  };
  std::unordered_set<std::uint64_t> invalid_edges;
  invalid_edges.reserve(static_cast<std::size_t>(invalid_patch_u.size()));
  for (R_xlen_t i = 0; i < invalid_patch_u.size(); ++i) {
    const int patch_u = invalid_patch_u[i];
    const int patch_v = invalid_patch_v[i];
    if (patch_u == NA_INTEGER || patch_v == NA_INTEGER ||
        patch_u < 1 || patch_v < 1 || patch_u == patch_v ||
        patch_seen.find(patch_u) == patch_seen.end() ||
        patch_seen.find(patch_v) == patch_seen.end()) {
      Rcpp::stop("Invalid distance-filter patch pair.");
    }
    invalid_edges.insert(pack_pair(
      static_cast<std::uint32_t>(patch_u),
      static_cast<std::uint32_t>(patch_v)
    ));
  }
  std::unordered_map<std::uint64_t, unsigned char> invalid_edge_directions;
  invalid_edge_directions.reserve(invalid_edges.size());

  std::vector<std::uint64_t> retained_edges;
  retained_edges.reserve(
    static_cast<std::size_t>(base_col_idx.size() / 2 + overlay_u.size())
  );
  for (int u = 0; u < base_node_count; ++u) {
    if ((u & 4095) == 0) Rcpp::checkUserInterrupt();
    for (int edge = base_row_ptr[u]; edge < base_row_ptr[u + 1]; ++edge) {
      const int v = base_col_idx[edge] - 1;
      const std::uint64_t patch_pair = pack_pair(
        static_cast<std::uint32_t>(id2patch[u]),
        static_cast<std::uint32_t>(id2patch[v])
      );
      if (invalid_edges.find(patch_pair) != invalid_edges.end()) {
        const bool forward = id2patch[u] < id2patch[v];
        invalid_edge_directions[patch_pair] |= forward ? 1U : 2U;
      } else if (v > u) {
        retained_edges.push_back(pack_pair(
          static_cast<std::uint32_t>(u),
          static_cast<std::uint32_t>(v)
        ));
      }
    }
  }
  for (R_xlen_t i = 0; i < overlay_u.size(); ++i) {
    const int u = overlay_u[i] - 1;
    const int v = overlay_v[i] - 1;
    const std::uint64_t patch_pair = pack_pair(
      static_cast<std::uint32_t>(id2patch[u]),
      static_cast<std::uint32_t>(id2patch[v])
    );
    if (invalid_edges.find(patch_pair) != invalid_edges.end()) {
      invalid_edge_directions[patch_pair] = 3U;
    } else {
      retained_edges.push_back(pack_pair(
        static_cast<std::uint32_t>(u),
        static_cast<std::uint32_t>(v)
      ));
    }
  }
  for (std::uint64_t invalid_edge : invalid_edges) {
    const auto observed = invalid_edge_directions.find(invalid_edge);
    if (observed == invalid_edge_directions.end() || observed->second != 3U) {
      Rcpp::stop(
        "A requested invalid distance edge is absent or not stored symmetrically."
      );
    }
  }

  std::sort(retained_edges.begin(), retained_edges.end());
  retained_edges.erase(
    std::unique(retained_edges.begin(), retained_edges.end()),
    retained_edges.end()
  );
  if (retained_edges.size() > static_cast<std::size_t>(
        std::numeric_limits<int>::max() / 2)) {
    Rcpp::stop("The filtered distance graph is too large for integer CSR storage.");
  }

  std::vector<int> degree(static_cast<std::size_t>(n_nodes), 0);
  for (std::uint64_t packed : retained_edges) {
    const std::uint32_t u = static_cast<std::uint32_t>(packed >> 32);
    const std::uint32_t v = static_cast<std::uint32_t>(packed & 0xffffffffULL);
    ++degree[u];
    ++degree[v];
  }
  std::vector<int> row_ptr(static_cast<std::size_t>(n_nodes) + 1, 0);
  for (R_xlen_t i = 0; i < n_nodes; ++i) {
    row_ptr[static_cast<std::size_t>(i) + 1] =
      row_ptr[static_cast<std::size_t>(i)] + degree[static_cast<std::size_t>(i)];
  }
  std::vector<int> col_idx(static_cast<std::size_t>(row_ptr.back()), 0);
  std::vector<int> cursor(row_ptr.begin(), row_ptr.end() - 1);
  for (std::uint64_t packed : retained_edges) {
    const std::uint32_t u = static_cast<std::uint32_t>(packed >> 32);
    const std::uint32_t v = static_cast<std::uint32_t>(packed & 0xffffffffULL);
    col_idx[static_cast<std::size_t>(cursor[u]++)] = static_cast<int>(v) + 1;
    col_idx[static_cast<std::size_t>(cursor[v]++)] = static_cast<int>(u) + 1;
  }
  for (R_xlen_t i = 0; i < n_nodes; ++i) {
    std::sort(
      col_idx.begin() + row_ptr[static_cast<std::size_t>(i)],
      col_idx.begin() + row_ptr[static_cast<std::size_t>(i) + 1]
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("pu_id") = pu_id,
    Rcpp::Named("id2patch") = Rcpp::clone(id2patch),
    Rcpp::Named("row_ptr") = Rcpp::wrap(row_ptr),
    Rcpp::Named("col_idx") = Rcpp::wrap(col_idx),
    Rcpp::Named("invalid_undirected_edges") =
      static_cast<double>(invalid_edges.size()),
    Rcpp::Named("input_adjacency_entries") =
      static_cast<double>(base_col_idx.size() + 2 * overlay_u.size()),
    Rcpp::Named("output_adjacency_entries") =
      static_cast<double>(col_idx.size())
  );
}


// Expand the prior graph through fragment ancestry, deduplicate packed edges,
// and export a deterministic one-based-neighbor CSR provisional graph.
// [[Rcpp::export]]
Rcpp::List stage6_build_provisional_fragment_graph_cpp(
    const Rcpp::IntegerVector& current_patch_ids,
    const Rcpp::IntegerVector& current_origin_patch_ids,
    const Rcpp::IntegerVector& old_patch_ids,
    const Rcpp::IntegerVector& old_row_ptr,
    const Rcpp::IntegerVector& old_col_idx,
    int pu_id) {
  const R_xlen_t n_current = current_patch_ids.size();
  const R_xlen_t n_old = old_patch_ids.size();
  if (current_origin_patch_ids.size() != n_current) {
    Rcpp::stop("Current patch and origin vectors must have equal lengths.");
  }
  if (pu_id == NA_INTEGER || pu_id < 1) {
    Rcpp::stop("The provisional graph PU ID must be a positive integer.");
  }
  if (n_current > static_cast<R_xlen_t>(std::numeric_limits<int>::max())) {
    Rcpp::stop("The provisional graph has too many nodes for integer CSR storage.");
  }

  std::unordered_set<int> current_patch_ids_seen;
  current_patch_ids_seen.reserve(static_cast<std::size_t>(n_current));
  std::unordered_map<int, std::vector<std::uint32_t> > current_rows_by_origin;
  current_rows_by_origin.reserve(static_cast<std::size_t>(n_current));
  std::vector<int> origin_order;
  origin_order.reserve(static_cast<std::size_t>(n_current));
  std::unordered_set<int> origins_seen;
  origins_seen.reserve(static_cast<std::size_t>(n_current));
  for (R_xlen_t i = 0; i < n_current; ++i) {
    const int patch_id = current_patch_ids[i];
    const int origin_id = current_origin_patch_ids[i];
    if (patch_id == NA_INTEGER || patch_id < 1 ||
        !current_patch_ids_seen.insert(patch_id).second) {
      Rcpp::stop("Current provisional patch IDs must be unique positive integers.");
    }
    if (origin_id == NA_INTEGER) {
      Rcpp::stop("Current provisional origin patch IDs must not be missing.");
    }
    current_rows_by_origin[origin_id].push_back(static_cast<std::uint32_t>(i));
    if (origins_seen.insert(origin_id).second) origin_order.push_back(origin_id);
  }

  if (old_row_ptr.size() != n_old + 1 || old_row_ptr[0] != 0) {
    Rcpp::stop("The old provisional graph has an invalid CSR row pointer.");
  }
  std::unordered_set<int> old_patch_ids_seen;
  old_patch_ids_seen.reserve(static_cast<std::size_t>(n_old));
  int previous_offset = 0;
  for (R_xlen_t i = 0; i < n_old; ++i) {
    const int patch_id = old_patch_ids[i];
    if (patch_id == NA_INTEGER || patch_id < 1 ||
        !old_patch_ids_seen.insert(patch_id).second) {
      Rcpp::stop("Old provisional graph patch IDs must be unique positive integers.");
    }
    const int offset = old_row_ptr[i + 1];
    if (offset == NA_INTEGER || offset < previous_offset ||
        offset > old_col_idx.size()) {
      Rcpp::stop("The old provisional graph has malformed CSR offsets.");
    }
    previous_offset = offset;
  }
  if (old_row_ptr[n_old] != old_col_idx.size()) {
    Rcpp::stop("The old provisional graph CSR offset does not match its adjacency length.");
  }
  for (R_xlen_t edge = 0; edge < old_col_idx.size(); ++edge) {
    const int neighbor = old_col_idx[edge];
    if (neighbor == NA_INTEGER || neighbor < 1 || neighbor > n_old) {
      Rcpp::stop("The old provisional graph contains an invalid neighbor row.");
    }
  }

  const auto checked_add = [](std::uint64_t total, std::uint64_t increment) {
    if (increment > std::numeric_limits<std::uint64_t>::max() - total) {
      Rcpp::stop("Provisional graph candidate-edge count overflowed.");
    }
    return total + increment;
  };

  // Count candidate edges first so the packed edge vector is allocated once.
  std::uint64_t candidate_count = 0;
  for (R_xlen_t old_u = 0; old_u < n_old; ++old_u) {
    if ((old_u & 4095) == 0) Rcpp::checkUserInterrupt();
    const auto descendants_u = current_rows_by_origin.find(old_patch_ids[old_u]);
    if (descendants_u == current_rows_by_origin.end()) continue;
    const int first = old_row_ptr[old_u];
    const int last = old_row_ptr[old_u + 1];
    for (int edge = first; edge < last; ++edge) {
      const int old_v = old_col_idx[edge] - 1;
      if (old_v <= old_u) continue;
      const auto descendants_v = current_rows_by_origin.find(old_patch_ids[old_v]);
      if (descendants_v == current_rows_by_origin.end()) continue;
      const std::uint64_t n_pairs =
        static_cast<std::uint64_t>(descendants_u->second.size()) *
        static_cast<std::uint64_t>(descendants_v->second.size());
      candidate_count = checked_add(candidate_count, n_pairs);
    }
  }
  for (int origin_id : origin_order) {
    const std::uint64_t group_size = static_cast<std::uint64_t>(
      current_rows_by_origin[origin_id].size()
    );
    candidate_count = checked_add(
      candidate_count,
      group_size < 2 ? 0 : group_size * (group_size - 1) / 2
    );
  }
  if (candidate_count > static_cast<std::uint64_t>(
        std::vector<std::uint64_t>().max_size())) {
    Rcpp::stop("The provisional graph candidate-edge vector is too large to allocate.");
  }

  std::vector<std::uint64_t> packed_edges;
  packed_edges.reserve(static_cast<std::size_t>(candidate_count));
  const auto append_edge = [&packed_edges](std::uint32_t u, std::uint32_t v) {
    if (u == v) return;
    const std::uint32_t low = std::min(u, v);
    const std::uint32_t high = std::max(u, v);
    packed_edges.push_back(
      (static_cast<std::uint64_t>(low) << 32) |
      static_cast<std::uint64_t>(high)
    );
  };

  // Expand inherited edges through the current descendants of each endpoint.
  for (R_xlen_t old_u = 0; old_u < n_old; ++old_u) {
    if ((old_u & 4095) == 0) Rcpp::checkUserInterrupt();
    const auto descendants_u = current_rows_by_origin.find(old_patch_ids[old_u]);
    if (descendants_u == current_rows_by_origin.end()) continue;
    const int first = old_row_ptr[old_u];
    const int last = old_row_ptr[old_u + 1];
    for (int edge = first; edge < last; ++edge) {
      const int old_v = old_col_idx[edge] - 1;
      if (old_v <= old_u) continue;
      const auto descendants_v = current_rows_by_origin.find(old_patch_ids[old_v]);
      if (descendants_v == current_rows_by_origin.end()) continue;
      for (std::uint32_t current_u : descendants_u->second) {
        for (std::uint32_t current_v : descendants_v->second) {
          append_edge(current_u, current_v);
        }
      }
    }
  }

  // Add provisional sibling edges among fragments of the same origin patch.
  for (int origin_id : origin_order) {
    const std::vector<std::uint32_t>& rows = current_rows_by_origin[origin_id];
    for (std::size_t i = 0; i < rows.size(); ++i) {
      for (std::size_t j = i + 1; j < rows.size(); ++j) {
        append_edge(rows[i], rows[j]);
      }
    }
  }
  if (packed_edges.size() != candidate_count) {
    Rcpp::stop("The provisional graph candidate-edge count is internally inconsistent.");
  }

  std::sort(packed_edges.begin(), packed_edges.end());
  packed_edges.erase(
    std::unique(packed_edges.begin(), packed_edges.end()),
    packed_edges.end()
  );
  const std::uint64_t unique_edges = packed_edges.size();
  if (unique_edges > static_cast<std::uint64_t>(
        std::numeric_limits<int>::max() / 2)) {
    Rcpp::stop("The provisional graph adjacency is too large for integer CSR storage.");
  }

  std::vector<int> degree(static_cast<std::size_t>(n_current), 0);
  for (std::uint64_t packed : packed_edges) {
    const std::uint32_t u = static_cast<std::uint32_t>(packed >> 32);
    const std::uint32_t v = static_cast<std::uint32_t>(packed & 0xffffffffULL);
    ++degree[u];
    ++degree[v];
  }
  std::vector<int> row_ptr(static_cast<std::size_t>(n_current) + 1, 0);
  for (R_xlen_t i = 0; i < n_current; ++i) {
    row_ptr[static_cast<std::size_t>(i) + 1] =
      row_ptr[static_cast<std::size_t>(i)] + degree[static_cast<std::size_t>(i)];
  }
  std::vector<int> col_idx(static_cast<std::size_t>(row_ptr.back()), 0);
  std::vector<int> cursor(row_ptr.begin(), row_ptr.end() - 1);
  for (std::uint64_t packed : packed_edges) {
    const std::uint32_t u = static_cast<std::uint32_t>(packed >> 32);
    const std::uint32_t v = static_cast<std::uint32_t>(packed & 0xffffffffULL);
    col_idx[static_cast<std::size_t>(cursor[u]++)] = static_cast<int>(v) + 1;
    col_idx[static_cast<std::size_t>(cursor[v]++)] = static_cast<int>(u) + 1;
  }
  for (R_xlen_t i = 0; i < n_current; ++i) {
    std::sort(
      col_idx.begin() + row_ptr[static_cast<std::size_t>(i)],
      col_idx.begin() + row_ptr[static_cast<std::size_t>(i) + 1]
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("pu_id") = pu_id,
    Rcpp::Named("id2patch") = Rcpp::clone(current_patch_ids),
    Rcpp::Named("row_ptr") = Rcpp::wrap(row_ptr),
    Rcpp::Named("col_idx") = Rcpp::wrap(col_idx),
    Rcpp::Named("candidate_undirected_edges") =
      static_cast<double>(candidate_count),
    Rcpp::Named("unique_undirected_edges") =
      static_cast<double>(unique_edges),
    Rcpp::Named("output_adjacency_entries") =
      static_cast<double>(col_idx.size())
  );
}


// Label rook-connected components from one-based raster cells. Input cells are
// sorted and deduplicated before traversal, fixing component discovery order.
// [[Rcpp::export]]
Rcpp::List stage6_label_patch_components_cpp(
    const Rcpp::IntegerVector& patch_cells_input,
    const Rcpp::IntegerVector& rook_start,
    const Rcpp::IntegerVector& rook_end,
    const Rcpp::IntegerVector& rook_to,
    int n_cells) {
  if (n_cells < 1 || rook_start.size() != n_cells || rook_end.size() != n_cells) {
    Rcpp::stop("Stage 6 component kernel received an invalid rook index.");
  }

  std::vector<int> patch_cells;
  patch_cells.reserve(static_cast<std::size_t>(patch_cells_input.size()));
  for (int cell : patch_cells_input) {
    if (cell != NA_INTEGER && cell >= 1 && cell <= n_cells) {
      patch_cells.push_back(cell);
    }
  }
  std::sort(patch_cells.begin(), patch_cells.end());
  patch_cells.erase(
    std::unique(patch_cells.begin(), patch_cells.end()),
    patch_cells.end()
  );

  const std::size_t n_patch_cells = patch_cells.size();
  if (n_patch_cells == 0) {
    return Rcpp::List::create(
      Rcpp::Named("cells") = Rcpp::IntegerVector(),
      Rcpp::Named("component_id") = Rcpp::IntegerVector(),
      Rcpp::Named("component_count") = 0
    );
  }

  std::unordered_map<int, int> position_by_cell;
  position_by_cell.reserve(n_patch_cells);
  for (std::size_t i = 0; i < n_patch_cells; ++i) {
    position_by_cell[patch_cells[i]] = static_cast<int>(i);
  }

  std::vector<int> component_id(n_patch_cells, 0);
  std::vector<int> queue(n_patch_cells, 0);
  int component_count = 0;

  for (std::size_t seed = 0; seed < n_patch_cells; ++seed) {
    if (component_id[seed] != 0) continue;
    ++component_count;
    std::size_t head = 0;
    std::size_t tail = 1;
    queue[0] = static_cast<int>(seed);
    component_id[seed] = component_count;

    while (head < tail) {
      const int current_position = queue[head++];
      const int current_cell = patch_cells[static_cast<std::size_t>(current_position)];
      const int start = rook_start[current_cell - 1];
      if (start == NA_INTEGER) continue;
      const int end = rook_end[current_cell - 1];
      if (end == NA_INTEGER || start < 1 || end < start || end > rook_to.size()) {
        Rcpp::stop("Stage 6 component kernel received malformed rook offsets.");
      }
      for (int edge = start - 1; edge < end; ++edge) {
        const int neighbor_cell = rook_to[edge];
        const auto found = position_by_cell.find(neighbor_cell);
        if (found == position_by_cell.end()) continue;
        const int neighbor_position = found->second;
        if (component_id[static_cast<std::size_t>(neighbor_position)] != 0) continue;
        component_id[static_cast<std::size_t>(neighbor_position)] = component_count;
        queue[tail++] = neighbor_position;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("cells") = Rcpp::wrap(patch_cells),
    Rcpp::Named("component_id") = Rcpp::wrap(component_id),
    Rcpp::Named("component_count") = component_count
  );
}


static Rcpp::List analyze_fragmentation_components_impl(
    const Rcpp::IntegerVector& index_pid,
    const Rcpp::IntegerVector& index_cell,
    const Rcpp::IntegerVector& changed_patch_ids,
    const Rcpp::IntegerVector& current_patch_id_by_cell,
    const Rcpp::NumericVector& cell_area_by_cell,
    const Rcpp::IntegerVector& rook_start,
    const Rcpp::IntegerVector& rook_end,
    const Rcpp::IntegerVector& rook_to,
    int n_cells,
    bool profile,
    bool compact_state_only) {
  using clock_type = std::chrono::steady_clock;
  typedef clock_type::time_point time_point;
  const auto elapsed_seconds = [](const time_point& started) {
    return std::chrono::duration<double>(clock_type::now() - started).count();
  };

  if (n_cells < 1 ||
      (!compact_state_only && current_patch_id_by_cell.size() != n_cells) ||
      cell_area_by_cell.size() != n_cells || rook_start.size() != n_cells ||
      rook_end.size() != n_cells) {
    Rcpp::stop("Fragmentation analysis received inconsistent raster-vector lengths.");
  }
  if (index_pid.size() != index_cell.size()) {
    Rcpp::stop("Fragmentation compact-index vectors must have equal lengths.");
  }

  // Stage 6 retains the complete serialized-state validation. Stage 7.2 has
  // already validated its authoritative compact index and validates only the
  // changed segments below, avoiding a full live-index scan per touched species.
  if (!compact_state_only) {
    int previous_pid = 0;
    std::unordered_set<int> indexed_cells;
    indexed_cells.reserve(static_cast<std::size_t>(index_cell.size()));
    for (R_xlen_t i = 0; i < index_pid.size(); ++i) {
      if ((i & 1048575) == 0) Rcpp::checkUserInterrupt();
      const int pid = index_pid[i];
      const int cell = index_cell[i];
      if (pid == NA_INTEGER || pid < 1 || (i > 0 && pid < previous_pid)) {
        Rcpp::stop("Fragmentation compact-index patch IDs must be sorted positive integers.");
      }
      if (cell == NA_INTEGER || cell < 1 || cell > n_cells ||
          !indexed_cells.insert(cell).second) {
        Rcpp::stop("Fragmentation compact-index cells must be unique valid positive integers.");
      }
      previous_pid = pid;
    }
  }

  std::unordered_set<int> changed_seen;
  changed_seen.reserve(static_cast<std::size_t>(changed_patch_ids.size()));
  for (int pid : changed_patch_ids) {
    if (pid == NA_INTEGER || pid < 1 || !changed_seen.insert(pid).second) {
      Rcpp::stop("Changed fragmentation patch IDs must be unique positive integers.");
    }
  }

  const int* pid_begin = index_pid.begin();
  const int* pid_end = index_pid.end();
  const R_xlen_t n_origins = changed_patch_ids.size();
  std::vector<std::pair<R_xlen_t, R_xlen_t> > segments;
  segments.reserve(static_cast<std::size_t>(n_origins));
  R_xlen_t maximum_output_size = 0;
  for (R_xlen_t origin = 0; origin < n_origins; ++origin) {
    const int patch_id = changed_patch_ids[origin];
    const int* lower = std::lower_bound(pid_begin, pid_end, patch_id);
    const int* upper = std::upper_bound(lower, pid_end, patch_id);
    if (lower == upper) {
      Rcpp::stop("Changed fragmentation patch ID %d is missing from the compact index.", patch_id);
    }
    const R_xlen_t first = static_cast<R_xlen_t>(lower - pid_begin);
    const R_xlen_t last = static_cast<R_xlen_t>(upper - pid_begin);
    if (compact_state_only) {
      if ((first > 0 && index_pid[first - 1] >= patch_id) ||
          (last < index_pid.size() && index_pid[last] <= patch_id)) {
        Rcpp::stop(
          "Changed fragmentation patch segment violates compact-index ordering."
        );
      }
    }
    if (last - first > R_XLEN_T_MAX - maximum_output_size) {
      Rcpp::stop("Fragmentation analysis output exceeds R vector storage.");
    }
    maximum_output_size += last - first;
    segments.emplace_back(first, last);
  }

  std::vector<int> output_cells;
  std::vector<int> output_component_id;
  output_cells.reserve(static_cast<std::size_t>(maximum_output_size));
  output_component_id.reserve(static_cast<std::size_t>(maximum_output_size));
  std::vector<int> origin_offsets(static_cast<std::size_t>(n_origins) + 1L, 0);
  std::vector<int> component_origin_patch_id;
  std::vector<int> component_id;
  std::vector<double> component_area_km2;
  std::vector<int> component_n_cells;
  std::vector<int> component_first_cell;
  std::vector<int> component_count_by_origin;
  component_count_by_origin.reserve(static_cast<std::size_t>(n_origins));

  double lookup_seconds = 0.0;
  double label_seconds = 0.0;
  double area_seconds = 0.0;

  for (R_xlen_t origin = 0; origin < n_origins; ++origin) {
    if ((origin & 255) == 0) Rcpp::checkUserInterrupt();
    const int patch_id = changed_patch_ids[origin];
    time_point subphase_started;
    if (profile) subphase_started = clock_type::now();

    std::vector<int> patch_cells;
    patch_cells.reserve(static_cast<std::size_t>(
      segments[static_cast<std::size_t>(origin)].second -
      segments[static_cast<std::size_t>(origin)].first
    ));
    for (R_xlen_t position = segments[static_cast<std::size_t>(origin)].first;
         position < segments[static_cast<std::size_t>(origin)].second;
         ++position) {
      const int cell = index_cell[position];
      if (index_pid[position] != patch_id || cell == NA_INTEGER ||
          cell < 1 || cell > n_cells) {
        Rcpp::stop(
          "Changed fragmentation patch segment contains an invalid entry."
        );
      }
      if (compact_state_only ||
          current_patch_id_by_cell[cell - 1] == patch_id) {
        patch_cells.push_back(cell);
      }
    }
    std::sort(patch_cells.begin(), patch_cells.end());
    const auto unique_end = std::unique(patch_cells.begin(), patch_cells.end());
    if (unique_end != patch_cells.end()) {
      Rcpp::stop("Fragmentation patch ID %d contains duplicate current cells.", patch_id);
    }
    if (patch_cells.empty()) {
      Rcpp::stop("No current cells remain for surviving fragmentation patch ID %d.", patch_id);
    }
    if (profile) lookup_seconds += elapsed_seconds(subphase_started);

    if (profile) subphase_started = clock_type::now();
    const std::size_t n_patch_cells = patch_cells.size();
    std::unordered_map<int, int> position_by_cell;
    position_by_cell.reserve(n_patch_cells);
    for (std::size_t i = 0; i < n_patch_cells; ++i) {
      position_by_cell[patch_cells[i]] = static_cast<int>(i);
    }
    std::vector<int> labels(n_patch_cells, 0);
    std::vector<int> queue(n_patch_cells, 0);
    int component_count = 0;
    for (std::size_t seed = 0; seed < n_patch_cells; ++seed) {
      if (labels[seed] != 0) continue;
      ++component_count;
      std::size_t head = 0;
      std::size_t tail = 1;
      queue[0] = static_cast<int>(seed);
      labels[seed] = component_count;
      while (head < tail) {
        const int current_position = queue[head++];
        const int current_cell = patch_cells[static_cast<std::size_t>(current_position)];
        const int start = rook_start[current_cell - 1];
        const int end = rook_end[current_cell - 1];
        if (start == NA_INTEGER) continue;
        if (end == NA_INTEGER || start < 1 || end < start || end > rook_to.size()) {
          Rcpp::stop("Fragmentation patch ID %d encountered malformed rook offsets.", patch_id);
        }
        for (int edge = start - 1; edge < end; ++edge) {
          const int neighbor_cell = rook_to[edge];
          if (neighbor_cell == NA_INTEGER || neighbor_cell < 1 || neighbor_cell > n_cells) {
            Rcpp::stop("Fragmentation patch ID %d encountered an invalid rook neighbor.", patch_id);
          }
          const auto found = position_by_cell.find(neighbor_cell);
          if (found == position_by_cell.end()) continue;
          const int neighbor_position = found->second;
          if (labels[static_cast<std::size_t>(neighbor_position)] != 0) continue;
          labels[static_cast<std::size_t>(neighbor_position)] = component_count;
          queue[tail++] = neighbor_position;
        }
      }
    }
    if (profile) label_seconds += elapsed_seconds(subphase_started);

    if (profile) subphase_started = clock_type::now();
    std::vector<double> areas(static_cast<std::size_t>(component_count), 0.0);
    std::vector<int> counts(static_cast<std::size_t>(component_count), 0);
    std::vector<int> first_cells(
      static_cast<std::size_t>(component_count),
      std::numeric_limits<int>::max()
    );
    for (std::size_t i = 0; i < n_patch_cells; ++i) {
      const int cell = patch_cells[i];
      const double area = cell_area_by_cell[cell - 1];
      if (!std::isfinite(area) || area <= 0.0) {
        Rcpp::stop("Fragmentation patch ID %d contains a nonpositive or nonfinite cell area.", patch_id);
      }
      const int component = labels[i] - 1;
      areas[static_cast<std::size_t>(component)] += area;
      counts[static_cast<std::size_t>(component)] += 1;
      if (cell < first_cells[static_cast<std::size_t>(component)]) {
        first_cells[static_cast<std::size_t>(component)] = cell;
      }
    }
    for (int component = 0; component < component_count; ++component) {
      const double area = areas[static_cast<std::size_t>(component)];
      if (!std::isfinite(area) || area <= 0.0) {
        Rcpp::stop("Fragmentation patch ID %d produced an invalid component area.", patch_id);
      }
      component_origin_patch_id.push_back(patch_id);
      component_id.push_back(component + 1);
      component_area_km2.push_back(area);
      component_n_cells.push_back(counts[static_cast<std::size_t>(component)]);
      component_first_cell.push_back(first_cells[static_cast<std::size_t>(component)]);
    }
    if (profile) area_seconds += elapsed_seconds(subphase_started);

    output_cells.insert(output_cells.end(), patch_cells.begin(), patch_cells.end());
    output_component_id.insert(
      output_component_id.end(), labels.begin(), labels.end()
    );
    if (output_cells.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      Rcpp::stop("Fragmentation analysis output exceeds the integer offset contract.");
    }
    origin_offsets[static_cast<std::size_t>(origin) + 1L] =
      static_cast<int>(output_cells.size());
    component_count_by_origin.push_back(component_count);
  }

  return Rcpp::List::create(
    Rcpp::Named("origin_patch_ids") = Rcpp::clone(changed_patch_ids),
    Rcpp::Named("origin_offsets") = Rcpp::wrap(origin_offsets),
    Rcpp::Named("cells") = Rcpp::wrap(output_cells),
    Rcpp::Named("component_id_by_cell") = Rcpp::wrap(output_component_id),
    Rcpp::Named("component_origin_patch_id") = Rcpp::wrap(component_origin_patch_id),
    Rcpp::Named("component_id") = Rcpp::wrap(component_id),
    Rcpp::Named("component_area_km2") = Rcpp::wrap(component_area_km2),
    Rcpp::Named("component_n_cells") = Rcpp::wrap(component_n_cells),
    Rcpp::Named("component_first_cell") = Rcpp::wrap(component_first_cell),
    Rcpp::Named("component_count_by_origin") = Rcpp::wrap(component_count_by_origin),
    Rcpp::Named("timing_seconds") = Rcpp::NumericVector::create(
      Rcpp::Named("patch_cell_lookup_seconds") = lookup_seconds,
      Rcpp::Named("component_label_seconds") = label_seconds,
      Rcpp::Named("fragment_area_seconds") = area_seconds
    )
  );
}


// ---- Dense Stage 6 and compact Stage 7.2 fragmentation exports -------------

// Stage 6 retains the original dense-state contract and cross-checks the live
// compact index against the dense patch-by-cell state.
// [[Rcpp::export]]
Rcpp::List stage6_analyze_fragmentation_components_cpp(
    const Rcpp::IntegerVector& index_pid,
    const Rcpp::IntegerVector& index_cell,
    const Rcpp::IntegerVector& changed_patch_ids,
    const Rcpp::IntegerVector& current_patch_id_by_cell,
    const Rcpp::NumericVector& cell_area_by_cell,
    const Rcpp::IntegerVector& rook_start,
    const Rcpp::IntegerVector& rook_end,
    const Rcpp::IntegerVector& rook_to,
    int n_cells,
    bool profile = false) {
  return analyze_fragmentation_components_impl(
    index_pid, index_cell, changed_patch_ids, current_patch_id_by_cell,
    cell_area_by_cell, rook_start, rook_end, rook_to, n_cells, profile, false
  );
}


// Stage 7.2 uses its already validated compact index as the live cell state and
// avoids allocating or scanning the dense patch-by-cell representation.
// [[Rcpp::export]]
Rcpp::List stage72_analyze_fragmentation_components_cpp(
    const Rcpp::IntegerVector& index_pid,
    const Rcpp::IntegerVector& index_cell,
    const Rcpp::IntegerVector& changed_patch_ids,
    const Rcpp::NumericVector& cell_area_by_cell,
    const Rcpp::IntegerVector& rook_start,
    const Rcpp::IntegerVector& rook_end,
    const Rcpp::IntegerVector& rook_to,
    int n_cells,
    bool profile = false) {
  return analyze_fragmentation_components_impl(
    index_pid, index_cell, changed_patch_ids, Rcpp::IntegerVector(0),
    cell_area_by_cell, rook_start, rook_end, rook_to, n_cells, profile, true
  );
}


// ---- Frontier scoring -------------------------------------------------------

// Dense reference export: scan every active species for each frontier cell and
// accumulate log-sum-exp terms in caller-provided canonical species order.
// [[Rcpp::export]]
Rcpp::List stage6_score_frontier_cpp(
    const Rcpp::IntegerVector& frontier_cells,
    const Rcpp::List& patch_ids_by_species,
    const Rcpp::List& scores_by_species,
    const bool return_diagnostics = false) {
  const R_xlen_t n_frontier = frontier_cells.size();
  const R_xlen_t n_species = patch_ids_by_species.size();
  if (scores_by_species.size() != n_species) {
    Rcpp::stop("Stage 6 frontier kernel received mismatched species lists.");
  }

  Rcpp::NumericVector max_log_term(n_frontier, R_NegInf);
  Rcpp::NumericVector sum_exp_terms(n_frontier, 0.0);
  std::vector<unsigned char> has_score(static_cast<std::size_t>(n_frontier), 0);
  Rcpp::Function base_exp = Rcpp::Environment::base_env()["exp"];
  double finite_contributions = 0;
  double scan_seconds = 0;
  double accumulation_seconds = 0;

  for (R_xlen_t species_index = 0; species_index < n_species; ++species_index) {
    Rcpp::IntegerVector patch_ids = patch_ids_by_species[species_index];
    Rcpp::NumericVector score_by_patch_id = scores_by_species[species_index];
    const R_xlen_t n_cells = patch_ids.size();
    const R_xlen_t n_scores = score_by_patch_id.size();

    std::vector<R_xlen_t> valid_positions;
    std::vector<double> valid_scores;
    valid_positions.reserve(static_cast<std::size_t>(n_frontier));
    valid_scores.reserve(static_cast<std::size_t>(n_frontier));
    const auto scan_started = std::chrono::steady_clock::now();
    for (R_xlen_t frontier_index = 0; frontier_index < n_frontier; ++frontier_index) {
      const int cell_id = frontier_cells[frontier_index];
      if (cell_id == NA_INTEGER || cell_id < 1 || cell_id > n_cells) {
        Rcpp::stop("Stage 6 frontier kernel received an invalid frontier cell ID.");
      }

      const int patch_id = patch_ids[cell_id - 1];
      if (patch_id == NA_INTEGER || patch_id < 1 || patch_id > n_scores) {
        continue;
      }

      const double log_score = score_by_patch_id[patch_id - 1];
      if (!std::isfinite(log_score)) {
        continue;
      }
      valid_positions.push_back(frontier_index);
      valid_scores.push_back(log_score);
    }
    scan_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - scan_started
    ).count();
    finite_contributions += static_cast<double>(valid_positions.size());

    const auto accumulation_started = std::chrono::steady_clock::now();
    std::vector<R_xlen_t> old_max_positions;
    std::vector<double> old_max_differences;
    std::vector<R_xlen_t> new_max_positions;
    std::vector<double> new_max_values;
    std::vector<double> new_max_differences;
    old_max_positions.reserve(valid_positions.size());
    old_max_differences.reserve(valid_positions.size());
    new_max_positions.reserve(valid_positions.size());
    new_max_values.reserve(valid_positions.size());
    new_max_differences.reserve(valid_positions.size());

    for (std::size_t i = 0; i < valid_positions.size(); ++i) {
      const R_xlen_t position = valid_positions[i];
      const double log_score = valid_scores[i];
      if (!has_score[static_cast<std::size_t>(position)]) {
        max_log_term[position] = log_score;
        sum_exp_terms[position] = 1.0;
        has_score[static_cast<std::size_t>(position)] = 1;
      } else if (log_score <= max_log_term[position]) {
        old_max_positions.push_back(position);
        old_max_differences.push_back(log_score - max_log_term[position]);
      } else {
        new_max_positions.push_back(position);
        new_max_values.push_back(log_score);
        new_max_differences.push_back(max_log_term[position] - log_score);
      }
    }

    if (!old_max_positions.empty()) {
      Rcpp::NumericVector differences = Rcpp::wrap(old_max_differences);
      Rcpp::NumericVector exponentials = base_exp(differences);
      for (std::size_t i = 0; i < old_max_positions.size(); ++i) {
        const R_xlen_t position = old_max_positions[i];
        sum_exp_terms[position] = sum_exp_terms[position] + exponentials[i];
      }
    }
    if (!new_max_positions.empty()) {
      Rcpp::NumericVector differences = Rcpp::wrap(new_max_differences);
      Rcpp::NumericVector exponentials = base_exp(differences);
      for (std::size_t i = 0; i < new_max_positions.size(); ++i) {
        const R_xlen_t position = new_max_positions[i];
        volatile double scaled_sum = sum_exp_terms[position] * exponentials[i];
        sum_exp_terms[position] = scaled_sum + 1.0;
        max_log_term[position] = new_max_values[i];
      }
    }
    accumulation_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - accumulation_started
    ).count();
  }
  const auto output_started = std::chrono::steady_clock::now();
  Rcpp::LogicalVector has_score_output(n_frontier);
  double scored_frontier_cells = 0;
  for (R_xlen_t i = 0; i < n_frontier; ++i) {
    has_score_output[i] = has_score[static_cast<std::size_t>(i)] != 0;
    if (has_score[static_cast<std::size_t>(i)] != 0) {
      scored_frontier_cells += 1;
    }
  }
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("max_log_term") = max_log_term,
    Rcpp::Named("sum_exp_terms") = sum_exp_terms,
    Rcpp::Named("has_score") = has_score_output
  );
  if (return_diagnostics) {
    const double output_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - output_started
    ).count();
    output["diagnostic_counts"] = Rcpp::NumericVector::create(
      Rcpp::Named("score_calls") = 1,
      Rcpp::Named("frontier_cells") = static_cast<double>(n_frontier),
      Rcpp::Named("active_species") = static_cast<double>(n_species),
      Rcpp::Named("species_frontier_probes") =
        static_cast<double>(n_frontier) * static_cast<double>(n_species),
      Rcpp::Named("sparse_membership_probes") =
        static_cast<double>(n_frontier) * static_cast<double>(n_species),
      Rcpp::Named("dense_probes_avoided") = 0,
      Rcpp::Named("finite_contributions") = finite_contributions,
      Rcpp::Named("scored_frontier_cells") = scored_frontier_cells
    );
    output["diagnostic_timing_seconds"] = Rcpp::NumericVector::create(
      Rcpp::Named("kernel_scan_seconds") = scan_seconds,
      Rcpp::Named("kernel_accumulation_seconds") = accumulation_seconds,
      Rcpp::Named("kernel_output_seconds") = output_seconds
    );
  }
  return output;
}


// Sparse export: visit only cell-species memberships, while preserving the same
// deterministic accumulation order and output contract as dense scoring.
// [[Rcpp::export]]
Rcpp::List stage6_score_frontier_sparse_cpp(
    const Rcpp::IntegerVector& frontier_cells,
    const Rcpp::IntegerVector& cell_offsets,
    SEXP species_ids_input,
    const int canonical_species_count,
    const Rcpp::IntegerVector& active_species_ids,
    const Rcpp::List& patch_ids_by_species,
    const Rcpp::List& scores_by_species,
    const bool return_diagnostics = false) {
  const R_xlen_t n_frontier = frontier_cells.size();
  const R_xlen_t n_active_species = patch_ids_by_species.size();
  if (scores_by_species.size() != n_active_species ||
      active_species_ids.size() != n_active_species) {
    Rcpp::stop("Sparse frontier kernel received mismatched active-species lists.");
  }
  if (canonical_species_count < 1 ||
      cell_offsets.size() < 2 ||
      cell_offsets[0] != 0) {
    Rcpp::stop("Sparse frontier kernel received an invalid cell-species index.");
  }
  const int n_cells = static_cast<int>(cell_offsets.size() - 1);
  const bool raw_species_ids = TYPEOF(species_ids_input) == RAWSXP;
  const bool integer_species_ids = TYPEOF(species_ids_input) == INTSXP;
  if (!raw_species_ids && !integer_species_ids) {
    Rcpp::stop("Sparse frontier species IDs must be raw or integer.");
  }
  const R_xlen_t membership_count = Rf_xlength(species_ids_input);
  if (cell_offsets[n_cells] != membership_count) {
    Rcpp::stop("Sparse frontier offsets do not match membership storage.");
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    if (cell_offsets[cell] < 0 ||
        cell_offsets[cell + 1] < cell_offsets[cell]) {
      Rcpp::stop("Sparse frontier offsets must be nondecreasing.");
    }
  }

  std::vector<int> active_index_by_canonical(
    static_cast<std::size_t>(canonical_species_count), -1
  );
  for (R_xlen_t active_index = 0;
       active_index < n_active_species;
       ++active_index) {
    const int species_id = active_species_ids[active_index];
    if (species_id == NA_INTEGER ||
        species_id < 1 ||
        species_id > canonical_species_count ||
        active_index_by_canonical[static_cast<std::size_t>(species_id - 1)] >= 0) {
      Rcpp::stop("Sparse frontier active species IDs are invalid or duplicated.");
    }
    active_index_by_canonical[static_cast<std::size_t>(species_id - 1)] =
      static_cast<int>(active_index);
  }

  std::vector<std::vector<R_xlen_t> > positions_by_species(
    static_cast<std::size_t>(n_active_species)
  );
  double sparse_membership_probes = 0;
  const auto scan_started = std::chrono::steady_clock::now();
  if (raw_species_ids) {
    const Rcpp::RawVector species_ids(species_ids_input);
    for (R_xlen_t frontier_index = 0; frontier_index < n_frontier; ++frontier_index) {
      const int cell_id = frontier_cells[frontier_index];
      if (cell_id == NA_INTEGER || cell_id < 1 || cell_id > n_cells) {
        Rcpp::stop("Sparse frontier kernel received an invalid frontier cell ID.");
      }
      const int first = cell_offsets[cell_id - 1];
      const int last = cell_offsets[cell_id];
      for (int position = first; position < last; ++position) {
        const int canonical_index = static_cast<int>(species_ids[position]);
        if (canonical_index < 0 || canonical_index >= canonical_species_count) {
          Rcpp::stop("Sparse frontier raw species ID is outside the canonical range.");
        }
        const int active_index =
          active_index_by_canonical[static_cast<std::size_t>(canonical_index)];
        if (active_index >= 0) {
          positions_by_species[static_cast<std::size_t>(active_index)].push_back(
            frontier_index
          );
          sparse_membership_probes += 1;
        }
      }
    }
  } else {
    const Rcpp::IntegerVector species_ids(species_ids_input);
    for (R_xlen_t frontier_index = 0; frontier_index < n_frontier; ++frontier_index) {
      const int cell_id = frontier_cells[frontier_index];
      if (cell_id == NA_INTEGER || cell_id < 1 || cell_id > n_cells) {
        Rcpp::stop("Sparse frontier kernel received an invalid frontier cell ID.");
      }
      const int first = cell_offsets[cell_id - 1];
      const int last = cell_offsets[cell_id];
      for (int position = first; position < last; ++position) {
        const int species_id = species_ids[position];
        if (species_id == NA_INTEGER ||
            species_id < 1 ||
            species_id > canonical_species_count) {
          Rcpp::stop("Sparse frontier integer species ID is outside the canonical range.");
        }
        const int active_index =
          active_index_by_canonical[static_cast<std::size_t>(species_id - 1)];
        if (active_index >= 0) {
          positions_by_species[static_cast<std::size_t>(active_index)].push_back(
            frontier_index
          );
          sparse_membership_probes += 1;
        }
      }
    }
  }
  double scan_seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - scan_started
  ).count();

  Rcpp::NumericVector max_log_term(n_frontier, R_NegInf);
  Rcpp::NumericVector sum_exp_terms(n_frontier, 0.0);
  std::vector<unsigned char> has_score(static_cast<std::size_t>(n_frontier), 0);
  Rcpp::Function base_exp = Rcpp::Environment::base_env()["exp"];
  double finite_contributions = 0;
  double accumulation_seconds = 0;

  for (R_xlen_t species_index = 0;
       species_index < n_active_species;
       ++species_index) {
    Rcpp::IntegerVector patch_ids = patch_ids_by_species[species_index];
    Rcpp::NumericVector score_by_patch_id = scores_by_species[species_index];
    const R_xlen_t n_cells_for_species = patch_ids.size();
    const R_xlen_t n_scores = score_by_patch_id.size();
    const std::vector<R_xlen_t>& candidate_positions =
      positions_by_species[static_cast<std::size_t>(species_index)];

    std::vector<R_xlen_t> valid_positions;
    std::vector<double> valid_scores;
    valid_positions.reserve(candidate_positions.size());
    valid_scores.reserve(candidate_positions.size());
    for (R_xlen_t frontier_index : candidate_positions) {
      const int cell_id = frontier_cells[frontier_index];
      if (cell_id > n_cells_for_species) {
        Rcpp::stop("Sparse frontier species vector has the wrong cell count.");
      }
      const int patch_id = patch_ids[cell_id - 1];
      if (patch_id == NA_INTEGER || patch_id < 1 || patch_id > n_scores) continue;
      const double log_score = score_by_patch_id[patch_id - 1];
      if (!std::isfinite(log_score)) continue;
      valid_positions.push_back(frontier_index);
      valid_scores.push_back(log_score);
    }
    finite_contributions += static_cast<double>(valid_positions.size());

    const auto accumulation_started = std::chrono::steady_clock::now();
    std::vector<R_xlen_t> old_max_positions;
    std::vector<double> old_max_differences;
    std::vector<R_xlen_t> new_max_positions;
    std::vector<double> new_max_values;
    std::vector<double> new_max_differences;
    old_max_positions.reserve(valid_positions.size());
    old_max_differences.reserve(valid_positions.size());
    new_max_positions.reserve(valid_positions.size());
    new_max_values.reserve(valid_positions.size());
    new_max_differences.reserve(valid_positions.size());

    for (std::size_t i = 0; i < valid_positions.size(); ++i) {
      const R_xlen_t position = valid_positions[i];
      const double log_score = valid_scores[i];
      if (!has_score[static_cast<std::size_t>(position)]) {
        max_log_term[position] = log_score;
        sum_exp_terms[position] = 1.0;
        has_score[static_cast<std::size_t>(position)] = 1;
      } else if (log_score <= max_log_term[position]) {
        old_max_positions.push_back(position);
        old_max_differences.push_back(log_score - max_log_term[position]);
      } else {
        new_max_positions.push_back(position);
        new_max_values.push_back(log_score);
        new_max_differences.push_back(max_log_term[position] - log_score);
      }
    }
    if (!old_max_positions.empty()) {
      Rcpp::NumericVector differences = Rcpp::wrap(old_max_differences);
      Rcpp::NumericVector exponentials = base_exp(differences);
      for (std::size_t i = 0; i < old_max_positions.size(); ++i) {
        const R_xlen_t position = old_max_positions[i];
        sum_exp_terms[position] += exponentials[i];
      }
    }
    if (!new_max_positions.empty()) {
      Rcpp::NumericVector differences = Rcpp::wrap(new_max_differences);
      Rcpp::NumericVector exponentials = base_exp(differences);
      for (std::size_t i = 0; i < new_max_positions.size(); ++i) {
        const R_xlen_t position = new_max_positions[i];
        volatile double scaled_sum = sum_exp_terms[position] * exponentials[i];
        sum_exp_terms[position] = scaled_sum + 1.0;
        max_log_term[position] = new_max_values[i];
      }
    }
    accumulation_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - accumulation_started
    ).count();
  }

  const auto output_started = std::chrono::steady_clock::now();
  Rcpp::LogicalVector has_score_output(n_frontier);
  double scored_frontier_cells = 0;
  for (R_xlen_t i = 0; i < n_frontier; ++i) {
    has_score_output[i] = has_score[static_cast<std::size_t>(i)] != 0;
    if (has_score[static_cast<std::size_t>(i)] != 0) scored_frontier_cells += 1;
  }
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("max_log_term") = max_log_term,
    Rcpp::Named("sum_exp_terms") = sum_exp_terms,
    Rcpp::Named("has_score") = has_score_output
  );
  if (return_diagnostics) {
    const double output_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - output_started
    ).count();
    const double dense_probes =
      static_cast<double>(n_frontier) * static_cast<double>(n_active_species);
    output["diagnostic_counts"] = Rcpp::NumericVector::create(
      Rcpp::Named("score_calls") = 1,
      Rcpp::Named("frontier_cells") = static_cast<double>(n_frontier),
      Rcpp::Named("active_species") = static_cast<double>(n_active_species),
      Rcpp::Named("species_frontier_probes") = dense_probes,
      Rcpp::Named("sparse_membership_probes") = sparse_membership_probes,
      Rcpp::Named("dense_probes_avoided") =
        std::max(0.0, dense_probes - sparse_membership_probes),
      Rcpp::Named("finite_contributions") = finite_contributions,
      Rcpp::Named("scored_frontier_cells") = scored_frontier_cells
    );
    output["diagnostic_timing_seconds"] = Rcpp::NumericVector::create(
      Rcpp::Named("kernel_scan_seconds") = scan_seconds,
      Rcpp::Named("kernel_accumulation_seconds") = accumulation_seconds,
      Rcpp::Named("kernel_output_seconds") = output_seconds
    );
  }
  return output;
}


// ---- Incremental compact-index publication ---------------------------------

// Splice sorted replacement segments into the patch-sorted compact index. The
// final vectors own their storage and preserve unique one-based cell IDs.
// [[Rcpp::export]]
Rcpp::List stage6_splice_patch_index_cpp(
    const Rcpp::IntegerVector& old_pid,
    const Rcpp::IntegerVector& old_cell,
    const Rcpp::IntegerVector& expected_patch_ids_input,
    const Rcpp::IntegerVector& pre_update_patch_ids_input,
    const Rcpp::IntegerVector& replacement_origin_ids_input,
    const Rcpp::IntegerVector& replacement_origin,
    const Rcpp::IntegerVector& replacement_pid,
    const Rcpp::IntegerVector& replacement_cell) {
  const R_xlen_t old_size = old_pid.size();
  if (old_cell.size() != old_size) {
    Rcpp::stop("patch_index pid and cell vectors must have equal lengths.");
  }
  if (replacement_origin.size() != replacement_pid.size() ||
      replacement_origin.size() != replacement_cell.size()) {
    Rcpp::stop("Compact-index replacement vectors must have equal lengths.");
  }

  std::unordered_set<int> old_cells_seen;
  old_cells_seen.reserve(static_cast<std::size_t>(old_size));
  std::vector<int> existing_patch_ids;
  std::unordered_map<int, std::pair<R_xlen_t, R_xlen_t> > bounds;
  int previous_pid = 0;
  for (R_xlen_t i = 0; i < old_size; ++i) {
    const int pid = old_pid[i];
    const int cell = old_cell[i];
    if (pid == NA_INTEGER || cell == NA_INTEGER || pid < 1 || cell < 1 ||
        (i > 0 && pid < previous_pid)) {
      Rcpp::stop("patch_index is not a valid sorted patch-to-cell index.");
    }
    if (!old_cells_seen.insert(cell).second) {
      Rcpp::stop("patch_index contains duplicated cell IDs.");
    }
    if (i == 0 || pid != previous_pid) {
      existing_patch_ids.push_back(pid);
      bounds[pid] = std::make_pair(i, i);
    } else {
      bounds[pid].second = i;
    }
    previous_pid = pid;
  }

  std::vector<int> expected_patch_ids = sorted_unique(
    expected_patch_ids_input, "expected_patch_ids"
  );
  std::vector<int> pre_update_patch_ids = sorted_unique(
    pre_update_patch_ids_input, "pre_update_patch_ids"
  );

  for (int pid : pre_update_patch_ids) {
    if (!contains_sorted(existing_patch_ids, pid)) {
      Rcpp::stop("Live pre-update patch IDs are missing from the compact index: %d.", pid);
    }
  }

  std::vector<int> replacement_origins;
  replacement_origins.reserve(replacement_origin_ids_input.size());
  std::unordered_set<int> replacement_origin_seen;
  for (int origin : replacement_origin_ids_input) {
    if (origin == NA_INTEGER || origin < 1 ||
        !replacement_origin_seen.insert(origin).second) {
      Rcpp::stop("replacement_rows_by_origin must have unique positive-integer names.");
    }
    replacement_origins.push_back(origin);
  }
  std::unordered_map<int, std::vector<R_xlen_t> > replacement_rows;
  for (R_xlen_t i = 0; i < replacement_origin.size(); ++i) {
    const int origin = replacement_origin[i];
    const int pid = replacement_pid[i];
    const int cell = replacement_cell[i];
    if (origin == NA_INTEGER || origin < 1 || pid == NA_INTEGER || pid < 1 ||
        cell == NA_INTEGER || cell < 1) {
      Rcpp::stop("Compact-index replacements must contain positive integers.");
    }
    if (replacement_origin_seen.find(origin) == replacement_origin_seen.end()) {
      Rcpp::stop("A compact-index replacement row has an unknown origin patch ID.");
    }
    replacement_rows[origin].push_back(i);
  }
  std::sort(replacement_origins.begin(), replacement_origins.end());
  replacement_origins.erase(
    std::unique(replacement_origins.begin(), replacement_origins.end()),
    replacement_origins.end()
  );
  for (int origin : replacement_origins) {
    if (!contains_sorted(existing_patch_ids, origin)) {
      Rcpp::stop("Compact-index replacements reference missing origin patch IDs: %d.", origin);
    }
    if (!contains_sorted(pre_update_patch_ids, origin)) {
      Rcpp::stop("Compact-index replacements must originate from live pre-update patches.");
    }
  }

  std::vector<int> affected_ids = replacement_origins;
  for (int pid : existing_patch_ids) {
    if (!contains_sorted(pre_update_patch_ids, pid)) affected_ids.push_back(pid);
  }
  for (int pid : pre_update_patch_ids) {
    if (!contains_sorted(expected_patch_ids, pid)) affected_ids.push_back(pid);
  }
  std::sort(affected_ids.begin(), affected_ids.end());
  affected_ids.erase(std::unique(affected_ids.begin(), affected_ids.end()), affected_ids.end());

  std::unordered_set<int> affected_set(affected_ids.begin(), affected_ids.end());
  std::unordered_map<int, std::vector<std::pair<int, int> > > retained_by_origin;
  std::vector<std::pair<int, int> > appended;
  std::vector<int> removed_cells;
  std::vector<int> relabel_cells;
  std::vector<int> relabel_patch_ids;
  R_xlen_t final_size = 0;
  old_cells_seen.clear();
  old_cells_seen.rehash(0);

  for (int origin : existing_patch_ids) {
    const std::pair<R_xlen_t, R_xlen_t> segment = bounds[origin];
    if (affected_set.find(origin) == affected_set.end()) {
      final_size += segment.second - segment.first + 1;
      continue;
    }

    std::unordered_set<int> old_origin_cells;
    old_origin_cells.reserve(static_cast<std::size_t>(segment.second - segment.first + 1));
    for (R_xlen_t i = segment.first; i <= segment.second; ++i) {
      old_origin_cells.insert(old_cell[i]);
    }
    std::unordered_set<int> replacement_cells;
    std::vector<std::pair<int, int> > retained_origin_rows;

    const auto replacement_it = replacement_rows.find(origin);
    if (replacement_it != replacement_rows.end()) {
      for (R_xlen_t row : replacement_it->second) {
        const int pid = replacement_pid[row];
        const int cell = replacement_cell[row];
        if (!contains_sorted(expected_patch_ids, pid)) continue;
        if (old_origin_cells.find(cell) == old_origin_cells.end() ||
            !replacement_cells.insert(cell).second) {
          Rcpp::stop("Invalid compact-index replacement for origin patch %d.", origin);
        }
        if (pid == origin) {
          retained_origin_rows.emplace_back(cell, pid);
        } else {
          appended.emplace_back(pid, cell);
        }
      }
    }

    std::sort(retained_origin_rows.begin(), retained_origin_rows.end());
    retained_by_origin[origin] = retained_origin_rows;
    final_size += static_cast<R_xlen_t>(retained_origin_rows.size());
    for (R_xlen_t i = segment.first; i <= segment.second; ++i) {
      if (replacement_cells.find(old_cell[i]) == replacement_cells.end()) {
        removed_cells.push_back(old_cell[i]);
      }
    }
  }

  std::sort(appended.begin(), appended.end());
  const int max_pre_update = pre_update_patch_ids.empty() ? 0 : pre_update_patch_ids.back();
  for (const auto& row : appended) {
    if (row.first <= max_pre_update) {
      Rcpp::stop("New fragment patch IDs must exceed every live pre-update patch ID.");
    }
  }
  final_size += static_cast<R_xlen_t>(appended.size());

  Rcpp::IntegerVector final_pid = Rcpp::no_init(final_size);
  Rcpp::IntegerVector final_cell = Rcpp::no_init(final_size);
  R_xlen_t output_position = 0;
  for (int origin : existing_patch_ids) {
    const std::pair<R_xlen_t, R_xlen_t> segment = bounds[origin];
    if (affected_set.find(origin) == affected_set.end()) {
      for (R_xlen_t i = segment.first; i <= segment.second; ++i) {
        final_pid[output_position] = old_pid[i];
        final_cell[output_position] = old_cell[i];
        ++output_position;
      }
    } else {
      const std::vector<std::pair<int, int> >& retained = retained_by_origin[origin];
      for (const auto& row : retained) {
        final_pid[output_position] = row.second;
        final_cell[output_position] = row.first;
        ++output_position;
      }
    }
  }
  for (const auto& row : appended) {
    final_pid[output_position] = row.first;
    final_cell[output_position] = row.second;
    ++output_position;
  }
  if (output_position != final_size) {
    Rcpp::stop("Compact-index kernel wrote an unexpected number of rows.");
  }

  std::vector<int> represented;
  int last_pid = 0;
  std::unordered_set<int> final_cells_seen;
  final_cells_seen.reserve(static_cast<std::size_t>(final_size));
  for (R_xlen_t i = 0; i < final_size; ++i) {
    if (i > 0 && final_pid[i] < final_pid[i - 1]) {
      Rcpp::stop("Incremental patch index is not sorted by patch ID.");
    }
    if (i == 0 || final_pid[i] != last_pid) represented.push_back(final_pid[i]);
    if (!final_cells_seen.insert(final_cell[i]).second) {
      Rcpp::stop("Incremental patch index contains duplicated cell IDs.");
    }
    last_pid = final_pid[i];
  }
  if (represented != expected_patch_ids) {
    Rcpp::stop("Incremental patch index does not match the expected surviving state.");
  }

  // Preserve the R reference's replacement-list and row order in relabel output.
  for (R_xlen_t row = 0; row < replacement_origin.size(); ++row) {
    const int pid = replacement_pid[row];
    if (contains_sorted(expected_patch_ids, pid) && pid != replacement_origin[row]) {
      relabel_cells.push_back(replacement_cell[row]);
      relabel_patch_ids.push_back(pid);
    }
  }

  std::sort(removed_cells.begin(), removed_cells.end());
  removed_cells.erase(std::unique(removed_cells.begin(), removed_cells.end()), removed_cells.end());

  return Rcpp::List::create(
    Rcpp::Named("pid") = final_pid,
    Rcpp::Named("cell") = final_cell,
    Rcpp::Named("removed_cells") = Rcpp::wrap(removed_cells),
    Rcpp::Named("relabel_cells") = Rcpp::wrap(relabel_cells),
    Rcpp::Named("relabel_patch_ids") = Rcpp::wrap(relabel_patch_ids)
  );
}
