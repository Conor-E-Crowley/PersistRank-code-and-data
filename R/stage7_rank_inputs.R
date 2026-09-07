# Canonical benchmark rank-table construction.
#
# Reads one rank-raster value vector plus the Stage 6 initial domain and cell
# areas. Returns the deterministic descending rank order used by exact-target
# reconstruction. Writes nothing. Equal ranks are resolved by raster cell ID.

build_stage7_rank_table <- function(rank_values, initially_alive,
                                    cell_area_by_cell) {
  assert(length(rank_values) == length(initially_alive),
         "Rank values and initial-domain vectors differ in length.")
  assert(length(cell_area_by_cell) == length(initially_alive),
         "Cell areas and initial-domain vectors differ in length.")
  cells <- which(as.logical(initially_alive))
  table <- data.table::data.table(
    cell = as.integer(cells),
    rank_value = as.numeric(rank_values[cells]),
    cell_area_km2 = as.numeric(cell_area_by_cell[cells])
  )
  assert(all(is.finite(table$rank_value)),
         "The rank map contains non-finite values in the initial domain.")
  assert(all(is.finite(table$cell_area_km2) & table$cell_area_km2 > 0),
         "Initial-domain cells require positive finite areas.")
  data.table::setorder(table, -rank_value, cell)
  table[, `:=`(
    cum_cells = seq_len(.N),
    cum_area_km2 = cumsum(cell_area_km2)
  )]
  table[]
}
