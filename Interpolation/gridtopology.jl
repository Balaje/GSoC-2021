"""
Function to implement GridTopology for BoundaryTriangulation
"""

function GridTopology(grid::BoundaryTriangulation)
  UnstructuredGridTopology(grid)
end

function UnstructuredGridTopology(grid::BoundaryTriangulation)
  cell_to_vertices, vertex_to_node, = _generate_cell_to_vertices_from_grid(grid)
  _generate_grid_topology_from_grid(grid,cell_to_vertices,vertex_to_node)
end

function UnstructuredGridTopology(grid::BoundaryTriangulation, cell_to_vertices::Table, vertex_to_node::AbstractVector)
  _generate_grid_topology_from_grid(grid,cell_to_vertices,vertex_to_node)
end

function _generate_grid_topology_from_grid(grid::BoundaryTriangulation, cell_to_vertices, vertex_to_node)

  node_to_coords = get_node_coordinates(grid)
  if vertex_to_node == 1:num_nodes(grid)
    vertex_to_coords = node_to_coords
  else
    vertex_to_coords = node_to_coords[vertex_to_node]
  end

  cell_to_type = get_cell_type(grid)
  polytopes = map(get_polytope, get_reffes(grid))

  Geometry.UnstructuredGridTopology(
    vertex_to_coords,
    cell_to_vertices,
    collect(cell_to_type),
    polytopes,
    NonOriented())
end

function _generate_cell_to_vertices_from_grid(grid::BoundaryTriangulation)
  if is_first_order(grid)
    cell_to_vertices = Table(get_cell_node_ids(grid))
    vertex_to_node = collect(1:num_nodes(grid))
    node_to_vertex = vertex_to_node
  else
    cell_to_nodes = get_cell_node_ids(grid)
    cell_to_cell_type = get_cell_type(grid)
    reffes = get_reffes(grid)
    cell_type_to_lvertex_to_lnode = map(get_vertex_node, reffes)
    cell_to_vertices, vertex_to_node, node_to_vertex = Geometry._generate_cell_to_vertices(
      cell_to_nodes,
      cell_to_cell_type,
      cell_type_to_lvertex_to_lnode,
      num_nodes(grid))
  end
  (cell_to_vertices, vertex_to_node, node_to_vertex)
end
