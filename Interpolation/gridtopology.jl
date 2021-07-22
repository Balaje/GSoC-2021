""" Functions to implement GridTopology """

# Restricted Triangulation
function GridTopology(trian::RestrictedTriangulation)
  parent_trian = trian.parent_trian
  parent_cell_to_vertices, parent_vertex_to_node, = _generate_cell_to_vertices_from_grid(parent_trian)
  cell_to_vertices = Table(lazy_map(Reindex(parent_cell_to_vertices),trian.cell_to_parent_cell))
  vertex_to_node = parent_vertex_to_node

  @notimplementedif (! is_regular(parent_trian)) "Extrtacting the GridTopology form a Grid only implemented for the regular case"
  node_to_coords = get_node_coordinates(trian)
  if vertex_to_node == 1:num_nodes(trian)
    vertex_to_coords = node_to_coords
  else
    vertex_to_coords = node_to_coords[vertex_to_node]
  end
  cell_to_type = collect(get_cell_type(trian))
  polytopes = map(get_polytope, get_reffes(trian))

  UnstructuredGridTopology(vertex_to_coords, cell_to_vertices, cell_to_type, polytopes, OrientationStyle(parent_trian))
end

# Boundary Triangulation
GridTopology(trian::BoundaryTriangulation) = GridTopology(trian.face_trian)
