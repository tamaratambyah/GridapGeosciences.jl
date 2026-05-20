
struct AtlasGrid{Dc,
                 Dp,
                 Tp,
                 A<:AbstractVector{<:AbstractVector{<:Point{Dp,Tp}}},
                 B<:AbstractVector{<:AbstractVector{<:Integer}},
                 C<:AbstractVector{<:Integer},
                 D<:Gridap.Geometry.OrientationStyle} <: Gridap.Geometry.Grid{Dc,Dp}
  cell_coordinates::A
  cell_node_ids::B
  cell_to_chart::C
  reffe::Gridap.ReferenceFEs.LagrangianRefFE{Dc}
  orientation_style::D
  function AtlasGrid(
    cell_coordinates::AbstractVector{<:AbstractVector{<:Point{Dp,Tp}}},
    cell_node_ids::AbstractVector{<:AbstractVector{<:Integer}},
    cell_to_chart::AbstractVector{<:Integer},
    reffe::Gridap.ReferenceFEs.LagrangianRefFE{Dc},
    orientation_style::Gridap.Geometry.OrientationStyle
  ) where {Dc,Dp,Tp}
    @assert length(cell_coordinates) == length(cell_node_ids) == length(cell_to_chart) "The number of cells must be the same in cell_coordinates, cell_node_ids and cell_to_chart."
    A = typeof(cell_coordinates)
    B = typeof(cell_node_ids)
    C = typeof(cell_to_chart)
    D = typeof(orientation_style)
    new{Dc,Dp,Tp,A,B,C,D}(
      cell_coordinates,
      cell_node_ids,
      cell_to_chart,
      reffe,
      orientation_style
    )
  end
end

Gridap.Geometry.get_reffes(g::AtlasGrid) = Fill(g.reffe, length(g.cell_to_chart))
Gridap.Geometry.get_cell_type(g::AtlasGrid) = Fill(1, length(g.cell_to_chart))

function Gridap.Geometry.get_node_coordinates(g::AtlasGrid) 
    Gridap.Helpers.@notimplemented "This function is not (cannot be) implemented for AtlasGrids. Use get_cell_coordinates instead."
end 

Gridap.Geometry.get_cell_node_ids(g::AtlasGrid) = g.cell_node_ids

function Gridap.Geometry.get_cell_coordinates(g::AtlasGrid)
  g.cell_coordinates
end