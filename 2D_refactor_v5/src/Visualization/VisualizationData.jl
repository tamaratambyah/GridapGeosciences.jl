""" Visualise at quad points

Quadrature degree = order(reffe) + 1
"""

function Gridap.Visualization.visualization_data(
  ::ManifoldName,
  trian::Triangulation, filebase::AbstractString;
  order=-1, nsubcells=-1, celldata=Dict(), cellfields=Dict())

  println("my vis")
  if order == -1 && nsubcells == -1
    # Use the given cells as visualization cells
    f = (reffe) -> UnstructuredGrid(reffe)
  elseif order != -1 && nsubcells == -1
    # Use cells of given order as visualization cells
    f = (reffe) -> UnstructuredGrid(LagrangianRefFE(Float64,get_polytope(reffe),order))
  elseif order == -1 && nsubcells != -1
    # Use linear sub-cells with nsubcells per direction
    f = (reffe) -> UnstructuredGrid(compute_reference_grid(reffe,nsubcells))
  else
    @unreachable "order and nsubcells kw-arguments can not be given at the same time"
  end

  ref_grids = map(f, get_reffes(trian))
  visgrid = Gridap.Visualization.VisualizationGrid(trian,ref_grids)

  reffee = get_reffes(trian)[1]
  quad = CellQuadrature(trian,get_order(reffee)+1)
  cell_quad = quad.cell_point

  cdata = Gridap.Visualization._prepare_cdata(celldata,visgrid.sub_cell_to_cell)
  pdata = Gridap.Visualization._prepare_pdata(trian,cellfields,cell_quad)

  (Gridap.Visualization.VisualizationData(visgrid,filebase;celldata=cdata,nodaldata=pdata),)
end
