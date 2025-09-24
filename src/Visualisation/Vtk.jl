"""
To visualise on the sphere (and other mapped domains), we trick vtk to plot
cellvalues on a visualisation mesh that is different to the one used for
evaluation of cell/node data.
This is achieved by evaluating a cellwise geo_map in _vtkpoints.
If no geo_map is provided, then evaluate the regular visualisation grid.
Need to dispatch through writevtk/createvtk to create_vtk_file
"""

function Gridap.Visualization.writevtk(args...;geo_map=nothing,
  compress=false,append=true,ascii=false,vtkversion=:default,kwargs...)
  map(Gridap.Visualization.visualization_data(args...;kwargs...)) do visdata
    write_vtk_file(
      visdata.grid,visdata.filebase,celldata=visdata.celldata,nodaldata=visdata.nodaldata,
      geo_map=geo_map,
      compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
    )
  end
end

"""
"""
function Gridap.Visualization.createvtk(args...;geo_map=nothing,
  compress=false,append=true,ascii=false,vtkversion=:default,kwargs...)
  v = Gridap.Visualization.visualization_data(args...;kwargs...)
  @notimplementedif length(v) != 1
  visdata = first(v)
  Gridap.Visualization.create_vtk_file(
    visdata.grid,visdata.filebase,celldata=visdata.celldata,nodaldata=visdata.nodaldata,
    geo_map=geo_map,
    compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
  )
end

function Gridap.Visualization.write_vtk_file(
  trian::Grid, filebase; celldata=Dict(), nodaldata=Dict(),
  geo_map=nothing,
  compress=false, append=true, ascii=false, vtkversion=:default
)
  vtkfile = Gridap.Visualization.create_vtk_file(
    trian, filebase, celldata=celldata, nodaldata=nodaldata,
    geo_map=geo_map,
    compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
  )
  outfiles = Gridap.Visualization.vtk_save(vtkfile)
end

function Gridap.Visualization.create_vtk_file(
  trian::Grid, filebase; celldata=Dict(), nodaldata=Dict(),
  geo_map=nothing,
  compress=false, append=true, ascii=false, vtkversion=:default
)
  println("my vis")

  # compute the regular visualisation points
  points = Gridap.Visualization._vtkpoints(trian)

  ## if geo_map provided, map the points to ambient space
  if geo_map != nothing
    points = mapped_vtkpoints(trian,geo_map)
  end


  cells = Gridap.Visualization._vtkcells(trian)
  vtkfile = Gridap.Visualization.vtk_grid(
    filebase, points, cells,
    compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
  )

  if num_cells(trian)>0
    for (k,v) in celldata
      component_names = Gridap.Visualization._data_component_names(v)
      Gridap.Visualization.vtk_cell_data(vtkfile, Gridap.Visualization._prepare_data(v), k; component_names)
    end
    for (k,v) in nodaldata
      component_names = Gridap.Visualization._data_component_names(v)
      Gridap.Visualization.vtk_point_data(vtkfile, Gridap.Visualization._prepare_data(v), k; component_names)
    end
  end

  return vtkfile
end


function mapped_vtkpoints(trian,geo_map::AbstractArray)
  println("mapped vkpoints")

  # apply the geo_map to cell_coords on trian, then convert to nodes
  cellx = get_cell_coordinates(trian)

  display(cellx)
  # println(length(cellx))
  # geo_map = geo_map[1:length(cellx)]

  cellx_mapped = lazy_map(evaluate,geo_map,cellx)

  # println(cellx_mapped)
  x_mapped, cell_to_offset = Gridap.Visualization._prepare_node_to_coords(cellx_mapped)

  T = eltype(x_mapped)
  D = num_components(T)

  xflat = collect(x_mapped)
  reshape(reinterpret(Float64,xflat),(D,length(x_mapped)))
end
