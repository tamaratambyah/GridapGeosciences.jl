function Gridap.Visualization.write_vtk_file(
  parts::AbstractArray,
  grid::AbstractArray{<:Gridap.Geometry.Grid}, filebase; celldata, nodaldata,
  geo_map,
  compress=false,append=true,ascii=false,vtkversion=:default
  )
  pvtk = Gridap.Visualization.create_vtk_file(
    parts,grid,filebase;celldata=celldata,nodaldata=nodaldata,
    geo_map=geo_map,
    compress=compress,append=append,ascii=ascii,vtkversion=vtkversion
  )
  map(Gridap.Visualization.vtk_save,pvtk)
end

function Gridap.Visualization.create_vtk_file(
  parts::AbstractArray,
  grid::AbstractArray{<:Gridap.Geometry.Grid},
  filebase;
  celldata, nodaldata,
  geo_map = map(p -> nothing, parts),
  compress=false,append=true,ascii=false,vtkversion=:default
)
  nparts = length(parts)
  map(parts,grid,celldata,nodaldata,geo_map) do part,g,c,n,gm
    Gridap.Visualization.create_pvtk_file(
      g,filebase;
      part=part,nparts=nparts,
      celldata=c,nodaldata=n,
      geo_map=gm,
      compress=compress,append=append,ascii=ascii,vtkversion=vtkversion
    )
  end
end

function Gridap.Visualization.writevtk(
  arg::GridapDistributed.DistributedModelOrTriangulation,args...;
  geo_map= map(p -> nothing, local_views(arg)),
  compress=false,append=true,ascii=false,vtkversion=:default,kwargs...
)
  parts=get_parts(arg)
  map(GridapDistributed.visualization_data(arg,args...;kwargs...)) do visdata
    Gridap.Visualization.write_vtk_file(
      parts,visdata.grid,visdata.filebase,celldata=visdata.celldata,nodaldata=visdata.nodaldata,
      geo_map=geo_map,
      compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
    )
  end
end

function Gridap.Visualization.createvtk(
  arg::GridapDistributed.DistributedModelOrTriangulation,args...;
  geo_map=nothing,
  compress=false,append=true,ascii=false,vtkversion=:default,kwargs...
)
  v = Gridap.Visualization.visualization_data(arg,args...;kwargs...)
  parts=get_parts(arg)
  @Gridap.Helpers.notimplementedif length(v) != 1
  visdata = first(v)
  Gridap.Visualization.create_vtk_file(
    parts,visdata.grid,visdata.filebase,celldata=visdata.celldata,nodaldata=visdata.nodaldata,
    geo_map=geo_map,
    compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
  )
end

function Gridap.Visualization.create_pvtk_file(
  trian::Gridap.Geometry.Grid, filebase; part, nparts, ismain=(part==1), celldata=Dict(), nodaldata=Dict(),
  geo_map=nothing,
  compress=false, append=true, ascii=false, vtkversion=:default
)
  println("my distributed vis")


  # compute the regular visualisation points
  points = Gridap.Visualization._vtkpoints(trian)

  # println(typeof(geo_map)<:AbstractArray)
  ## if geo_map provided, map the points to ambient space
  if geo_map != nothing
    points = mapped_vtkpoints(trian,geo_map)
  end

  cells = Gridap.Visualization._vtkcells(trian)
  vtkfile = Gridap.Visualization.pvtk_grid(
    filebase, points, cells;part=part, nparts=nparts, ismain=ismain,
    compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
  )

  if num_cells(trian) > 0
    for (k, v) in celldata
      # component_names are actually always nothing as there are no field in ptvk atm
      component_names = Gridap.Visualization._data_component_names(v)
      vtkfile[k, Gridap.Visualization.VTKCellData(), component_names=component_names] = Gridap.Visualization._prepare_data(v)
    end
    for (k, v) in nodaldata
      component_names = Gridap.Visualization._data_component_names(v)
      vtkfile[k, Gridap.Visualization.VTKPointData(), component_names=component_names] = Gridap.Visualization._prepare_data(v)
    end
  end
  return vtkfile
end
