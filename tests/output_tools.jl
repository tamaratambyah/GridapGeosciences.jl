using CSV
using DataFrames
using Gridap
using GridapGeosciences
using GridapDistributed


function latlon_geo_map_func(trian::Triangulation)
  panel_ids = get_panel_ids(trian)
  latlon_geo_map_func(panel_ids)
end

function latlon_geo_map_func(panel_ids::AbstractArray{Int})
  println("latolon serial geo map")
  return lazy_map(p -> Cartesian2SphereicalMap() ∘ MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
end

function latlon_geo_map_func(trian::GridapDistributed.DistributedTriangulation)
  println("distributed latlon geo map")

  model = get_background_model(trian)
  owned_panel_ids = get_owned_panel_ids(model)
  return geo_map_func(owned_panel_ids)
end

function latlon_geo_map_func(owned_panel_ids::AbstractArray)
  println("distributed latlon geo map")

  @assert typeof(owned_panel_ids) <: DebugArray || typeof(owned_panel_ids) <: MPIArray "\n Not distributed panel ids"

  cell_geo_map = map(owned_panel_ids) do pid
    return lazy_map(p -> Cartesian2SphereicalMap() ∘MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pid)
  end
  return cell_geo_map
end

"""
save the mesh of the model
"""
function save_mesh(dir::String,model::DiscreteModel)
  trian = Triangulation(model)
  panel_ids = get_panel_ids(model)
  xs, ys = get_latlon_coords(trian,panel_ids)
  _save_mesh(dir,xs,ys)
end

function save_mesh(dir::String,dmodel::DistributedParametricDiscreteModel)
  xs_array = Float64[]
  ys_array = Float64[]

  dpanel_ids = get_owned_panel_ids(dmodel)
  dtrian = Triangulation(dmodel)

  map(local_views(dtrian),dpanel_ids) do trian,panel_ids
    xs, ys = get_latlon_coords(trian,panel_ids)
    push!(xs_array,xs...)
    push!(ys_array,ys...)
  end

  _save_mesh(dir,xs_array,ys_array)
end

function get_latlon_coords(trian::Triangulation,panel_ids::AbstractArray)
  pts = get_cell_ref_coordinates(trian)
  cmaps  = get_cell_map(trian)
  cell_geo_map = latlon_geo_map_func(panel_ids)
  ambient_cmaps = lazy_map(∘,cell_geo_map,cmaps)

  latlon_coords = lazy_map(evaluate,ambient_cmaps,pts)

  coords = reduce(vcat,collect(latlon_coords))
  xs = map(x->x[1],coords)
  ys = map(x->x[2],coords)
  return xs, ys
end


function _save_mesh(dir::String,xs_array::Vector{<:Float64},ys_array::Vector{<:Float64})
  output_file = joinpath(dir,"coords.csv")
  initialize_csv(output_file,"xs", "ys")
  append_to_csv(output_file; xs=xs_array, ys=ys_array)
end

"""
save cell fields
"""
function save_cellfields(dir::String,trian::Triangulation,t::Float64,cellfields::AbstractArray,labels::Vector{String})
  pts = get_cell_points(trian)

  cdata = map(x->reduce(vcat,collect(x(pts))), cellfields)
  _save_cellfields(dir,t,cdata,labels)
end


function save_cellfields(dir::String,trian::GridapDistributed.DistributedTriangulation,t::Float64,cellfields::AbstractArray,labels::Vector{String})
  cdata_array = [ Float64[] for i in 1:length(cellfields) ]
  pts = get_cell_points(trian)

  for (i,field) in enumerate(cellfields)
    d = field(pts)
    map(d) do cdata
      data = reduce(vcat,collect(cdata))
      push!(cdata_array[i],data...)
    end
  end

  _save_cellfields(dir,t,cdata_array,labels)
end

function _save_cellfields(dir::String,t::Float64,cdata,labels::Vector{String})
  _labs = map(x->Symbol(x),labels)
  output_file = joinpath(dir,"output_t$t.csv")
  initialize_csv(output_file,_labs...)
  fields = map((x,y)-> x=>y, _labs, cdata )
  append_to_csv(output_file;  fields...)
end


"""Write scalar diagnostics to csv.
Diagnostics should be added as kwargs: field=values"""
function write_to_csv(csv_file_path; kwargs...)
  df = DataFrame(kwargs... )
  header = names(df)
  CSV.write(csv_file_path, df, header=header)
end

"""Append line to existing csv file.
row should be given as kwargs: field=value"""
function append_to_csv(csv_file_path; kwargs...)
  df = DataFrame(kwargs )
  CSV.write(csv_file_path, df, delim=",", append=true)
end

"""Create a .cvs file header. header should be specified as args"""
function initialize_csv(csv_file_path, args...)
  header = [String(_) for _ in args]
  CSV.write(csv_file_path,[], writeheader=true, header=header)
end

"""Wrapper to get a vector from a csv field. The fieldname argument should be passed as a symbol,
i.e. :<fieldname>"""
function get_scalar_field_from_csv(csv_file_path, fieldname)
  t = CSV.read(csv_file_path, DataFrame)
  getproperty(t, fieldname)
end

function get_vector_field_from_csv(csv_file_path, fieldname)
  t = CSV.read(csv_file_path, DataFrame)
  getproperty(t, fieldname)
end
