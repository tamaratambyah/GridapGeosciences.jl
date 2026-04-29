using Gridap
using GridapGeosciences
using Gridap.Geometry
using DrWatson
using Plots

include("convergence_tools.jl")
dir = plotsdir("plot_latlon")
!isdir(dir) && mkdir(dir)
n_ref_lvls = 2
panel_models  = get_refined_models(n_ref_lvls,true)


panel_model = panel_models[2]
panel_ids = get_panel_ids(panel_model)

cmap = get_cell_map(get_grid(panel_model))
ref_points = get_cell_ref_coordinates(panel_model)
coords = lazy_map(evaluate,cmap,ref_points)
cell_geo_map = lazy_map(p -> ForwardMap(p), panel_ids)

xyz = lazy_map(evaluate,cell_geo_map,coords)


# # plot()
# # for _xyz in xyz
# #   ids = collect([1,3,4,2,1])

# #   x = map(x->x[1],_xyz)
# #   y = map(x->x[2],_xyz)
# #   z = map(x->x[3],_xyz)
# #   plot!(x[ids],y[ids],z[ids],label=false,marker=:circle)
# # end
# # plot!(show=true)


colors = palette(:tab10)

### Plot individual panel
for pid in collect(1:6)
  p1 = panel_ids .== pid
  xyz1 = xyz[p1]
  plot()
  for i in 1:length(xyz1)
    cellx = xyz1[i]
    out = similar(cellx,VectorValue{2,Float64})

    cellx[1:4]

    x = map(x->x[1],cellx)
    y = map(x->x[2],cellx)
    z = map(x->x[3],cellx)
    r = sqrt.(x.^2 + y.^2 + z.^2)
    θ = rem2pi.(atan.(y, x),RoundDown)
    # _θ = 2*π .+  rem2pi.(atan.(y, x),RoundUp)
    ϕ = asin.(z./r)

    # if there are negative ys
    if any(y .< 0 ) && any(x.>0)
      θ = 2*π .+  rem2pi.(atan.(y, x),RoundUp)
    end
    out = map(θ,ϕ) do θ,ϕ
      VectorValue(θ,ϕ)
    end

    _θ = map(x->x[1],out)
    _ϕ = map(x->x[2],out)
    ids = collect([1,3,4,2,1])
    plot!(_θ[ids],_ϕ[ids],label=false,color=colors[i],marker=:circle)
  end
  plot!(show=true,title="panel $pid")
  savefig(dir*"/panel_$pid")
end

### map all points of the sphere
fi = lazy_map(p->Cartesian2SphericalMap(),panel_ids)
latlon_cell_geo_map = lazy_map(∘, fi, cell_geo_map)
θϕ = lazy_map(evaluate,latlon_cell_geo_map,coords)
# function lazy_collect(cache,a)
#   for i in eachindex(a)
#     getindex!(cache,a,i)
#   end
# end
# cache = array_cache(θϕ)
# bm2() = lazy_collect(cache,θϕ)
# using BenchmarkTools
# @benchmark bm2()


plot()
ids = collect([1,3,4,2,1])
for (i,θϕ) in enumerate(θϕ)
  θ = map(x->x[1],θϕ)
  ϕ = map(x->x[2],θϕ)
  plot!(θ[ids],ϕ[ids],label=false,color=colors[panel_ids[i]],marker=:circle)
end
plot!(show=true)
savefig(dir*"/all_panels")


### plot on vtk
fX(XYZ::VectorValue{3}) = XYZ[1]*XYZ[2]*XYZ[3]
f = panel_to_cartesian(fX)

panel_models  = get_refined_models(n_ref_lvls,true)
panel_model = panel_models[1]
panel_ids = get_panel_ids(panel_model)

Ω_panel = Triangulation(panel_model)
f_panel_cf = ParametricCellField(f,Ω_panel,panel_ids)

# cell_geo_map = geo_map_func(Ω_panel)
# fi = lazy_map(p->Cartesian2SphericalMap(p),panel_ids)
# latlon_cell_geo_map = lazy_map(∘, fi, cell_geo_map)
latlon_cell_geo_map = latlon_geo_map_func(Ω_panel)

panel_cfs = [f_panel_cf,panel_ids]
labels = ["u","pid"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/latlon_model",cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map,order=1)


##### Distributed test

using MPI
using PartitionedArrays
using GridapDistributed

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=6)
panel_model = omodel.parametric_dmodel
Ω_panel = Triangulation(panel_model)
panel_ids = get_panel_ids(panel_model)
f_panel_cf = ParametricCellField(f,Ω_panel,panel_ids)

cell_geo_map = geo_map_func(Ω_panel)
latlon_cell_geo_map = latlon_geo_map_func(Ω_panel)

panel_cfs = [f_panel_cf,get_owned_panel_ids(panel_model)]
labels = ["u","pid"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/distributed_latlon_model",cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map,order=1)
# writevtk(Ω_panel,dir*"/distributed_model",cellfields=cellfields,append=false,geo_map=cell_geo_map,order=1)

# p1 = map(panel_ids) do p
#   p .== 5
# end

# xyz5 = map(local_views(panel_model),p1) do model, p1
#   cmap = get_cell_map(get_grid(model))
#   ref_points = get_cell_ref_coordinates(model)
#   coords = lazy_map(evaluate,cmap,ref_points)
#   ab5 = coords[p1]
#   cell_geo_map = lazy_map(p -> ForwardMap(5), collect(1:sum(p1)))
#   # fi = lazy_map(p->Cartesian2SphericalMap(p),collect(1:sum(p1)))
#   # latlon_cell_geo_map = lazy_map(∘, fi, cell_geo_map)

#   return lazy_map(evaluate,cell_geo_map,ab5)
# end

# θϕ5 = map(local_views(panel_model),p1) do model, p1
#   cmap = get_cell_map(get_grid(model))
#   ref_points = get_cell_ref_coordinates(model)
#   coords = lazy_map(evaluate,cmap,ref_points)
#   ab5 = coords[p1]
#   cell_geo_map = lazy_map(p -> ForwardMap(5), collect(1:sum(p1)))
#   fi = lazy_map(p->Cartesian2SphericalMap(5),collect(1:sum(p1)))
#   latlon_cell_geo_map = lazy_map(∘, fi, cell_geo_map)

#   return lazy_map(evaluate,latlon_cell_geo_map,ab5)
# end

# xx = xyz5.item./1
# _θϕ = θϕ5.item./1

# cell=513
# cellx = xx[cell]
# θϕ = _θϕ[cell]
# x = map(x->x[1],cellx)
# y = map(x->x[2],cellx)
# z = map(x->x[3],cellx)
# θ = map(x->x[1],θϕ)
# diffe = abs(minimum(θ) - maximum(θ))

# idx = abs.(y) .<1e-16
# if any(idx)
#   y[idx].=0.0
# end


# diffe > 3.5 &&  any(y .< 0 )

# if  diffe > 3.5 && any(y .< 0 ) # minimum(θ) == 0.0 &&
#   θ[θ .== 0.0] .= 2*π
# end
# if p == 5 && diffe > 3.5  && any(y .> 0 ) #&& maximum(θ) == 2*π
# # 2pi -> 0.0
# θ[θ .== 2*π] .= 0.0
# end

# minimum(θ) == 0.0 #0.0 -> 2pi
# maximum(θ) == 2*π



# yy = [Vector{Float64}(undef,4 ) for i in 1:length(_θϕ)]
# for (i,cell_θϕ) in enumerate(_θϕ)
#   x = map(x->x[1],cell_θϕ)
#   y = map(x->x[2],cell_θϕ)
#   diffe = abs(minimum(x) - maximum(x))
#   if diffe > 4

#   println("cell $i")
#   println("   θ = $x",)
#   println("   y = $(xx[i])",)
#   end
#   # println("   ϕ = $y",)
#   yy[i] = y
# end








plot()
for cellx in xx
  x = map(x->x[1],cellx)
  y = map(x->x[2],cellx)
  z = map(x->x[3],cellx)
  # plot!(x[ids],y[ids],z[ids],label=false,marker=:circle)
  r = sqrt.(x.^2 + y.^2 + z.^2)
  θ = rem2pi.(atan.(y, x),RoundDown)
  ϕ = asin.(z./r)

    # if there are negative ys
    if any(y .<= 0 ) && any(x.>=0)
      θ = 2*π .+  rem2pi.(atan.(y, x),RoundUp)
    end
    out = map(θ,ϕ) do θ,ϕ
      VectorValue(θ,ϕ)
    end

    _θ = map(x->x[1],out)
    _ϕ = map(x->x[2],out)
    ids = collect([1,3,4,2,1])
    plot!(_θ[ids],_ϕ[ids],label=false,marker=:circle)
end
plot!(show=true)
