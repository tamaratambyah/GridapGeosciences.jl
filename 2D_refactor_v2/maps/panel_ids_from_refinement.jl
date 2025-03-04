using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Test

function transfer_panel_ids!(panel_id::Vector{Int64},
        n2o_fine::Vector{Int64},n2o_coarse::Vector{Int64},n2o_other::Vector{Int64})

  for i in unique(n2o_fine)
    ids = findall(x->x==i,n2o_other)
    panel_id[ids] .= n2o_coarse[i]
  end

end

function get_n2o_maps(models::Vector{<:DiscreteModel})
  nmodels = length(models)
  n2o_maps = [Int64[] for i in 1:nmodels]

  # the coarsest model is simply the cell id
  copy!(n2o_maps[1], collect(1:num_cells(models[1])) )

  for i = 2:nmodels
    glue = get_adaptivity_glue(models[i])
    n2o_faces_map = glue.n2o_faces_map[end]

    copy!(n2o_maps[i], n2o_faces_map)
  end

  n2o_maps
end

function get_panel_ids(models::Vector{<:DiscreteModel})
  nmodels = length(models)

  n2o_maps = get_n2o_maps(models) # get the n2o_maps for each model in models
  panel_ids = [Int64[] for i in 1:nmodels] # global array of panel_ids for each cell

  # models 1,2 are simply the n2o_maps
  copy!(panel_ids[1],n2o_maps[1])
  copy!(panel_ids[2],n2o_maps[2])


  # this "3" is hardcoded since 1,2 are dealt with above
  # lvls == represents number of times to apply transfer_panel_ids
  lvls = collect(1:nmodels) - 3*ones(Int64,nmodels)

  # initialise work arrays
  n2o_fine = Int64[]
  n2o_other = Int64[]
  n2o_coarse = Int64[]
  panel_id = Int64[]

  for model_id = 3:nmodels

    copy!(panel_id,n2o_maps[model_id]) # start with the n2o_map for the model

    # apply transfer_panel_ids the appropriate number of times
    for nlvl in 0:lvls[model_id]
      i, j = nlvl, nlvl+1
      copy!(n2o_fine,n2o_maps[model_id-i])
      copy!(n2o_other,panel_id)
      copy!(n2o_coarse, n2o_maps[model_id-j])

      transfer_panel_ids!(panel_id,n2o_fine,n2o_coarse,n2o_other)

    end

    copy!(panel_ids[model_id],panel_id)
  end

  panel_ids

end

### Test with CartesianDiscreteModel
coarse_model = CartesianDiscreteModel((0,1,0,1),(2,2))
model = Gridap.Adaptivity.refine(coarse_model)
fine_model = Gridap.Adaptivity.refine(model)
super_fine_model = Gridap.Adaptivity.refine(fine_model)
super_duper_fine_model = Gridap.Adaptivity.refine(super_fine_model)

# list of models going coarse -> fine
models = [coarse_model, model, fine_model, super_fine_model, super_duper_fine_model]
panel_ids = get_panel_ids(models)

# test there are no panel_ids > num_cells(coarse_model)
for i in 1:length(panel_ids)
  tmp = ( panel_ids[i] .<= num_cells(coarse_model) )
  @test sum(tmp) == length(panel_ids[i]  )
end


### Test with cube_model
include("cube_surface_1_cell_per_panel_2D.jl")
coarse_model = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel_2D()...)
model = Gridap.Adaptivity.refine(coarse_model)
fine_model = Gridap.Adaptivity.refine(model)
super_fine_model = Gridap.Adaptivity.refine(fine_model)
super_duper_fine_model = Gridap.Adaptivity.refine(super_fine_model)

# list of models going coarse -> fine
models = [coarse_model, model, fine_model, super_fine_model, super_duper_fine_model]
panel_ids = get_panel_ids(models)

# test there are no panel_ids > num_cells(coarse_model)
for i in 1:length(panel_ids)
  tmp = ( panel_ids[i] .<= num_cells(coarse_model) )
  @test sum(tmp) == length(panel_ids[i]  )
end



#################### debugging #####################################


# glue = get_adaptivity_glue(model)
# fine_glue = get_adaptivity_glue(fine_model)
# super_fine_glue = get_adaptivity_glue(super_fine_model)
# super_duper_fine_glue = get_adaptivity_glue(super_duper_fine_model)

# coarse_n2o_faces_map = collect(1:num_cells(coarse_model))
# n2o_faces_map = glue.n2o_faces_map[end]
# fine_n2o_faces_map = fine_glue.n2o_faces_map[end]
# super_fine_n2o_faces_map = super_fine_glue.n2o_faces_map[end]
# super_duper_fine_n2o_faces_map = super_duper_fine_glue.n2o_faces_map[end]

# # list of n2o_maps:  going coarse -> fine
# n2o_maps = [coarse_n2o_faces_map, n2o_faces_map, fine_n2o_faces_map,
#             super_fine_n2o_faces_map, super_duper_fine_n2o_faces_map]


# model_id = 4

# copy!(n2o_fine,n2o_maps[model_id])
# copy!(n2o_other,n2o_maps[model_id])
# copy!(n2o_coarse, n2o_maps[model_id-1])

# panel_id = Vector{eltype(n2o_fine)}(undef,length(n2o_fine))

# transfer_panel_ids!(panel_id,n2o_fine,n2o_coarse,n2o_other)

# copy!(n2o_fine,n2o_maps[model_id-1])
# copy!(n2o_other,panel_id)
# copy!(n2o_coarse, n2o_maps[model_id-2])
# transfer_panel_ids!(panel_id,n2o_fine,n2o_coarse,n2o_other)

# copy!(panel_ids[model_id],panel_id)



# for i in unique(super_fine_n2o_faces_map)
#   ids = findall(x->x==i,super_fine_n2o_faces_map)
#   panel_id[ids] .= fine_n2o_faces_map[i]
# end

# _panel_id = copy(panel_id)
# for i in unique(fine_n2o_faces_map)
#   ids = findall(x->x==i,_panel_id)
#   panel_id[ids] .= n2o_faces_map[i]
# end


# for i in unique(fine_n2o_faces_map)
#   ids = findall(x->x==i,fine_n2o_faces_map)
#   panel_id[ids] .= n2o_faces_map[i]
# end
