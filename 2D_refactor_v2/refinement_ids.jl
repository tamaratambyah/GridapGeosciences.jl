using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity

function transfer_panel_ids!(panel_id,n2o_fine,n2o_coarse,n2o_other)
  for i in unique(n2o_fine)
    ids = findall(x->x==i,n2o_other)
    panel_id[ids] .= n2o_coarse[i]
  end
end

coarse_model = CartesianDiscreteModel((0,1,0,1),(2,2))
model = Gridap.Adaptivity.refine(coarse_model)
fine_model = Gridap.Adaptivity.refine(model)
super_fine_model = Gridap.Adaptivity.refine(fine_model)
super_duper_fine_model = Gridap.Adaptivity.refine(super_fine_model)

glue = get_adaptivity_glue(model)
fine_glue = get_adaptivity_glue(fine_model)
super_fine_glue = get_adaptivity_glue(super_fine_model)
super_duper_fine_glue = get_adaptivity_glue(super_duper_fine_model)

coarse_n2o_faces_map = collect(1:num_cells(coarse_model))
n2o_faces_map = glue.n2o_faces_map[end]
fine_n2o_faces_map = fine_glue.n2o_faces_map[end]
super_fine_n2o_faces_map = super_fine_glue.n2o_faces_map[end]
super_duper_fine_n2o_faces_map = super_duper_fine_glue.n2o_faces_map[end]

# list of n2o_maps:  going coarse -> fine
n2o_maps = [coarse_n2o_faces_map, n2o_faces_map, fine_n2o_faces_map,
            super_fine_n2o_faces_map, super_duper_fine_n2o_faces_map]

nmodels = length(n2o_maps)



panel_ids = [Int64[] for i in 1:nmodels]
copy!(panel_ids[1],n2o_maps[1])
copy!(panel_ids[2],n2o_maps[2])

n2o_fine = Int64[]
n2o_other = Int64[]
n2o_coarse = Int64[]

nlvl_of_refinement = [3,3,4]

for model_id = 3:nmodels

  copy!(n2o_fine, n2o_maps[model_id])
  copy!(n2o_other, n2o_maps[model_id])
  copy!(n2o_coarse, n2o_maps[model_id-1])

  panel_id = Vector{eltype(n2o_fine)}(undef,length(n2o_fine))

  transfer_panel_ids!(panel_id,n2o_fine,n2o_coarse,n2o_other)

  copy!(n2o_fine,n2o_maps[model_id-1])
  copy!(n2o_other,panel_id)
  copy!(n2o_coarse, n2o_maps[model_id-2])

  transfer_panel_ids!(panel_id,n2o_fine,n2o_coarse,n2o_other)
  copy!(panel_ids[model_id],panel_id)

end
model_id = 5
copy!(n2o_fine,n2o_maps[model_id-2])
copy!(n2o_other,panel_ids[model_id])
copy!(n2o_coarse, n2o_maps[model_id-3])

transfer_panel_ids!(panel_ids[model_id],n2o_fine,n2o_coarse,n2o_other)




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










for i in unique(fine_n2o_faces_map)
  ids = findall(x->x==i,fine_n2o_faces_map)
  panel_id[ids] .= n2o_faces_map[i]
end
