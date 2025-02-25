using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity


coarse_model = CartesianDiscreteModel((0,1,0,1),(2,2))
model = Gridap.Adaptivity.refine(coarse_model)
fine_model = Gridap.Adaptivity.refine(model)
super_fine_model = Gridap.Adaptivity.refine(fine_model)

glue = get_adaptivity_glue(model)
fine_glue = get_adaptivity_glue(fine_model)
super_fine_glue = get_adaptivity_glue(super_fine_model)

coarse_n2o_faces_map = collect(1:num_cells(coarse_model))

n2o_faces_map = glue.n2o_faces_map[end]
o2n_faces_map = glue.o2n_faces_map

fine_n2o_faces_map = fine_glue.n2o_faces_map[end]
fine_o2n_faces_map = fine_glue.o2n_faces_map


panel_id = Vector{Int}(undef,length(fine_n2o_faces_map))

for i in unique(fine_n2o_faces_map)
  ids = findall(x->x==i,fine_n2o_faces_map)
  panel_id[ids] .= n2o_faces_map[i]
end

_panel_id = copy(panel_id)
for i in unique(n2o_faces_map)
  ids = findall(x->x==i,_panel_id)
  panel_id[ids] .= coarse_n2o_faces_map[i]
end



super_fine_n2o_faces_map = super_fine_glue.n2o_faces_map[end]
super_fine_o2n_faces_map = super_fine_glue.o2n_faces_map

#  collect1d(fine_o2n_faces_map)

panel_id = Vector{Int}(undef,length(super_fine_n2o_faces_map))

for i in unique(super_fine_n2o_faces_map)
  ids = findall(x->x==i,super_fine_n2o_faces_map)
  panel_id[ids] .= fine_n2o_faces_map[i]
end

_panel_id = copy(panel_id)
for i in unique(fine_n2o_faces_map)
  ids = findall(x->x==i,_panel_id)
  panel_id[ids] .= n2o_faces_map[i]
end










for i in unique(fine_n2o_faces_map)
  ids = findall(x->x==i,fine_n2o_faces_map)
  panel_id[ids] .= n2o_faces_map[i]
end
