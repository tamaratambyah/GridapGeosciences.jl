
using Gridap
using GridapGeosciences
using DrWatson
using FillArrays

function MapCubeToSphere(xyz)
  radius = 1
  x,y,z = xyz
  xₛ = x*sqrt(1.0-y^2/2-z^2/2+y^2*z^2/(3.0))
  yₛ = y*sqrt(1.0-z^2/2-x^2/2+x^2*z^2/(3.0))
  zₛ = z*sqrt(1.0-x^2/2-y^2/2+x^2*y^2/(3.0))
  radius*Point(xₛ,yₛ,zₛ)
end

n = 4
domain = (-1,1,-1,1,-1,1)
cells  = (n,n,n)
model  = CartesianDiscreteModel(domain,cells)

writevtk(model,datadir("models")*"/cube",append=false)

panel_to_entity=[23,26,24,25,21,22]
entity_to_panel=Dict(23=>1,26=>2,24=>3,25=>4,21=>5,22=>6)
ptr_ncells_panel=zeros(Int64,length(panel_to_entity)+1)

labels = get_face_labeling(model)
bgface_to_mask = Gridap.Geometry.get_face_mask(labels,"boundary",2) # extract index of boundary
Γface_to_bgface = findall(bgface_to_mask) # get the index of those on the boundary
face_to_entity  = labels.d_to_dface_to_entity[3][Γface_to_bgface] # get the faces


# assign the face cells to a panel
for entity in face_to_entity
  panel=entity_to_panel[entity]
  ptr_ncells_panel[panel+1]=ptr_ncells_panel[panel+1]+1 # number of face cells on each panel
end
ptr_ncells_panel[1]=1

for i=1:length(ptr_ncells_panel)-1
  ptr_ncells_panel[i+1]=ptr_ncells_panel[i+1]+ptr_ncells_panel[i]
end

Γface_to_bgface_panelwise=similar(Γface_to_bgface)
for (i,entity) in enumerate(face_to_entity)
  println(entity)
  panel=entity_to_panel[entity]
  Γface_to_bgface_panelwise[ptr_ncells_panel[panel]]=Γface_to_bgface[i]
  ptr_ncells_panel[panel]=ptr_ncells_panel[panel]+1
end

Γface_to_bgface_panelwise

cube_surface_trian = BoundaryTriangulation(model,Γface_to_bgface_panelwise)
writevtk(cube_surface_trian,datadir("models")*"/cube_surface_trian",append=false)

order = 1

# Generate high-order FE map and ordering
vector_reffe=ReferenceFE(lagrangian,VectorValue{3,Float64},order)
V = FESpace(cube_surface_trian,vector_reffe; conformity=:H1)
vh = interpolate(MapCubeToSphere,V)
writevtk(cube_surface_trian,datadir("models")*"/cube_surface_trian",
          cellfields=["vh"=>vh],append=false)

evaluate((vh),get_cell_points(cube_surface_trian))


scalar_reffe=ReferenceFE(QUAD,lagrangian,Float64,order)
xref=Gridap.ReferenceFEs.get_node_coordinates(scalar_reffe)
xrefₖ=Fill(xref,num_cells(cube_surface_trian))
vhx=lazy_map(evaluate,Gridap.CellData.get_data(vh),xrefₖ)


V = FESpace(cube_surface_trian,scalar_reffe; conformity=:H1)
node_coordinates = Vector{Point{3,Float64}}(undef,num_free_dofs(V))
cell_node_ids    = get_cell_dof_ids(V)

dof_vector = Vector{Point{3,Float64}}(undef,num_free_dofs(V))
cell_vector = lazy_map(evaluate,Gridap.CellData.get_data(vh),xrefₖ)

cache_cell_node_ids = array_cache(cell_node_ids)
cache_cell_vector   = array_cache(cell_vector)
for k=1:length(cell_node_ids)
    current_node_ids = getindex!(cache_cell_node_ids,cell_node_ids,k)
    current_values   = getindex!(cache_cell_vector,cell_vector,k)
    for (i,id) in enumerate(current_node_ids)
    dof_vector[current_node_ids[i]]=current_values[i]
    end
end
node_coordinates = copy(dof_vector)



cell_types  = collect(Fill(1,num_cells(cube_surface_trian)))
cell_reffes = [scalar_reffe]





cube_surface_grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                                      Gridap.Arrays.Table(cell_node_ids),
                                                      cell_reffes,
                                                      cell_types,
                                                      Gridap.Geometry.Oriented())

cube_surface_model = Gridap.Geometry.compute_active_model(cube_surface_trian)
writevtk(cube_surface_model,datadir("models")*"/cube_surface_model",
          append=false)

topology = Gridap.Geometry.get_grid_topology(cube_surface_model)
labeling = Gridap.Geometry.get_face_labeling(cube_surface_model)

cubed_sphere_linear_model = Gridap.Geometry.UnstructuredDiscreteModel(cube_surface_grid,topology,labeling)

writevtk(cubed_sphere_linear_model,datadir("models")*"/cube_sphere_linear_model",
          append=false)




m1=Fill(Gridap.Fields.GenericField(MapCubeToSphere),num_cells(cube_surface_trian))
m2=get_cell_map(cube_surface_trian)
map=lazy_map(∘,m1,m2)


_model = CubedSphereDiscreteModel(n)
_trian = Triangulation(_model)

pts = get_cell_points(Measure(_trian,4))

cmaps = Gridap.Arrays.collect1d(get_cell_map(_trian))
_Jt = lazy_map(∇,cmaps)
_Jtx = lazy_map(evaluate,_Jt,Gridap.CellData.get_data(pts))./1

Jt    = lazy_map(∇,map)



# Get the index of the panel for each element
fl = get_face_labeling(cubed_sphere_linear_model)
panel_id  = fl.d_to_dface_to_entity[3]
sign_flip = [panel_id[i] == 25 || panel_id[i] == 21 || panel_id[i] == 24 for i=1:length(panel_id)]
fsign_flip = lazy_map(Gridap.Fields.ConstantField,sign_flip)

a = (1,2,3,4,5,6)
b = TensorValue{2,3}(a)
b[1,1]
b[1,2]
b[1,3]

b[2,1]
b[2,2]
b[2,3]


using StaticArrays
mat = MMatrix{2,3}(a)
_mat = SMatrix{2,3}(a)

function _unit_outward_normal(v::Gridap.Fields.MultiValue{Tuple{2,3}},sign_flip::Bool)
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  n = VectorValue(n1,n2,n3)
  (-1)^sign_flip*n/norm(n)
end

lazy_map(Operation(_unit_outward_normal),Jt,fsign_flip)

nc=num_cells(cubed_sphere_linear_model)
tface_to_mface=Gridap.Fields.IdentityVector(nc)
tface_to_mface_map=Fill(Gridap.Fields.GenericField(identity),nc)
mface_to_tface=tface_to_mface
Gridap.Geometry.FaceToFaceGlue(tface_to_mface,tface_to_mface_map,mface_to_tface)
