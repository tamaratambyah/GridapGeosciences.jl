using Gridap.FESpaces
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.Arrays
using FillArrays

## Gridap Fixes

function Gridap.ReferenceFEs._Nedelec_face_moments(p, fshfs, c_fips, fcips, fwips)
  nc = length(c_fips)
  cfshfs = fill(fshfs, nc)
  cvals = lazy_map(evaluate,cfshfs,c_fips)

  fvs = Gridap.ReferenceFEs._nfaces_vertices(Float64,p,num_dims(p)-1)
  fts = [hcat([vs[2]-vs[1]...],[vs[3]-vs[1]...]) for vs in fvs]

  # Ref facet FE functions evaluated at the facet integration points (in ref facet)
  cvals = [fwips[i].*cvals[i] for i in 1:nc]

  fns = Gridap.ReferenceFEs.get_facet_normal(p)
  # @santiagobadia : Temporary hack for making it work for structured hex meshes
  ft = eltype(fns)
  cvals = [ Gridap.ReferenceFEs._broadcast_extend(ft,Tm,b) for (Tm,b) in zip(fts,cvals)]
  cvals = [ Gridap.ReferenceFEs._broadcast_cross(ft,n,b) for (n,b) in zip(fns,cvals)]
  return cvals
end


function Gridap.FESpaces.get_sign_flip(model::DiscreteModel,
                       cell_reffe::AbstractArray{<:GenericRefFE{Nedelec}})
    Gridap.Arrays.lazy_map(Gridap.FESpaces.SignFlipMap(model),
            cell_reffe,
            IdentityVector(Int32(num_cells(model))))
end


function Gridap.FESpaces.get_cell_dof_basis(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:GenericRefFE{Nedelec}},
  ::CurlConformity,
  sign_flip=Gridap.FESpaces.get_sign_flip(model, cell_reffe))
  cell_map  = get_cell_map(Triangulation(model))
  phi       = cell_map[1]
  reffe     = cell_reffe[1]
  Dc        = num_dims(reffe)
  et        = eltype(return_type(get_prebasis(reffe)))
  pt        = Point{Dc,et}
  Dp        = first(size(return_type(phi,zero(pt))))
  cell_dofs = lazy_map(get_dof_basis,cell_reffe)
  cell_ownids = lazy_map(get_face_own_dofs,cell_reffe)
  cell_map = get_cell_map(Triangulation(model))
  cell_Jt = lazy_map(Broadcasting(∇),cell_map)
  cell_x = lazy_map(get_nodes,cell_dofs)
  cell_Jtx  = lazy_map(evaluate,cell_Jt,cell_x)
  k = TransformNedelecDofBasisFixed{Dc,Dp}()
  lazy_map(k,cell_dofs,cell_Jtx,cell_ownids,sign_flip)
end

struct TransformNedelecDofBasisFixed{Dc,Dp} <: Map end;

  function Gridap.Fields.return_cache(::TransformNedelecDofBasisFixed{Dc,Dp},
                        dofs,
                        Jtx,
                        ownids,
                        ::AbstractVector{Bool}) where {Dc,Dp}
    # Assumes the same element type in all the mesh
    nodes, nf_nodes, nf_moments = get_nodes(dofs),get_face_nodes_dofs(dofs),get_face_moments(dofs)
    face_moments = [ similar(i,VectorValue{Dp,typeof(dofs.nodes[1][1])})  for i in nf_moments ]
    db = MomentBasedDofBasis(nodes,face_moments,nf_nodes)
    db
  end

function Gridap.Fields.evaluate!(cache,
                   ::TransformNedelecDofBasisFixed{Dc,Dp},
                   dofs,
                   Jtx,
                   ownids,
                   sign_flip::AbstractVector{Bool}) where {Dc,Dp}
  basis = cache # Assumes the same element type in all the mesh

  for face in 1:length(ownids)
    face_dofs_ids = ownids[face]
    face_point_ids = dofs.face_nodes[face]
    face_moments_out = basis.face_moments[face]
    face_moments_in = dofs.face_moments[face]
    for p in 1:length(face_point_ids)
      # F = transpose(Jtx[p])
      F = transpose(Jtx[face_point_ids[p]])
      for i in 1:length(face_dofs_ids)
        sign = (-1)^sign_flip[face_dofs_ids[i]]
        face_moments_out[p,i] = sign * (F⋅face_moments_in[p,i])
      end
    end
  end
  basis
end


struct CoVariantPiolaMapFixed <: Map end

function Gridap.Fields.evaluate!(
  cache,
  ::Broadcasting{typeof(∇)},
  a::Fields.BroadcastOpFieldArray{CoVariantPiolaMapFixed})
  v, Jt, sign_flip = a.args
  # Assuming J comes from an affine map
  ∇v = Broadcasting(∇)(v)
  k = CoVariantPiolaMapFixed()
  Broadcasting(Operation(k))(∇v,Jt,sign_flip)
end

function Gridap.Arrays.lazy_map(
  ::Broadcasting{typeof(gradient)},
  a::LazyArray{<:Fill{Broadcasting{Operation{CoVariantPiolaMapFixed}}}})
  v, Jt, sign_flip = a.args
  ∇v = lazy_map(Broadcasting(∇),v)
  k = CoVariantPiolaMapFixed()
  lazy_map(Broadcasting(Operation(k)),∇v,Jt,sign_flip)
end

function Gridap.Fields.evaluate!(cache,::CoVariantPiolaMapFixed,
                   v::Number,
                   Jt::Number,
                   sign_flip::Bool)
   (((-1)^sign_flip)*v)⋅(transpose(pinvJt(Jt)))# we multiply by the right side to compute the gradient correctly
end

function Gridap.Fields.evaluate!(cache,
                   k::CoVariantPiolaMapFixed,
                   v::AbstractVector{<:Field},
                   phi::Field,
                   sign_flip::AbstractVector{<:Field})
  Jt = ∇(phi)
  Broadcasting(Operation(k))(v,Jt,sign_flip)
end

function Gridap.Arrays.lazy_map(
  k::CoVariantPiolaMapFixed,
  cell_ref_shapefuns::AbstractArray{<:AbstractArray{<:Field}},
  cell_map::AbstractArray{<:Field},
  sign_flip::AbstractArray{<:AbstractArray{<:Field}})

  cell_Jt = lazy_map(∇,cell_map)
  lazy_map(Broadcasting(Operation(k)),cell_ref_shapefuns,cell_Jt,sign_flip)
end


function Gridap.FESpaces.get_cell_shapefuns(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:GenericRefFE{Nedelec}},
  ::CurlConformity,
  sign_flip=Gridap.FESpaces.get_sign_flip(model, cell_reffe))

  cell_reffe_shapefuns = lazy_map(get_shapefuns,cell_reffe)
  cell_map = get_cell_map(Triangulation(model))
  k = CoVariantPiolaMapFixed()
  lazy_map(k,
           cell_reffe_shapefuns,
           cell_map,
           lazy_map(Broadcasting(constant_field), sign_flip))
end

## GridapDistributed Fixes
function Gridap.FESpaces.FESpace(model::GridapDistributed.DistributedDiscreteModel,
                 reffe::GenericRefFE{Nedelec};
                 conformity=nothing,
                 split_own_and_ghost=false,
                 constraint=nothing,
                 kwargs...)
  cell_reffes = map(local_views(model)) do m
    Fill(reffe,num_cells(m))
  end
  GridapDistributed._common_fe_space_constructor(model,cell_reffes;conformity,split_own_and_ghost,kwargs...)
end

function FESpaces.FESpace(model::GridapDistributed.DistributedDiscreteModel,
                          reffe::Tuple{Nedelec,Any,Any};
                          conformity=nothing,
                          split_own_and_ghost=false,
                          constraint=nothing,
                          kwargs...)

  cell_reffes = map(local_views(model)) do m
    basis,reffe_args,reffe_kwargs = reffe
    cell_reffe = ReferenceFE(m,basis,reffe_args...;reffe_kwargs...)
  end
  GridapDistributed._common_fe_space_constructor(model,cell_reffes;conformity,split_own_and_ghost,kwargs...)
end
