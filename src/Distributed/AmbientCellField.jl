"""
AmbientCellField

Is the distributed implementation of AmbientCellField.
In such function, we call AmbientCellField on the local model and then
recompute the triangulation to ensure proper handling of ghost cells in octree periodic meshes.

"""

function AmbientCellField(f::Function,
  trian::GridapDistributed.DistributedTriangulation)

  panel_model = trian.model

  fields = map(local_views(panel_model)) do lmodel
    AmbientCellField(f,Triangulation(lmodel))
  end

  trians = map(local_views(panel_model)) do lmodel
    Triangulation(lmodel)
  end

  _trian = GridapDistributed.DistributedTriangulation(trians,panel_model)
  GridapDistributed.DistributedCellField(fields,_trian)
end
