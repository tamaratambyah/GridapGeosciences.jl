"""
CellField

Is the distributed implementation of CellField.
In such function, we call CellField on the local model and then
recompute the triangulation to ensure proper handling of ghost cells in octree periodic meshes.

"""
function GridapDistributed.CellData.CellField(f::Function,
  trian::GridapDistributed.DistributedTriangulation)

  model = trian.model

  fields = map(local_views(model)) do lmodel
    CellField(f,Triangulation(lmodel))
  end

  trians = map(local_views(model)) do lmodel
    Triangulation(lmodel)
  end

  _trian = GridapDistributed.DistributedTriangulation(trians,model)
  GridapDistributed.DistributedCellField(fields,_trian)
end
