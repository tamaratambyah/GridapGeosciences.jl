"""
AmbientCellField

A AmbientCellField is returns an GenericCellField object, where the cell_field is an
array of cell-wise functions.
The user must define a function that given an inverse map, returns a function
that takes coordinates in ambient space as input, and returns the value of the function.
This means AmbientCellField is different to CellField, where the user passes
a function that takes points in physical space and returns the function evaluated in physical space.

Example usage is:

```
function ambient_sgrad(f::Function,inverse_map::Field)
  function _gradf(x)
    αβ = inverse_map(x)
    m = ForwardMap(inverse_map.panel,inverse_map.radius)
    sgrad(panel_f(f),m)(αβ)
  end
end


f_cf = AmbientCellField(ambient_sgrad,Ω_ambient)
```


The input function `f` takes an InverseMap object, and returns another function `_f`.
Such function takes points in the ambient coordinate system of the manifold, and returns the
ambient field evaluated at parametric points.
Note, this will only work for manifolds where the inverse map can be defined analytically


"""


function AmbientCellField(f::Function,
  trian::BodyFittedTriangulation{Dc,Dp,<:CubedSphereAmbientDiscreteModel}) where {Dc,Dp}

  ambient_model = get_background_model(trian)
  panel_model = ambient_model.panel_model
  panel_ids = get_panel_ids(panel_model)
  @check length(panel_ids) == num_cells(trian) "\n Incorrect panel ids"

  inv_map_generator = get_inverse_map_generator(ambient_model)
  inverse_maps = lazy_map(inv_map_generator,panel_ids)
  cell_field = lazy_map(m->GenericField(f(m)),inverse_maps)
  CellData.GenericCellField(cell_field,trian,PhysicalDomain())
end
