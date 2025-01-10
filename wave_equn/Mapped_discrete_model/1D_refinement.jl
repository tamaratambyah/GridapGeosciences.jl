function Gridap.Adaptivity.refine(model::MyMappedDiscreteModel{Dc}, cell_partition::Int=2) where Dc
  println("my refined model")
  partition = Tuple(fill(cell_partition,Dc))
  return Gridap.Adaptivity.refine(model,partition)
end

### refinement
refined_new_model = Gridap.Adaptivity.refine(new_model,2) # equal to cart_model
