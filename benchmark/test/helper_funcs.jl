
function save_ouput(dir,filename,t,out,order)
  keys = ["model","tbench","flops","ops","t","ts","order"]
  vals = [filename,t,out[1],out[2],out[3],out[4],order]
  d = Dict(keys .=> vals)
  safesave(datadir(dir, (filename*"_order$(order).jld2")), d)
end

function get_quadrature(degree)
  quad = Quadrature(QUAD, tensor_product, degree)
  quad.coordinates, quad.weights
end
