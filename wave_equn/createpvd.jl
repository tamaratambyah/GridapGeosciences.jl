function make_pvd(dirs::Vector{String},simName::String,outname::String,convert2seconds)
  for dir in dirs
    println(dir)
    _make_pvd(dir,simName,outname,convert2seconds)
  end

end

function make_pvd(dir::String,simName::String,outname::String,convert2seconds)
  _make_pvd(dir,simName,outname,convert2seconds)
end

function _make_pvd(dir_loc,simName,outname,convert2seconds)

  files = filter(x->endswith(x, ".vtu"), readdir((dir_loc)))


  top = ["<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",
        "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n",
    "<Collection>\n"]

  tail = ["</Collection>\n",
          "</VTKFile>\n"]

  output = String[]
  times = []

  for f in files
    str = f
    time = str[length(simName)+1:length(str)-4]
    println(time)
    t = parse(Float64,time)/(convert2seconds)
    t2 = string(t)
    # println(t)
    my_string = "<DataSet timestep=\"$t2\" part=\"0\" file=\"$str\"/>\n"
    push!(output,my_string)
    push!(times,parse(Float64,time))

  end

  p = sortperm(times)
  output_sorted = copy(output[p])

  final = []
  append!(final,top)
  append!(final,output_sorted)
  append!(final,tail)

  open(dir_loc*"/$outname.pvd", "w") do file
      for str in final
        write(file, str)
      end
  end


end
