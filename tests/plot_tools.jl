using ScatteredInterpolation
using GLMakie

include("output_tools.jl")

function plot_latlon(dir::String,n,plotName::String,cf::String)

  ## open the mesh file and get the x,y points
  mesh_file = joinpath(dir,"coords.csv")
  xs = get_scalar_field_from_csv(mesh_file, "xs")
  ys = get_scalar_field_from_csv(mesh_file, "ys")
  points = transpose([xs ys])

  ## make the grid for interpolation
  x = range(minimum(xs), stop = maximum(xs), length = n)
  y = range(minimum(ys), stop = maximum(ys), length = n)
  X = repeat(x, n)[:]
  Y = repeat(y', n)[:]
  gridPoints = [X Y]'

  ## get all the output_t.csv files
  _files = filter(x->endswith(x, ".csv"), readdir((dir)))
  files = _files[2:end] # remove coords

  _dir = dir*"/figs"
  !isdir(_dir) && mkdir(_dir)

  ## for each file,
  ## 1. extract t from the file name
  ## 2. get the scalar field and interpolate onto the plotting grid
  ## 3. plot with Makie and save
  for (i,file) in enumerate(files)
    t = file[9:end-4]
    output_file = joinpath(dir,file)
    cdata = get_scalar_field_from_csv(output_file, cf)

    itp = ScatteredInterpolation.interpolate(NearestNeighbor(), points, cdata);
    interpolated = ScatteredInterpolation.evaluate(itp, gridPoints)
    gridded = reshape(interpolated, n, n)

    f = Figure()
    ax = Axis(f[1, 1], limits = (minimum(xs), maximum(xs), minimum(ys), maximum(ys)),
    xlabel = "longitude", ylabel = "latitude",
        title = plotName*": t = $t")
    GLMakie.surface!(ax,x, y, zeros(size(gridded)),color=gridded)
    Colorbar(f[1,2], limits = (minimum(cdata), maximum(cdata)), colormap = :viridis)
    save(_dir*"/output_$i.png", f)
  end

end
