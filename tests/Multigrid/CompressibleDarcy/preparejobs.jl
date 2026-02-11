
using Mustache
using DrWatson
using DataStructures

function jobname(args...)
  s = savename(args...;connector="_",equals="_",sort=false,ignores=(:dir,:return_vtk,:driver))
  s = replace(s,"return_vtk_" => "")
  s = replace(s,"dir_" => "")
  return s
end

function clean_params(d)
  o = OrderedDict()
  for k in keys(d)
    if k == :theta
      o[k] = ceil(Int, d[k]*100)
    elseif isa(d[k],Real)
      o[k] = Int(ceil(d[k]))
    elseif isa(d[k],Symbol) || isa(d[k],String) || isa(d[k],Number)
      o[k] = d[k]
    elseif isa(d[k],Tuple)
      o[k] = prod(d[k])
    end
  end
  return o
end

# queue, cpumem, ncputask, ncpunode, walltime
queue_info = Dict{Symbol,Any}(
  :normal => ("normal",4,"4:00:00",48,1),
  :hugemem => ("hugemem",30,"5:00:00",48,1),
  :megamem => ("megamem",60,"12:00:00",48,1),
  :saphire => ("normalsr",4.8,"2:00:00",104,4)
)

function jobdict(name,params,q,np)

  queue, cpumem, walltime, ncpunode, ncputask = queue_info[q]
  project = haskey(ENV,"PROJECT") ? ENV["PROJECT"] : "zg98"

  ntasks = prod(np)
  ncpus  = ntasks * ncputask
  mem = floor(Int, ncpus * cpumem)
  # nnodes = ceil(Int, ncpus / ncpunode)
  # ncpus = nnodes * ncpunode
  # ntaskspernode = min(ncpunode ÷ ncputask,ntasks)

  Dict(
    "project" => project,
    "q" => queue,
    "o" => datadir(name*".o"),
    "e" => datadir(name*".e"),
    "walltime" => walltime,
    # "nnodes" => nnodes,
    "ncpus" => ncpus,
    "mem" => "$(mem)G",
    # "ntasks" => ntasks,
    # "ncputask" => ncputask,
    # "ntaskspernode" => ntaskspernode,
    "name" => name,
    "np" => np,            # Number of processors
    "projectdir" => projectdir(),
    "datadir"    => datadir(),
    "simName"  => name,
    "c" => params[:c],
    "α" => params[:α],
    "n" => params[:n], # Number of cells
    "order" => params[:order],
    "iters" => params[:iters],
    "itu" => params[:itu],
    "itp" => params[:itp],
    "dir" => params[:dir],
    "return_vtk" => params[:return_vtk],
    "driver" => params[:driver],
  )
end


function generate_dictionaries(orders,alphas)
  dicts = OrderedDict[]
  for order in orders
    for α in alphas
      aux = OrderedDict{Symbol,Any}(
                  :c => 0.0,
                  :α => α,
                  :n => 4,
                  :order => order,
                  :iters => 1000, # iters of external FGMRESSolver
                  :itu => 1,# u_solver
                  :itp => 1, #, p_solver
                  :return_vtk => 0,
                  :dir => dir,
                  :driver => driver,
                )
        push!(dicts,aux)
    end
  end
  return dicts
end

############################################

jobsdir = projectdir("gmg_jobs")
!isdir(jobsdir) && mkdir(jobsdir)

dir = datadir("Multigrid")
!isdir(dir) && mkdir(dir)

wdir = projectdir("tests/Multigrid/CompressibleDarcy")
driver = wdir*"/convergence.jl"

np  = 16
queue = :normal
orders = [2]
alphas = [1, 10]

dicts = generate_dictionaries(orders,alphas)

template = read(wdir*"/template.sh",String)

for params in dicts
  fparams = clean_params(params)
  name = jobname(fparams)

  # Generate job file
  jobfile = projectdir(jobsdir*"/",name*".sh")
  open(jobfile,"w") do io
   render(io,template,jobdict(name,params,queue,np))
  end
end
