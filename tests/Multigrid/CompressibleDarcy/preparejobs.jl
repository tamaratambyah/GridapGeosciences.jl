
using Mustache
using DrWatson
using DataStructures

function jobname(args...)
  s = savename(args...;connector="_",equals="_",sort=false,
              ignores=(:dir,:return_vtk,:driver,:np))
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

function jobdict(name,params,q)
  np = params[:np] # number of processors

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
    "o" => jobsdir*"/$(name).o",
    "e" => jobsdir*"/$(name).e",
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
    "α" => params[:alpha],
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


function generate_dictionaries(orders,alphas,nlevels,exiters,itus,itps,nprocs,compress)
  dicts = OrderedDict[]
  for order in orders
    for c in compress
    for α in alphas
      for (nl,np) in zip(nlevels,nprocs)
        for (iters,itu,itp) in zip(exiters,itus,itps)
            aux = OrderedDict{Symbol,Any}(
                        :c => c,
                        :alpha => α,
                        :n => 2^(nl),
                        :order => order,
                        :iters => iters, # iters of external FGMRESSolver
                        :itu => itu,# u_solver
                        :itp => itp, #, p_solver
                        :return_vtk => 1,
                        :dir => dir,
                        :driver => driver,
                        :np => np,  # number of processors
                      )
            push!(dicts,aux)
          end
        end
      end
    end
  end
  return dicts
end

############################################

jobsdir = projectdir("gmg_jobs")
!isdir(jobsdir) && mkdir(jobsdir)

dir = datadir("Multigrid_compressible_p1")
!isdir(dir) && mkdir(dir)

wdir = projectdir("tests/Multigrid/CompressibleDarcy")
driver = wdir*"/convergence.jl"


queue = :normal
orders = [1]
compress = [1,10,100]
alphas = [1,10,100]
nlevels = [3,4,5,6]
nprocs  = [4,4,16,16]
exiters = [1000,20] # kyrol iters
itus = [0,1] # iterate u via gmg
itps = [0,1] # iterate p via cg

dicts = generate_dictionaries(orders,alphas,nlevels,exiters,itus,itps,nprocs,compress)

template = read(wdir*"/template.sh",String)

for params in dicts
  fparams = clean_params(params)
  name = jobname(fparams)

  # Generate job file
  jobfile = projectdir(jobsdir*"/",name*".sh")
  open(jobfile,"w") do io
   render(io,template,jobdict(name,params,queue))
  end
end

# submit job file
for job in readdir(jobsdir)
  println("Submitting job: ", job)
  run(`qsub $jobsdir/$job`)
end


println("--DONE--")
