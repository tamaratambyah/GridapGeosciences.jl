
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
  :normal => ("normal",4,"48:00:00",48,1),
  :hugemem => ("hugemem",30,"12:00:00",48,1),
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
    "mem" => "2990G", #"$(mem)G",
    # "ntasks" => ntasks,
    # "ncputask" => ncputask,
    # "ntaskspernode" => ntaskspernode,
    "name" => name,
    "np" => np,            # Number of processors
    "projectdir" => projectdir(),
    "datadir"    => datadir(),
    "simName"  => name,
    "n" => params[:n], # refinement level
    "order" => params[:order],
    "dir" => params[:dir],
    "return_vtk" => params[:return_vtk],
    "driver" => params[:driver],
  )
end


function generate_dictionaries(orders,nlevels,nprocs)
  dicts = OrderedDict[]
  for order in orders
      for (nl,np) in zip(nlevels,nprocs)
            aux = OrderedDict{Symbol,Any}(
                        :n => nl, ## on cubed sphere meshes
                        :order => order,
                        :return_vtk => 1,
                        :dir => dir,
                        :driver => driver,
                        :np => np,  # number of processors
                      )
            push!(dicts,aux)
          end
        end

  return dicts
end

############################################

jobsdir = projectdir("hodge_jobs_thick")
!isdir(jobsdir) && mkdir(jobsdir)

dir = datadir("HodgeLaplacianConvergence_30quad_thick")
!isdir(dir) && mkdir(dir)

wdir = projectdir("tests/HodgeLaplacian")
driver = wdir*"/HodgeLaplacian.jl"


queue = :megamem
orders = [0, 1]
nlevels = [1,2,3,4,5] # cubed sphere mesh
nprocs  = [4,4,16,32,48]

dicts = generate_dictionaries(orders,nlevels,nprocs)

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
