"""
load the last solution from a folder
"""
function load_last(ranks,X,sim_dir::String,simName::String)
  folders = readdir(sim_dir)

  @assert !isempty(folders) "/n No saved solutions"

  f = folders[end]
  t = parse(Float64,f[length(simName)+2:length(f)])
  i_am_main(ranks) && println("restart t = $(t)")
  x =  pload(joinpath(sim_dir,f),ranks)
  return t,FEFunction(X,x)
end

"""
unwrap the iterator, save solution, and findal solution
"""
function unwrap(it,ranks,solT,dir,tF,freq=10)
  sim_dir = dir*"/sim_data"
  final_dir = dir*"/final_solution"

  counter = 1
  while !isnothing(it)
    data, state = it
    t, xh = data

    i_am_main(ranks) && println("t = ", t)

    mod(counter,freq) == 0 && psave(sim_dir*"/solT_$t",xh)

    if t >= tF - Gridap.ODEs.ε
      i_am_main(ranks) && println("Saving final solution")
      psave(final_dir*"/solT_$t",xh)
    end

    counter = counter + 1
    it = iterate(solT, state)
  end

end
