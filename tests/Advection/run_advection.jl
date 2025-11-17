using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed

using DrWatson


using Gridap.CellData, Gridap.Geometry

# include("AdvectionDGUpwinding.jl")
# include("../convergence_tools.jl")
# include("advection_funcs.jl")

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

include("TransientAdvectionDGUpwinding.jl")
with_mpi() do distribute
  TransientAdvectionDGUpwinding.main(distribute,nprocs;octree=true)
end




# n_ref_lvls = 4
# ps = [1]#,2,3]
# ls = LUSolver()

# vX = panel_to_cartesian(tangent_vec(vecX))
# u = panel_to_cartesian(u0)
# uvX = panel_to_cartesian(u0vecX)

# dir = datadir("AdvectionDGConvergence")
# (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)


# models =  get_octree_refined_models(ranks,n_ref_lvls)
# # models  = get_distributed_refined_models(ranks,nprocs,n_ref_lvls,false)



# p_convergence_test(ranks,ps,models,AdvectionDGUpwinding.advection_dg_solver,dir,u,vX,uvX,ls,false)


# using CSV
# using DataFrames
# using SparseArrays
# using DrWatson
# using LinearAlgebra
# dir = datadir("AdvectionDGConvergence")
# Mcsv=CSV.read(dir*"/A_proc1.csv",DataFrame)
# M1=Array(sparse(Mcsv[!,:I],Mcsv[!,:J],Mcsv[!,:V]))

# Mcsv=CSV.read(dir*"/A_proc2.csv",DataFrame)
# M2=Array(sparse(Mcsv[!,:I],Mcsv[!,:J],Mcsv[!,:V]))

# M1 == M2
# norm(M1-M2)

# Mcsv=CSV.read(dir*"/A_proc1_octree.csv",DataFrame)
# octree_M1=Array(sparse(Mcsv[!,:I],Mcsv[!,:J],Mcsv[!,:V]))

# Mcsv=CSV.read(dir*"/A_proc2_octree.csv",DataFrame)
# octree_M2=Array(sparse(Mcsv[!,:I],Mcsv[!,:J],Mcsv[!,:V]))

# octree_M1 == octree_M2
# octree_M1 ≈ octree_M2
# norm(octree_M1-octree_M2)


# ## RHS vector

# b1 = CSV.read("b1.csv", DataFrame)
# b2 = CSV.read("b2.csv", DataFrame)
# norm(b1.col-b2.col)

# Mcsv=CSV.read(dir*"/b_proc1.csv",DataFrame)
# b1=Mcsv.col

# Mcsv=CSV.read(dir*"/b_proc2.csv",DataFrame)
# b2=Mcsv.col

# b1 == b2

# norm(b1-b2)

# Mcsv=CSV.read(dir*"/b_proc1_octree.csv",DataFrame)
# octree_b1=Mcsv.col

# Mcsv=CSV.read(dir*"/b_proc2_octree.csv",DataFrame)
# octree_b2=Mcsv.col

# octree_b1 == octree_b2
# octree_b1 ≈ octree_b2
# norm(octree_b1-octree_b2)


# p_fe = 1, slope = 3.776689394494557
# 	 Error: [2.75456060042055e-6, 4.157379495931199e-5, 0.0006132366256157495, 0.006919990522519908]
# 	    nc: [256.0, 64.0, 16.0, 4.0]
# 	    dx: [0.09045015681978345, 0.1809003136395669, 0.3618006272791338, 0.7236012545582676]
