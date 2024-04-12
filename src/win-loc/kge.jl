##to remove artifacts
using Pkg
using Dates
using PkgCleanup
PkgCleanup.artifacts() 
PkgCleanup.manifests() 
Pkg.gc(;collect_delay=Dates.Day(0))

using Rasters      #geht zzt in ubuntu
using ARC

using Base.Threads
nthreads()

##geospatial julia
#https://github.com/acgeospatial/Julia_Geospatial/tree/master/02_Notebooks

using ArchGDAL; const AG = ArchGDAL
f="D:\\Bodendaten\\Buek200\\bk200.shp"
ls=AG.listdrivers()
using DataFrames
dd = DataFrame(ls)
dd[("ARG")]
ls[("R")]
#dd[In("ARG")]
d=AG.read(f)
layer = ArchGDAL.getlayer(d, 0)

d=AG.load(f)

using Printf
using Statistics
function kge(obs::Array{Float64,1}, pred::Array{Float64,1}, obs_mean::Float64, obs_std::Float64, pred_mean::Float64, pred_std::Float64)::Float64
    r = cor(obs, pred)
    alpha = pred_std / obs_std
    beta = pred_mean - alpha * obs_mean
    kge_val = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
    return kge_val
end

observed = [1.2, 2.3, 3.4, 4.5, 5.6]
predicted = [1.1, 2.2, 3.3, 4.6, 5.7]
obs_mean = mean(observed)
obs_std = std(observed)
pred_mean = mean(predicted)
pred_std = std(predicted)

kge_val = kge(observed, predicted, obs_mean, obs_std, pred_mean, pred_std)
println("KGE value: ", kge_val)

function kge2(observed, simulated)
    r = cor(observed, simulated)
    α = std(simulated) / std(observed)
    β = mean(simulated) / mean(observed)
    return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
end


observed = [1.0, 2.0, 3.0, 4.0, 5.0]
simulated = [0.9, 1.8, 3.1, 4.2, 5.3]
kge_value = kge2(observed, simulated)
println("KGE value is $kge_value")

##3rd try

using CSV, Statistics
using DataFrames; 
#CSV.read(source, DataFrame)
#observed = CSV.read("observed_data.csv")[:, 1]
#simulated = CSV.read("simulated_data.csv")[:, 1]
infile="D:/Wasim/Goldbach/revision/v3/qout_gb"

observed = CSV.read(infile,DataFrame)[:,5]
simulated = CSV.read(infile,DataFrame)[:,6]

function kge(observed::Vector{T}, simulated::Vector{T}) where T<:AbstractFloat
    # Calculate the mean, standard deviation, and correlation coefficient of the observed and simulated data
    obs_mean, sim_mean = mean(observed), mean(simulated)
    obs_std, sim_std = std(observed), std(simulated)
    r = cor(observed, simulated)
    
    # Calculate the KGE using the formula
    sqrt((r - 1)^2 + (obs_std / sim_std - 1)^2 + (obs_mean / sim_mean - 1)^2)
end

kge_value = kge(observed, simulated)
println("KGE value is $kge_value")

kge2_value = kge2(observed, simulated)
println("KGE value is $kge2_value")


using Statistics

# Calculate means and standard deviations
observed_mean = mean(observed)
predicted_mean = mean(simulated)
observed_std = std(observed)
predicted_std = std(simulated)

# Calculate KGE
xkge = 1 - sqrt((observed_std - predicted_std)^2 + (observed_mean - predicted_mean)^2) / (2 * observed_std^2)
println("KGE value is $xkge")



#################27.02.23#########
/home/ubu/.julia
#file hist
npp ./logs/repl_history.jl 

#Microsoft.PowerShell.Core\FileSystem::\\wsl.localhost\Ubuntu-20.04\home\ubu\.julia:
#/calculating folder size - starting on: /home/ubu/.julia                                                       
#/home/ubu/.julia       7.433  GB 

# mode: julia
using OhMyREPL
# time: 2023-02-27 11:47:52 CET
# mode: julia
	using Rasters
	ts=Raster("tsoilrcm_stack.2017.nc",missingval=-9999)
	using Plots; ts[t=2]|>plot



    ##geht aber NUR in ubu20 wg der registry...
/mnt/d/wslconda/rainfarm/julia-1.8.5/bin/julia --threads auto -q
using Base.Threads
nthreads()

cd("/mnt/d/temp/saale/out_smf200/smf-v1")
using Rasters,Plots
ts=Raster("tsoilrcm_stack.2017.nc",missingval=-9999)
ts[t=2]|>plot

using Rasters,Plots
pwd()
ts=Raster("tsoilrcm_stack.2017.nc",missingval=-9999)
ts[t=2]|>plot

using BenchmarkTools 

init() = begin
    n = 250
    a1 = Vector{Vector{Float64}}(undef, n)
    for i in 1:n 
        a1[i] = [rand(), rand()]
    end
    a1
end

# without multithreading
round_serial(a1) = begin
    n = length(a1)
    a2 = Vector{Vector{Float64}}(undef, n)
    for i in 1:n 
        a2[i] = [round(x, digits=5) for x in a1[i]]
    end
    a2
end

# with multithreading
round_parallel(a1) = begin
    n = length(a1)
    a3 = Vector{Vector{Float64}}(undef, n)
    Threads.@threads for i in 1:n 
        a3[i] = [round(x, digits=5) for x in a1[i]]
    end
    a3
end

a1 = init()
a2 = @btime round_serial($a1)
a3 = @btime round_parallel($a1)
@assert isapprox(a3, a2)