using DataFrames, CSV, Statistics
#CSV.read(source, DataFrame)

function read_f(path::AbstractString, ext::AbstractString)
    files = readdir(path)
    dfs = DataFrame[]
    for file in files
        file_path = joinpath(path, file)
        if isfile(file_path) && endswith(file, ext)
            df = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
            push!(dfs, df)
        elseif isdir(file_path)
            dfs_in_subdir = read_f(file_path, ext)
            dfs = vcat(dfs, dfs_in_subdir)
        end
    end
    return dfs
end


#directory_path = ARGS[1]
directory_path = "D:/Wasim/Goldbach/revision/v3"
directory_path = "D:/Wasim/Goldbach/revision/"
#dfs = read_f(directory_path,"qout_gb")
dd = read_f(directory_path,"qout_gb")
size(dd)
#dfs = read_files(directory_path)
#observed = CSV.read("observed_data.csv")[:, 1]
#simulated = CSV.read("simulated_data.csv")[:, 1]
#infile="D:/Wasim/Goldbach/revision/v3/qout_gb"
#infile=ARGS[1]; #geht
#infile=ARGS[1]; #geht

#observed = CSV.read(infile,DataFrame)[:,5]
#simulated = CSV.read(infile,DataFrame)[:,6]

observed  = dd[:,5]
simulated = dd[:,6]



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