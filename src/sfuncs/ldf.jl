# 
function ldf(path::AbstractString, prefix::AbstractString)
        files = readdir(path)
        dfs = DataFrame[]
        for file in files
            if isfile(file) && occursin(Regex(prefix),file)&& (!occursin(r"txt|yrly|nc|png|svg",file))
            file_path = joinpath(path, file)
        println("reading",file_path)
        p1 = loaddf(file_path)
        push!(dfs, p1)
            end
        end
        return(dfs)
    end

    """
    r = cor(simulated, observed)
    α = std(simulated) / std(observed)
    β = mean(simulated) / mean(observed)
    KGE = 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
    """
    