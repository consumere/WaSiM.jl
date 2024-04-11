# 
function rtreeo()
    dir = pwd()
    dirs = readdir(dir)
    routes::Vector{String} = []
    for directory in dirs
        if isdir("$dir/" * directory)
            push!(routes, "$dir/$directory")
        else
            if ~(directory in routes)
                push!(routes, "$routes/$directory")
                #newread = dir * "/$directory"
                #push!(newread, "$dir/$directory")
                #[push!(routes, r) for r in newrs]
            end
        end
    end
    routes
end

