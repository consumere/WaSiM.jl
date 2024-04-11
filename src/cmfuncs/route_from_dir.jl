# 
function route_from_dir(dir::String)
        dirs = readdir(dir)
        routes::Vector{String} = []
        for directory in dirs
            if isfile("$dir/" * directory)
                push!(routes, "$dir/$directory")
            else
                if ~(directory in routes)
                    newread = dir * "/$directory"
                    newrs = route_from_dir(newread)
                    [push!(routes, r) for r in newrs]
                end
            end
        end
        routes
    end

    