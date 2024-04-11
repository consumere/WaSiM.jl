# 
function rglob(prefix::Regex)
        rootdir="."
        results = []
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if (occursin(prefix,filename))
                    push!(results, joinpath(looproot, filename)) 
                end
            end
        end
        return string.(results)
    end


    #go dir up
    