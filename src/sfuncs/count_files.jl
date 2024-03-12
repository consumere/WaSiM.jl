# 
function count_files(path::String, level::Int64)
    n = 0  # Number of files
    s = 0  # Total size in bytes
    
    for entry in readdir(path)
        if entry == "." || entry == ".."
            continue
        end
        
        subpath = joinpath(path, entry)
        stat_result = stat(subpath)
        
        if isfile(stat_result)
            n += 1
            s += stat_result.size
        elseif isdir(stat_result)
            println("$(repeat(" ", level * 2))[$entry]")
            subn, subs = count_files(subpath, level + 1)
            n += subn
            s += subs
        end
    end
    
    return n, s
end

