# 
function fdd(;cwd=pwd())
    dirs = readdir(cwd)
    
    if length(dirs) == 0 
        println("$cwd is empty!")
        return
    end
    
    if filter(x -> (isdir(x)),dirs) == []
        bn = basename(cwd)
        @info "no dirs in $bn !"
        dd()
        return
    end
    
    s = []
    for dir in dirs
        if isdir(dir)
            push!(s,joinpath(cwd, dir))
            size = 0
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    size += stat(joinpath(root, file)).size
                end
            end
	    @printf("%-40s %15.2f MB\n","$(cwd)\\$dir:",size/1024^2);
	end
    end
    return(s)
end


