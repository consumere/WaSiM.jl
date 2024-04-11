# 
function du()
        cwd = pwd()
        n = length(readdir(cwd))
    #    dirs = readdir(cwd)
        osize = 0
        for (root, dirs, files) in walkdir(cwd)
        for file in files
            osize += stat(joinpath(root, file)).size
        end
        end 
        println("$(n) files in directory")
        @printf("%-40s %15.2f MB\n","$(cwd):",osize/1024^2)
    end 

    