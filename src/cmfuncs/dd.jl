# 
function dd()
        cwd = pwd()
        osize = 0
        for (root, dirs, files) in walkdir(cwd)
        for file in files
            osize += stat(joinpath(root, file)).size
        end
        end 
        @printf("%-40s %15.3f GB\n","$(cwd):",osize/1024^3);
    end 

    