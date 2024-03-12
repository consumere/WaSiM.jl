# 
function fdd()
        cwd = pwd()
        dirs = readdir(".")
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
    #        else
    #	    @printf("%-40s\n","$(cwd)");
        end
        end
        return(s)
    end

    