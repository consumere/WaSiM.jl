# 
function jdd()
    cwd = pwd()
    dirs = readdir(".")
    for dir in dirs
        if isdir(dir)
            size = 0
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    size += stat(joinpath(root, file)).size
                end
            end
	    @printf("%-40s %15.2f MB\n","$(cwd)\\$dir:",size/1024^2);
        end
    end
end

