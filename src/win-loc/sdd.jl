using Printf
#@printf("%.3f", rand())

function jdd()
    cwd = pwd()
#    cwd = cd("../")
    dirs = readdir(".")
    for dir in dirs
        if isdir(dir)
#	println("check:",dir)	
            size = 0
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    size += stat(joinpath(root, file)).size
                end
		#sort!(files, by = x -> mtime(x), rev=true)
            end
#            println("$(cwd)/$dir: \t", size / 1024^2, "\t MB")
#	    s = @sprintf "this is a %s %15.1f" "test" 34.567;
	    #@printf("%s %15.1f MB\n","$(cwd)/$dir:",size/1024^2);
	    @printf("%-30s %15.2f MB\n","---/$dir:",size/1024^2);
        end

    end
end

function dd()
    cwd = pwd()
    osize = 0
    for (root, dirs, files) in walkdir(cwd)
     for file in files
         osize += stat(joinpath(root, file)).size
     end
    end 
    @printf("%-30s %15.3f GB\n",cwd,osize/1024^3);
end 

#println("time sorted folders...")
#println("folder sizes on julia...\nlocation:\n",pwd())
println("folder sizes on julia...\nlocation:")
dd()
jdd()


# using Dates
# function jyll(pattern::String)
    # files = filter(x -> occursin(pattern, x), readdir("."))
    # sort!(files, by = x -> mtime(x), rev=true)
    # for file in files[1:min(10, length(files))]
        # size = stat(file).size
        # time = stat(file).mtime
        # println("$(lpad(file, 40)): \t", (size / 1024^2), "\t MB")
    # end
# end


# using FileIO
# dirs = filter(isdir, readdir("."))
# println(dirs)

#using Pkg
#Pkg.add("FileIO")
#using FileIO
#using Base

# Find all directories under the current directory
#dirs = filter(isdir, readdir("."))


# # Iterate over the directories and print their size if it's greater than 10 MB
# for dir in dirs
    # path = joinpath(".", dir)
    # f = stat(path)
    # println(dirs)
    # if f.size >= 10^7
        # size = f.size
        # scale = 0
        # while size > 1024
            # size /= 1024
            # scale += 1
        # end
        # units = ["B", "KB", "MB", "GB", "TB", "PB"][scale + 1]
        # println(rpad(round(size, 3), 6, " "), " ", units, " ", f.ctime, "  ", dir)
    # end
# end
