# 
function du(;cwd=pwd())
    osize = 0
    n = 0
    for (root, dirs, files) in walkdir(cwd)
     for file in files
         osize += stat(joinpath(root, file)).size
         n += 1
     end
     for dir in dirs
        printstyled("check dir: $dir\n",color=:light_red)
     end
    end 
    println("$(n) files in directory")
    @printf("%-40s %15.2f MB\n","$(cwd):",osize/1024^2)
end 

