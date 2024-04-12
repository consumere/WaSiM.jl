using Printf
n=repeat(" - -",10)

#s="\n",n*" subfolders of ",basename(pwd())*n,"\n"
#printf "%-30s %s %30s\n" a b c


function jdd()
    cwd = pwd()
    dirs = readdir(".")
    #out = []
    for dir in dirs
        if isdir(dir)
            size = 0
            cnt = 0
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    size += stat(joinpath(root, file)).size
                    cnt += 1
                end
            end
	    @printf("%-40s (%d files) %15.2f MB\n","$(cwd)\\$dir:",cnt,size/1024^2);
	    #@printf("%-40s %.2f MB\t(%d files)\n","$(cwd)\\$dir:",size/1024^2,cnt);
        #push!(out,size/1024^2)
        end
    end
    #return(map(string,out))
end

function dd()
    cwd = pwd()
    dirs = readdir(".")
    osize = 0
    cnt = 0
    out = []
    for (root, dirs, files) in walkdir(cwd)
     for file in files
         osize += stat(joinpath(root, file)).size
         cnt += 1
     end
    end 
    @printf("%-40s (%d files) %15.4f GB\n","$(cwd):",cnt,osize/1024^3);
    #@printf("%-40s %.4f GB\t(%d files)\n","$(cwd):",osize/1024^3,cnt);
    #push!(out,osize/1024^2)
    #return(map(string,out))
end 

#println("folder sizes on julia...:\n",pwd(),":")
@printf("Folder sizes on Julia...:\n")
dd()
# n=repeat(" - -",10)
# println("\n",n*" subfolders of ",basename(pwd())*n,"\n")
#println(s)
@printf("\n%-30s %s%s %30s\n\n",n," subfolders of ",basename(pwd()),n)
jdd()

#v = [dd(),"MB",jdd(),"MB"]
#v = [dd(),jdd()]

#s = reduce(vcat, v)
#println(s)

#sort!(v, by = x -> x[2]); v
#broadcast(string,v)