using Printf
# using Glob
# #using Statistics #for sum
# fl = glob("**/*", ".") # get all files recursively from the current directory
# #fsum = []
# #local fsum = []

function sdrec()
cwd="."
global fsum = 0
global counter = 0 
for (root, dirs, files) in walkdir(cwd)
   for f in files
    if isfile(f)
        counter += 1
        fsum += filesize(f)
    end
#    return(counter,fsum)
   end
end
println(counter, " files from\t", pwd()) # print the number of files and the current directory
println("total size is ", @sprintf("%.2f", fsum / (1024^2)), " MB ") # print the total size in megabytes
end

sdrec()

# fsum = 0;
# for f in fl
#     if isfile(f)
#     local fsum += filesize(f) # add the file size to the sum
# end


# function jdd()
#     cwd = pwd()
#     dirs = readdir(".")
#     out = []
#     for dir in dirs
#         if isdir(dir)
#             size = 0
#             for (root, dirs, files) in walkdir(dir)
#                 for file in files
#                     size += stat(joinpath(root, file)).size
#                 end
#             end
# 	    @printf("%-40s %15.2f MB\n","$(cwd)\\$dir:",size/1024^2);
#         push!(out,size/1024^2)
#         end
#     end
#     return(map(string,out))
# end

# function dd()
#     cwd = pwd()
#     dirs = readdir(".")
#     osize = 0
#     out = []
#     for (root, dirs, files) in walkdir(cwd)
#      for file in files
#          osize += stat(joinpath(root, file)).size
#      end
#     end 
#     @printf("%-40s %15.2f GB\n","$(cwd):",osize/1024^3);
#     push!(out,osize/1024^2)
#     return(map(string,out))
# end 

# println("folder sizes on julia...:\n")

# dd()

# jdd()
