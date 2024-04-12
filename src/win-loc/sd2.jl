#no printf!

# function sdf()
    
#     function calculate_folder_size(directory)
#         size = 0
#         count = 0
#         for (root, dirs, files) in walkdir(directory)
#             for file in files
#                 size += stat(joinpath(root, file)).size
#                 count += 1
#             end
#         end
#         return size, count
#     end

#     function print_folder_size(directory, size, count)
#         size_gb = round(size / 1024^3, digits=3)
#         printstyled(rpad("$directory: $size_gb GB", 40), color=:green)
#         printstyled(lpad("($count files)", 20), "\n", color=:green)
#     end

#     printstyled("folder sizes on Julia...\n", color=:red)

#     cwd = pwd()
#     dirs = readdir(cwd)

#     rs, rcnt = calculate_folder_size(cwd)

#     print_folder_size(cwd,rs,rcnt)

#     n = repeat(" - -", 10)
#     printstyled(n * " subfolders of " * basename(cwd) * n, "\n", color=:yellow)

#     for dir in dirs
#         if isdir(dir)
#             size, count = calculate_folder_size(joinpath(cwd, dir))
#             print_folder_size(dir, size, count)
#         end
#     end
# end

#This version eliminates the need to traverse the directory structure twice for the current directory and its subdirectories. Instead, it calculates the size and count on-the-fly during a single traversal. This should result in improved performance, especially for larger directory structures.
function sdf()
    function calculate_folder_size_recursive(directory)
        size = 0
        count = 0
        for (root, dirs, files) in walkdir(directory)
            for file in files
                size += stat(joinpath(root, file)).size
                count += 1
            end
        end
        return size, count
    end

    function print_folder_size(directory, size, count)
        size_gb = round(size / 1024^3, digits=3)
        printstyled(rpad("$directory: $size_gb GB", 40), color=:green)
        printstyled(lpad("($count files)", 20), "\n", color=:green)
    end

    printstyled("folder sizes on Julia...\n", color=:red)

    cwd = pwd()
    rs, rcnt = calculate_folder_size_recursive(cwd)

    print_folder_size(cwd, rs, rcnt)

    n = repeat(" - -", 10)
    printstyled(n * " subfolders of " * basename(cwd) * n, "\n", color=:yellow)

    for (dirpath, dirs, files) in walkdir(cwd)
        for dir in dirs
            size, count = calculate_folder_size_recursive(joinpath(dirpath, dir))
            print_folder_size(joinpath(dirpath, dir), size, count)
        end
    end
end


sdf()

# function dd()
#     cwd = pwd()
#     dirs = readdir(".")
#     osize = 0
#     cnt = 0
#     out = []
#     for (root, dirs, files) in walkdir(cwd)
#      for file in files
#          osize += stat(joinpath(root, file)).size
#          cnt += 1
#      end
#     end 
#     #@printf("%-40s (%d files) %15.4f GB\n","$(cwd):",cnt,osize/1024^3);
#     os=round(osize/1024^3,digits=3)
#     printstyled(rpad("$(cwd): $os GB\n",20),lpad("($cnt files)\n",10),color=:green);
# end 

# function jdd()
#     cwd = pwd()
#     dirs = readdir(".")
#     for dir in dirs
#         if isdir(dir)
#             osize = 0
#             cnt = 0
#             for (root, dirs, files) in walkdir(dir)
#                 for file in files
#                     osize += stat(joinpath(root, file)).size
#                     cnt += 1
#                 end
#             end
# 	   # @printf("%-40s (%d files) %15.2f MB\n","$(cwd)\\$dir:",cnt,size/1024^2);
#         os=round(osize/1024^2,digits=3)
#         printstyled(rpad("$(dir): ",40),lpad("$os MB", 10),lpad("($cnt files)\n",30),color=:green);
#         end
#     end
# end

# printstyled("folder sizes on julia...\n",color=:red)
# # println("folder sizes on julia...:\n",pwd(),":")

# dd()
# n=repeat(" - -",10)
# #println("\n",n*" subfolders of ",basename(pwd())*n,"\n")
# printstyled(n*" subfolders of ",basename(pwd())*n,"\n",color=:yellow)

# jdd()



# cd(dirname(r))
# pwd()

# sdf()
# cd("..")
# sdf()
