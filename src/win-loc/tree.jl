#tree.jl



#using FileIO #no need
# function jfind()
#     cwd = pwd()
#     dirs = readdir(".")
#     prefix = "   "
#     for dir in dirs
#         if isdir(dir)
#             relpath = relative_path(cwd, joinpath(cwd, dir))
#       #      relpath = replace(relpath, "mnt/", "")
#             relpath = uppercasefirst(relpath)
#             relpath = replace(relpath, "/", ":")
#             relpath = replace(relpath, r"[^-][^/]*/", "--")
#             relpath = prefix * "|" * relpath
#             if relpath == prefix * "|."
#                 relpath = "$(uppercasefirst(cwd))\n$relpath"
#             end
#             println("$dir:\n$relpath")
#         end
#     end
# end
#cwd=pwd()

function fdi(cwd = ".", prefix=" ")
    res::Vector{String} = []
    for (looproot, dir, filenames) in walkdir(cwd)
        for relpath in dir
            push!(res, joinpath(looproot,relpath))
        end
    end
    return res
end

#show(fdi())
paths = fdi()

#    println("generating\t",basename(z),"...")
#    m=(filter(line -> occursin(r"^[LIN. R]|^[LOG. R]|^CO",line),readlines(open(z))))
for relpath in paths
    prefix = " "
            #relpath = uppercasefirst(relpath)
            relpath = replace(relpath, r"[^\\]"=> ":",count=1)
            #relpath = replace(relpath, r"[^\\]*."=>"--",count=1)
            
            ###helper to fast cdinto....
            relpath = replace(relpath, r"[^\\]*."=>"cd ",count=1)
            relpath = replace(relpath, "\\"=> "/")
            #relpath = replace(relpath, "/"=> ":")
            #relpath = replace(relpath, r"[^-][^/]*/"=>"--")
            #relpath = replace(relpath, r"[^-][^\\]*\\"=>"--")
            # relpath = prefix * "|" * relpath
            # if relpath == prefix * "|."
            #         relpath = "$(uppercasefirst(cwd))\n$relpath"
            #     end
    println(relpath)
end

# println("$dir:\n$relpath")
# end

# results = []
# for (looproot, dir, filenames) in walkdir(cwd)
#         push!(results, joinpath(dir,filenames))
# end
        
# isdir(cwd)


# splitdir(paths[5])
# x=paths[5]

# for (looproot, dir, filenames) in walkdir(x)
#     for i in filenames
#         if occursin(r"^r.*", i) println(i) end
#     end
# end


# open(x) do f
#     for i in eachline(f)
#         if occursin(r"^r.*", i) println(i) end
#     end
# end

# #das geht niur unter ubu20, wg netcdf package...
# using NetCDF
# for (looproot, dir, filenames) in walkdir(x)
#     for i in filenames
#         if occursin(r"^r.*", i) && (endswith(i,"nc")) ncinfo(i) end
#     end
# end

# for (looproot, dir, filenames) in walkdir(x)
#     for i in filenames
#         if occursin(r"^r.*.nc", i) ncinfo(i) end
#     end
# end


# l=[]
# for (looproot, dir, filenames) in walkdir(x)
#     for i in filenames
#         if occursin(r"^sb0.*.nc", i) 
#             ncinfo(i) 
#             push!(l, ncread(i))     #hier fehlt die var... 
#         end
#     end
# end


# m=(filter(line -> occursin(r"^Var",line),ncinfo(t)))
