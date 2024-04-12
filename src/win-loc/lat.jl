using Dates

function merge_vectors(vector1::Vector{Any}, vector2::Vector{Pair{String, DateTime}})
    merged_vector = []
    for (item1, item2) in zip(vector1, vector2)
        merged_item = (item1, item2.first, item2.second)
        push!(merged_vector, merged_item)
    end
    return merged_vector
end

function merge_vectors(vector1::Vector{String}, vector2::Dict{String, DateTime})
    merged_vector = []
    for (item1, item2) in zip(vector1, vector2)
        merged_item = (item1, item2.first, item2.second)
        push!(merged_vector, merged_item)
    end
    return merged_vector
end


function latx()
    """
    list_files_sorted_by_last_change
    """
    directory = pwd()
    files = readdir(directory)  
    file_times = Dict{String, Dates.DateTime}()
    for file in files
        file_path = joinpath(directory, file)
        stat_info = stat(file_path)
        file_times[file] = Dates.unix2datetime(stat_info.mtime)
    end
    #xf = sort(file_times, by = x -> x[2], rev=true)
    file_times = collect(file_times)
    xf = sort(file_times, by = x -> x[2], rev=true)
    xf = first(xf,11)
    sz = []                 #::Vector{Float64}
    for i in 1:length(xf)
       push!(sz,round(stat(xf[i][1]).size/1024^2,digits=4))
    end
    xf2 = merge_vectors(sz,xf)
    for (size, file, dt) in xf2
        datetime = Dates.format(dt, "yyyy-mm-dd HH:MM")
        printstyled(rpad(file,35),color=:yellow),
        printstyled("$datetime\t" ,color=:green),
        printstyled(rpad(size,7)," MB\n",color=:yellow)    
    end
end

latx()

# function latn()
#     """
#     list_files_sorted_by_last_change
#     """
#     directory = pwd()
#     files = readdir(directory)  
#     file_times = Dict{String, Dates.DateTime}()
#     for file in files
#         file_path = joinpath(directory, file)
#         stat_info = stat(file_path)
#         file_times[file] = Dates.unix2datetime(stat_info.mtime)
#     end
    
#     #f2 = merge_vectors(files,file_times)
#     sorted_files = sort(files, by=file -> file_times[file], rev=true)
#     #sorted_files = sort(f2, by=f2 -> file_times[f2], rev=true)
    
#     # #xf = sort(file_times, by = x -> x[2], rev=true)
#     # xf = sort(file_times, by = x -> x[2][1], rev=true)
#     # xf = first(xf,11)
#     # sz = []                 #::Vector{Float64}
#     # for i in 1:length(xf)
#     #    push!(sz,round(stat(xf[i][1]).size/1024^2,digits=4))
#     # end
#     # xf2 = merge_vectors(sz,xf)
#     # for (size, file, datetime) in xf2
#     #     printstyled(rpad(file,35),color=:yellow),printstyled(
#     #         "$datetime  $size MB\n",color=:green)
#     # end
#     xf = first(sorted_files,11)
#     printstyled("$xf\n",color=:green)
# end

#latn()


#julia $pypt/julia/lat.jl
#julia --startup-file=no C:\Users\Public\Documents\Python_Scripts\julia\lat.jl
#jl "C:/Users/Public/Documents/Python_Scripts/julia/lat.jl"

# function pretty_print_files(files::Vector{Pair{String, Dates.DateTime}})
#     for (file, datetime) in files
#         @printf("%-35s %s\n", file, datetime)
#     end
# end
# pretty_print_files(xf)


# function latraw()
#     """
#     list_files_sorted_by_last_change
#     """
#     directory = pwd()
#     files = readdir(directory)  
#     file_times = Dict{String, Dates.DateTime}()
#     for file in files
#         file_path = joinpath(directory, file)
#         stat_info = stat(file_path)
#         file_times[file] = Dates.unix2datetime(stat_info.mtime)
#     end
#     #sorted_files = sort(files, by=file -> file_times[file], rev=true)
        

#     #printstyled("$sorted_files  $file_times[file]\n",color=:magenta)
#     #xf = first((sorted_files,file_times),5)
#     # xf = first(file_times,5)
#     #xf = sort(xf, by=xf -> file_times[file], rev=true)

#     #xf = first(sorted_files,5)
#     # xc = []
#     # for i in xf 
#     #     v = Dates.unix2datetime(stat(xf[i]).mtime)
#     #     push!(v,xc)
#     # end

#     # Dates.unix2datetime(stat(xf[2]).mtime)
#     xf = sort(file_times, by = x -> x[2],rev=true)
#     xf = first(xf,11)

#     sz = []                 #::Vector{Float64}
#     for i in 1:length(xf)
#        push!(sz,round(stat(xf[i][1]).size/1024^2,digits=3))
#     end

#     #hcat(xf,sz)


#     function merge_vectors(vector1::Vector{Float64}, vector2::Vector{Pair{String, DateTime}})
#         merged_vector = []
#         for (item1, item2) in zip(vector1, vector2)
#             merged_item = (item1, item2.first, item2.second)
#             push!(merged_vector, merged_item)
#         end
#         return merged_vector
#     end
    
#     xf2 = merge_vectors(sz,xf)

#     for (size, file, datetime) in xf2
#         printstyled(rpad(file,35),color=:yellow),printstyled(
#             "$datetime  $size MB\n",color=:green)
#     end
#     #printstyled("$xf\n",color=:magenta)
#     #typeof(xf)
#     #return DataFrame(name=sorted_files,last_modified=file_times[file])
# end

