#--startup-file=no

# using Distributed
# @distributed
using Base.Threads
function count_files(path::String, level::Int64)
    n = 0  # Number of files
    s = 0  # Total size in bytes
    
    Threads.@threads for entry in readdir(path)
        if entry == "." || entry == ".."
            continue
        end
        
        subpath = joinpath(path, entry)
        stat_result = stat(subpath)
        
        if isfile(stat_result)
            n += 1
            s += stat_result.size
        elseif isdir(stat_result)
            println("$(repeat(" ", level * 2))[$entry]")
            subn, subs = count_files(subpath, level + 1)
            n += subn
            s += subs
        end
    end
    
    return n, s
end

function main(args)
    path = length(args) > 0 ? args[1] : raw"."  # Directory path
    
    n, s = count_files(path, 0)
    
    println("Number of files: $n")
    println("Total size: $(round(s / (1024 * 1024), digits=2)) MB")
end

# Run the main function
main(ARGS)


# global cntf = 0
# global tot = 0
# wd = pwd()

# function process_directory(root)
#     local cnt = 0
#     local sz = 0.0
#     s = filter(isfile, readdir(root))
#     for file in s
#         cnt += 1
#         sz += filesize(joinpath(root, file)) / 2^20
#     end
#     return cnt, sz
# end

# function parallel_sum_dirs()
#     dirs = filter(isdir, readdir(wd))
#     thread_counts = Vector{Int}(undef, length(dirs))
#     thread_sizes = Vector{Float64}(undef, length(dirs))
#     @distributed for (i, d) in enumerate(dirs)
#         thread_counts[i], thread_sizes[i] = process_directory(joinpath(wd, d))
#     end
#     return thread_counts, thread_sizes
# end
# cntf, tot = parallel_sum_dirs()
# printstyled("Size of $wd and its subdirectories:\n $cntf files\t$tot MB\n", color=:green)



