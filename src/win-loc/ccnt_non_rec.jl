

###nonrec
# function jlcntnr(){
# 	#faster than C
#      x="C:\Users\Public\Documents\Python_Scripts\julia\win\ccnt_non_rec.jl"
#      jl $x
# }

function count_files(path::String, level::Int64)
    n = 0  # Number of files
    s = 0  # Total size in bytes
    
    for entry in readdir(path)
        if entry == "." || entry == ".."
            continue
        end
        
        subpath = joinpath(path, entry)
        stat_result = stat(subpath)
        
        if isfile(stat_result)
            n += 1
            s += stat_result.size
        elseif isdir(stat_result)
            continue
            # println("$(repeat(" ", level * 2))[$entry]")
            # subn, subs = count_files(subpath, level + 1)
            # n += subn
            # s += subs
        end
    end
    
    return n, s
end

function main(args)
    path = length(args) > 0 ? args[1] : raw"."  # Directory path
    
    n, s = count_files(path, 0)
    
    pathname = length(args) > 0 ? realpath(args[1]) : pwd()

    println("Non-recursive count on: $pathname")
    println("Number of files: $n")
    println("Total size: $(round(s / (1024 * 1024), digits=2)) MB")
end

# Run the main function
main(ARGS)