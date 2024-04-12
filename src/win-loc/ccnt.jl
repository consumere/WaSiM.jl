# 	#faster than C
# jlcnt(){
#      x="C:\Users\Public\Documents\Python_Scripts\julia\win\ccnt.jl"
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
    if args[1]==raw".."
        path=dirname(pwd())
    end
    #pthome=pwd()
    #path = length(args) > 0 ? args[1] : pthome  # Directory path
    
    n, s = count_files(path, 0)
    
    printstyled("Directory: $path\n",color=:yellow)
    printstyled("Number of files: $n\n",color=:yellow)
    printstyled("Total size: $(round(s / (1024 * 1024), digits=2)) MB\n",color=:yellow)
end

# Run the main function
main(ARGS)

