function print_sorted_sizes(rootdir::String)
    sizes = Dict{String,Int}()
    for entry in readdir(rootdir)
        if entry != "." && entry != ".." # skip current and parent directory links
            fullpath = joinpath(rootdir, entry)
            if isfile(fullpath)
                sizes[fullpath] = stat(fullpath).size
            elseif isdir(fullpath)
                for (subdir, _, filenames) in walkdir(fullpath)
                    for filename in filenames
                        fullpath = joinpath(subdir, filename)
                        sizes[fullpath] = stat(fullpath).size
                    end
                end
            end
        end
    end
    sorted_sizes = sort(collect(sizes), by=x->x[2], rev=false)
    for (name, size) in sorted_sizes
        if size >= 1000000
            printstyled(rpad(name, 60, ' '), rpad(size ÷ 10^6, 10, ' '), "MB\n",color=:green)
        end
#        println(rpad(name, 60, ' '), rpad(size ÷ 10^6, 10, ' '), "MB")
    end
end

#        printstyled("$(f):\t $(sizes[f] / 1_000_000) MB\n",color=:green)

printstyled("folders smaller 1MB will be omitted...\n",color=:red)
print_sorted_sizes(".")
