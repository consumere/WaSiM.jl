function print_sorted_sizes(dir)
    folders = [joinpath(dir, f) for f in readdir(dir)]
    sizes = Dict()
    for f in folders
        if isdir(f)
            sizes[f] = get_folder_size(f)
        end
    end
    sorted_folders = sort(collect(keys(sizes)), by=x->sizes[x], rev=false)
    for f in sorted_folders
        if sizes[f] >= 1e9
            printstyled(rpad(f,60, ' '), rpad(sizes[f] ÷ 1024^3, 6, ' '), "GB\n",color=:yellow)
        elseif sizes[f] >= 1e6
            printstyled(rpad(f,60, ' '), rpad(sizes[f] ÷ 1024^2, 6, ' '), "MB\n",color=:green)
        elseif sizes[f] >= 1e3
            printstyled(rpad(f,60, ' '), rpad(sizes[f] ÷ 1024^1, 6, ' '), "KB\n",color=:magenta)
        end
    end
end

#printstyled("KB\n",color=:magenta)
# 	println("$(f): $(sizes[f] / 1_000_000) MB")

function get_folder_size(folder)
    files = readdir(folder)
    size = 0
    for file in files
        path = joinpath(folder, file)
        if isfile(path)
            size += stat(path).size
        elseif isdir(path)
            size += get_folder_size(path)
        end
    end
    return size
end

printstyled("folders smaller 1MB will be omitted...\n",color=:red)
print_sorted_sizes(pwd())

# # Round the value to two decimal places
# value = round(value; digits=2)