
function fd()
    cwd = pwd()
    dirs = readdir(cwd)
    printstyled("$cwd\n",color=:magenta)
    for dir in dirs
        if isdir(dir)
            printstyled("$dir\n",color=:green)
        end
    end
end

# check if command line arguments are provided
if length(ARGS) == 0
    println("Please provide the path to the directory to be sorted.")
    println("Usage: julia fds.jl /path/to/directory")
    fd()
    exit(1)
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
    os=round(osize/1024^3,digits=3)
    printstyled(rpad("$(cwd): $os GB\n",20),lpad("($cnt files)\n",10),color=:green);
end 



function print_sorted_sizes(dir::String)
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

function get_folder_size(folder::String)
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

# get the directory from the command line argument
dir = ARGS[1]
#printstyled("folders smaller 1MB will be omitted...\n",color=:red)
printstyled("total folder size:\n",color=:red)
dd()
print_sorted_sizes(dir)
