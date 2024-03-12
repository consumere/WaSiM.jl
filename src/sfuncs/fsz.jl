# 
function fsz(rootdir::String=".";rec=false)
    total_size = 0
    files = readdir(rootdir)  # Get a list of files in the current directory
    
    for file in files
        filepath = joinpath(rootdir, file)  # Get the full path of the file
        if isfile(filepath)
            size = stat(filepath).size  # Get the file size
            total_size += size
        end
    end
    
    nr=length(files)

    total_size_mb = total_size / (1024 * 1024)  # Convert size to megabytes
    total_size_mb = round(total_size_mb, digits=2)  # Round to 2 decimal places
    printstyled("Total size of $nr files in $(pwd()): $total_size_mb MB\n", color=:green)
    
    if rec
        dirs = readdir(rootdir)
        for dir in dirs
            if isdir(dir)
                cd(dir)
                fsz()
                cd("..")
            end
        end
    end
    #return total_size_mb
end

