# 
function print_folder_size(directory, size, count)
        size_gb = round(size / 1024^3, digits=3)
        printstyled(rpad("$directory: $size_gb GB", 40), color=:green)
        printstyled(lpad("($count files)", 20), "\n", color=:green)
    end

    printstyled("folder sizes on Julia...\n", color=:red)

    cwd = pwd()
    dirs = readdir(cwd)

    rs, rcnt = calculate_folder_size(cwd)

    print_folder_size(cwd,rs,rcnt)

    n = repeat(" - -", 10)
    printstyled(n * " subfolders of " * basename(cwd) * n, "\n", color=:yellow)

    for dir in dirs
        if isdir(dir)
            size, count = calculate_folder_size(joinpath(cwd, dir))
            print_folder_size(dir, size, count)
        end
    end
end

"""
Read wasim ts with DelimitedFiles.readdlm, skipto line 3 
no header column
"""
