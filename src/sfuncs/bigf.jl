# 
function bigf()

    printstyled("5 biggest folders...\n", color=:red)

    cwd = pwd()
    dirs = filter(isdir, readdir(; sort=false))
    
    rs, rcnt = [], []

    for i in dirs
        ts, tcnt = calculate_folder_size(i)
        push!(rs, round(ts / 1024^2, digits=3))  # Convert size to MB and round to 3 digits
        push!(rcnt, tcnt)
    end

    # Check if the matrix is empty
    if isempty(dirs)
        println("No directories found.")
        return
    end
    
    # Create an array of directories, sizes in MB, and file counts
    data = [(dir, rs[i], rcnt[i]) for (i, dir) in enumerate(dirs)]

    # Sort the array by size (second element) in descending order
    sorted_data = sort(data, by = x -> x[2], rev = true)

    if isempty(sorted_data)
        println("No directories with size information found.")
        return
    end

    # Extract the largest directories and display them in green
    largest_dirs = sorted_data[1:5]
    
    # Create a formatted string for printing
    dirs_to_display = join([string(dir[1], "\t\t Size: ", dir[2], " MB \t Files: ", dir[3], "\t") for dir in largest_dirs], "\n")
    #dirs_to_display = join([string(dir[1], " (Size: ", dir[2], " MB, Files: ", dir[3], ")") for dir in largest_dirs], "\n")

    printstyled("$dirs_to_display\n--> ", color=:green)
end

"""
import InteractiveUtils.clipboard
wslpath()|>clipboard
cb(wslpath())
"""
