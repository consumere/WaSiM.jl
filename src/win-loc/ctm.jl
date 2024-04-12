
# Define the root directory for the search
root_dir = "."  # Change this to the desired root directory

if length(ARGS) == 0
	println("need args! <file>...")
    exit()
end

fm=ARGS[1];

# Find all files containing "fm" and not "*intern*"
file_paths = filter(f -> (contains(f,fm) && !occursin("intern",f)), readdir(root_dir))

function grep_with_context(pattern, filename, context)
    lines = readlines(filename)   
    for (i, line) in enumerate(lines)
        if occursin(pattern, line)
                #println("$(filename):$(i): $(line)")
                for j in max(1, i - context):min(length(lines), i + context)
                    #printstyled("$(filename):$j: $(lines[j])\n",color=:red)
                    printstyled("$(filename):$j:",color=:red)
                    printstyled("$(lines[j])\n",color=:green)
                end
                println("-" ^ 80)
        end
    end
end
    
#grep_with_context("routing_model", file_paths[1], 2)
function grep_files(pattern, file_paths, context)
    for file_path in file_paths
        grep_with_context(pattern, file_path, context)
    end
end


grep_files("routing_model", file_paths, 2)