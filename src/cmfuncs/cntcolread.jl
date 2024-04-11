# 
function cntcolread(x::Vector{Any})
        # Get a list of all files in the directory
        files = x
    
        files = filter(inF->isfile(inF),files)
                            #if isfile(inF)
        file_columns = []
        
        for file in files
            # Open the file for reading
            open(file, "r") do io
                # Read the first line of the file
                line = readline(io)
                # Split the line on tabs
                columns = split(line, '\t')
                # Count the number of columns
                num_columns = length(columns)
                push!(file_columns, (file, num_columns))
            end
        end
        
        # Sort the file_columns array based on the number of columns in descending order
        sorted_files = sort(file_columns, by = x -> x[2], rev = true)
        
        for (file, num_columns) in sorted_files
            printstyled(
                rpad("File: $file",45),
            lpad(" | Columns: $num_columns\n",10),color=:green,bold=true)
        end
        return file_columns
    end

    