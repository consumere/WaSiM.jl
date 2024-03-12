# 
function kgedf()
        path = pwd()
        files = glob(r"_output.txt|_outputjl") #non-recurse
        results = DataFrame(File = String[], KGE = Float64[])  # Create an empty DataFrame to store the results
        
        for file in filter(file -> endswith(file, "_output.txt"), files)
            output = DelimitedFiles.readdlm(file, '\t', String)
            match = Grep.grep(r"KGE", output)
            
            if !isempty(match)
                fn = first(split(file, "_qout"))
                parsed_lines = Float64[]
                
                for line in match
                    line_parts = split(line)
                    kge_value = parse(Float64, line_parts[end])
                    push!(parsed_lines, kge_value)
                end
                
                sort!(parsed_lines, rev = true)  # Sort the parsed KGE values in descending order
                
                for kge_value in parsed_lines
                    push!(results, (File = fn, KGE = kge_value))  # Add the result to the DataFrame
                end
            end
        end
        
        sort!(results, :KGE, rev = true)  # Sort the DataFrame by the 'KGE' column in descending order
        return results
    end

    