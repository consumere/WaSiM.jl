# 
function kgewrite(;output_file="kge-table.txt")
    path = pwd()
    files = glob(r"_output.txt|_outputjl") #non-recurse
    
    # Create an empty vector to store the matched lines
    matched_lines = Vector{String}()
    
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file, '\t', String)
        match = Grep.grep(r"KGE", output)
        if !isempty(match)
            fn = first(split(file, "_qout"))
            for line in sort(match, by = x -> parse(Float64, split(x)[end]); rev = true)
                line = strip(line)  # remove leading and trailing whitespace
                line = join(split(line), " ")  # remove inner whitespaces
                push!(matched_lines, "$fn: $line")  # collect the matched line
            end
        end
    end
    
    #output_file = "matched_results.csv"
    writedlm(output_file, matched_lines, '\t')
    println("Matched results saved to $output_file")
end

